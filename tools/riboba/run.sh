#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
cd "$PROJECT_ROOT"

DATASET_ROOT="${DATASET_ROOT:-${TABRIBO_DATA_ROOT:-${PROJECT_ROOT}/data}/datasets/p_site_offset_prediction/landmark_riboseq}"
RESULT_ROOT="${RESULT_ROOT:-${PROJECT_ROOT}/results/p_site_offset_prediction/landmark_riboseq}"
NATIVE_OUTPUT_ROOT="${NATIVE_OUTPUT_ROOT:-${RESULT_ROOT}/native_outputs}"
OUTDIR="${OUTDIR:-${NATIVE_OUTPUT_ROOT}/riboba}"
SAMPLE_TABLE="${SAMPLE_TABLE:-${DATASET_ROOT}/samples.tsv}"
TEST="${TEST:-${DATASET_ROOT}/splits/within_sample/test.parquet}"
GENOME_BAM_DIR="${GENOME_BAM_DIR:-${DATASET_ROOT}/alignments/genome_bams}"
ANNOT_DIR="${ANNOT_DIR:-${DATASET_ROOT}/views/annotations}"
TRANSCRIPT_OUTDIR="${TRANSCRIPT_OUTDIR:-${DATASET_ROOT}/views/external_transcripts}"
PREDICTIONS="${PREDICTIONS:-${RESULT_ROOT}/within_sample/riboba/candidate_predictions.parquet}"
READ_PREDICTIONS="${READ_PREDICTIONS:-${RESULT_ROOT}/within_sample/riboba/read_predictions.parquet}"
SPLIT="${SPLIT:-p_site_offset_prediction_landmark_riboseq_within_sample}"
RUN_PIPELINE="${RUN_PIPELINE:-1}"
FORCE="${FORCE:-1}"

RIBOBA_REFS="${RIBOBA_REFS:-dmel scer}"
RIBOBA_FASTQ_MAX_READS="${RIBOBA_FASTQ_MAX_READS:-250000}"
RIBOBA_INPUT_SEED="${RIBOBA_INPUT_SEED:-1729}"
RIBOBA_THREADS="${RIBOBA_THREADS:-4}"
RIBOBA_RETRY_NO_ADD5="${RIBOBA_RETRY_NO_ADD5:-1}"
RIBOBA_FORCE_NO_ADD5="${RIBOBA_FORCE_NO_ADD5:-0}"
RIBOBA_OFFSET_MIN_GROUP_COUNT="${RIBOBA_OFFSET_MIN_GROUP_COUNT:-25}"
RIBOBA_TRUNC_MARGIN="${RIBOBA_TRUNC_MARGIN:-16}"
RIBOBA_SRC="${RIBOBA_SRC:-${PROJECT_ROOT}/tools/riboba/src}"
RIBOBA_PREPARE_SCRIPT="${RIBOBA_PREPARE_SCRIPT:-${RIBOBA_SRC}/R/prepare_bam.sh}"

if [[ -f "${PROJECT_ROOT}/scripts/local_env.sh" ]]; then
  source "${PROJECT_ROOT}/scripts/local_env.sh"
fi

R_LIBS_USER="${R_LIBS_USER:-${TABRIBO_EXTERNAL_RLIB_ROOT:-${PROJECT_ROOT}/external_tools/Rlib}}"
export R_LIBS_USER

LOG_DIR="${PROJECT_ROOT}/logs/baselines"
STATUS="${OUTDIR}/riboba_status.tsv"
PIPELINE_STATUS="${OUTDIR}/riboba_pipeline_status.tsv"
mkdir -p "$OUTDIR" "$LOG_DIR" "$(dirname "$PREDICTIONS")" "$(dirname "$READ_PREDICTIONS")"
printf 'sample_id\tref_key\tstatus\tmessage\n' > "$STATUS"
printf 'stage\tstatus\tmessage\n' > "$PIPELINE_STATUS"

read -r -a ref_keys <<< "$RIBOBA_REFS"
if [[ ${#ref_keys[@]} -eq 0 ]]; then
  printf 'overall\tfailed\tRIBOBA_REFS is empty\n' >> "$PIPELINE_STATUS"
  exit 1
fi

for ref_key in "${ref_keys[@]}"; do
  case "$ref_key" in
    dmel|scer) ;;
    *)
      printf 'overall\tfailed\tunsupported reference key: %s\n' "$ref_key" >> "$PIPELINE_STATUS"
      exit 1
      ;;
  esac
done

annotations_missing=0
for ref_key in "${ref_keys[@]}"; do
  if [[ ! -s "${ANNOT_DIR}/${ref_key}.riboba_flanked.gtf" ]]; then
    annotations_missing=1
  fi
done
if [[ "$annotations_missing" == "1" ]]; then
  python "${PROJECT_ROOT}/tools/shared/prepare_annotations.py" \
    --outdir "$ANNOT_DIR" \
    --transcript-outdir "$TRANSCRIPT_OUTDIR" \
    --refs "${ref_keys[@]}"
fi

if [[ "$FORCE" == "1" || ! -s "${OUTDIR}/riboba_input_status.tsv" ]]; then
  if ! python "${PROJECT_ROOT}/tools/riboba/prepare_inputs.py" \
      --sample-table "$SAMPLE_TABLE" \
      --test "$TEST" \
      --genome-bam-dir "$GENOME_BAM_DIR" \
      --outdir "$OUTDIR" \
      --max-background-reads "$RIBOBA_FASTQ_MAX_READS" \
      --seed "$RIBOBA_INPUT_SEED" \
      --refs "${ref_keys[@]}" \
      > "${LOG_DIR}/riboba_prepare_inputs.stdout.log" \
      2> "${LOG_DIR}/riboba_prepare_inputs.stderr.log"; then
    printf 'overall\tfailed\tinput preparation failed; see %s\n' \
      "${LOG_DIR}/riboba_prepare_inputs.stderr.log" >> "$PIPELINE_STATUS"
    exit 1
  fi
fi

mark_reference_samples() {
  local rnase_table="$1"
  local ref_key="$2"
  local status="$3"
  local message="$4"
  while IFS=$'\t' read -r sample_id _rest; do
    [[ "$sample_id" == "sample_id" || -z "$sample_id" ]] && continue
    printf '%s\t%s\t%s\t%s\n' "$sample_id" "$ref_key" "$status" "$message" >> "$STATUS"
  done < "$rnase_table"
}

successful_profiles=0
attempted_profiles=0

for ref_key in "${ref_keys[@]}"; do
  ref_root="${OUTDIR}/refs/${ref_key}"
  rnase_table="${ref_root}/rnase.tsv"
  fastq_dir="${OUTDIR}/fastq/${ref_key}"
  annotation_outdir="${ref_root}/annotation"
  prepare_outdir="${ref_root}/prepare_bam"
  profile_status="${ref_root}/profile_status.tsv"

  if [[ ! -s "$rnase_table" ]]; then
    continue
  fi
  ref_samples="$(awk -F '\t' 'NR > 1 {count += 1} END {print count + 0}' "$rnase_table")"
  attempted_profiles=$((attempted_profiles + ref_samples))

  case "$ref_key" in
    dmel)
      scientific_name="Drosophila.melanogaster"
      annotation_version="BDGP6.32"
      ;;
    scer)
      scientific_name="Saccharomyces.cerevisiae"
      annotation_version="R64-1-1"
      ;;
  esac

  genome_file="${DATASET_ROOT}/references/${ref_key}/genome.fa"
  gff_file="${ANNOT_DIR}/${ref_key}.riboba_flanked.gtf"
  rrna_file="${DATASET_ROOT}/references/${ref_key}/ncrna.fa"
  tx_info="${annotation_outdir}/transcript/tx_information.RData"
  cds_fasta="${annotation_outdir}/transcript/longest_txs.fa"
  lncrna_fasta="${annotation_outdir}/transcript/lncRNA_txs.fa"
  annotation_rebuilt=0
  annotation_complete=1
  for annotation_file in "$tx_info" "$cds_fasta" "$lncrna_fasta"; do
    if [[ ! -s "$annotation_file" ]]; then
      annotation_complete=0
    fi
  done

  if [[ "$FORCE" == "1" || "$annotation_complete" != "1" ]]; then
    clear_annotation=1
    if RIBOBA_ANNOT_DIR="$annotation_outdir" \
        RIBOBA_CLEAR_ANNOT_DIR="$clear_annotation" \
        RIBOBA_SCIENTIFIC_NAME="$scientific_name" \
        RIBOBA_ANNOTATION_VERSION="$annotation_version" \
        RIBOBA_GENOME_FILE="$genome_file" \
        RIBOBA_GFF_FILE="$gff_file" \
        Rscript "${PROJECT_ROOT}/tools/riboba/run_riboba_annotation.R" \
        > "${LOG_DIR}/riboba_annotation_${ref_key}.stdout.log" \
        2> "${LOG_DIR}/riboba_annotation_${ref_key}.stderr.log"; then
      annotation_rebuilt=1
    else
      mark_reference_samples \
        "$rnase_table" "$ref_key" "failed" \
        "native annotation build failed; see riboba_annotation_${ref_key}.stderr.log"
      continue
    fi
  fi

  prepare_needed="$FORCE"
  if [[ "$annotation_rebuilt" == "1" ]]; then
    prepare_needed=1
  fi
  while IFS=$'\t' read -r sample_id _rest; do
    [[ "$sample_id" == "sample_id" || -z "$sample_id" ]] && continue
    if [[ ! -s "${prepare_outdir}/bam/${sample_id}_cds_txsorted.bam" ]]; then
      prepare_needed=1
    fi
  done < "$rnase_table"

  if [[ "$prepare_needed" == "1" ]]; then
    if ! PATH="${PROJECT_ROOT}/.conda-envs/tabribo-local/bin:${PATH}" \
        bash "$RIBOBA_PREPARE_SCRIPT" \
          --workdir "$prepare_outdir" \
          --fastq-dir "$fastq_dir" \
          --annotate-dir "$annotation_outdir" \
          --genome-file "$genome_file" \
          --species "${scientific_name/./ }" \
          --rrna-fa "$rrna_file" \
          --cds-fa "${annotation_outdir}/transcript/longest_txs.fa" \
          --lncrna-fa "${annotation_outdir}/transcript/lncRNA_txs.fa" \
          --threads "$RIBOBA_THREADS" \
          > "${LOG_DIR}/riboba_prepare_bam_${ref_key}.stdout.log" \
          2> "${LOG_DIR}/riboba_prepare_bam_${ref_key}.stderr.log"; then
      mark_reference_samples \
        "$rnase_table" "$ref_key" "failed" \
        "native BAM preparation failed; see riboba_prepare_bam_${ref_key}.stderr.log"
      continue
    fi
  fi

  if [[ "$RUN_PIPELINE" != "1" ]]; then
    mark_reference_samples "$rnase_table" "$ref_key" "prepared" "native BAMs prepared; model run skipped"
    continue
  fi

  printf 'sample_id\tstatus\tmessage\n' > "$profile_status"
  mkdir -p "${ref_root}/results" "${ref_root}/profile_status"
  while IFS=$'\t' read -r sample_id rnase _rest; do
    [[ "$sample_id" == "sample_id" || -z "$sample_id" ]] && continue
    psite_output="${ref_root}/downstream/${sample_id}/${sample_id}.riboba_psite_profile.tsv"
    if [[ "$FORCE" != "1" && -s "$psite_output" ]]; then
      printf '%s\tok\texisting native RiboBA P-site profile; RNase=%s\n' \
        "$sample_id" "$rnase" >> "$profile_status"
      continue
    fi

    sample_profile_status="${ref_root}/profile_status/${sample_id}.tsv"
    printf 'sample_id\tstatus\tmessage\n' > "$sample_profile_status"
    pipeline_rc=0
    RIBOBA_SAMPLE="$sample_id" \
      RIBOBA_TX_INFO="$tx_info" \
      RIBOBA_BAM_DIR="${prepare_outdir}/bam" \
      RIBOBA_RNASE_TABLE="$rnase_table" \
      RIBOBA_RESULT="${ref_root}/results/${sample_id}.rds" \
      RIBOBA_TX_MAP="${ref_root}/riboba_tx_map.tsv" \
      RIBOBA_DOWNSTREAM_DIR="${ref_root}/downstream" \
      RIBOBA_PROFILE_STATUS="$sample_profile_status" \
      RIBOBA_RETRY_NO_ADD5="$RIBOBA_RETRY_NO_ADD5" \
      RIBOBA_FORCE_NO_ADD5="$RIBOBA_FORCE_NO_ADD5" \
      RIBOBA_EXPORT_DOWNSTREAM=1 \
      RIBOBA_OFFSET_MIN_GROUP_COUNT="$RIBOBA_OFFSET_MIN_GROUP_COUNT" \
      RIBOBA_ALLOW_SYNTHETIC_OFFSETS=0 \
      RIBOBA_TRUNC_MARGIN="$RIBOBA_TRUNC_MARGIN" \
      RIBOBA_ALLOW_PROFILE_FALLBACK=0 \
      RIBOBA_PROFILE_ONLY=1 \
      Rscript "${PROJECT_ROOT}/tools/riboba/run_riboba_pipeline.R" \
      > "${LOG_DIR}/riboba_pipeline_${ref_key}_${sample_id}.stdout.log" \
      2> "${LOG_DIR}/riboba_pipeline_${ref_key}_${sample_id}.stderr.log" || pipeline_rc=$?

    if [[ -s "$sample_profile_status" && "$(wc -l < "$sample_profile_status")" -gt 1 ]]; then
      awk 'NR > 1' "$sample_profile_status" >> "$profile_status"
    else
      printf '%s\tfailed\tnative profile process exited %s before status output; see riboba_pipeline_%s_%s.stderr.log\n' \
        "$sample_id" "$pipeline_rc" "$ref_key" "$sample_id" >> "$profile_status"
    fi
  done < "$rnase_table"

  awk -F '\t' -v OFS='\t' -v ref="$ref_key" \
    'NR > 1 {print $1, ref, $2, $3}' "$profile_status" >> "$STATUS"
  ref_successes="$(awk -F '\t' 'NR > 1 && $2 == "ok" {count += 1} END {print count + 0}' "$profile_status")"
  successful_profiles=$((successful_profiles + ref_successes))
done

if [[ "$RUN_PIPELINE" != "1" ]]; then
  printf 'overall\tok\tinput and BAM preparation completed; RUN_PIPELINE is not 1\n' >> "$PIPELINE_STATUS"
  exit 0
fi

if [[ "$successful_profiles" == "0" ]]; then
  printf 'overall\tfailed\tno sample produced a native RiboBA P-site profile (%s attempted)\n' \
    "$attempted_profiles" >> "$PIPELINE_STATUS"
  exit 1
fi

python "${PROJECT_ROOT}/tools/shared/convert_profile_predictions.py" \
  --method riboba \
  --test "$TEST" \
  --sample-table "$SAMPLE_TABLE" \
  --baseline-root "$OUTDIR" \
  --annotation-dir "$ANNOT_DIR" \
  --output "$PREDICTIONS"

prediction_rows="$(
  python -c \
    'import pyarrow.parquet as pq, sys; print(pq.ParquetFile(sys.argv[1]).metadata.num_rows)' \
    "$PREDICTIONS"
)"
if [[ "$prediction_rows" -gt 0 ]]; then
  python "${PROJECT_ROOT}/tabriboicl.py" evaluate \
    --test "$TEST" \
    --predictions "$PREDICTIONS" \
    --read-output "$READ_PREDICTIONS" \
    --model riboba \
    --split "$SPLIT"
else
  python -c \
    'import pandas as pd, sys; pd.DataFrame(columns=["read_id", "pred_offset", "pred_frame", "p_pred", "true_offset", "true_frame", "q_true"]).to_parquet(sys.argv[1], index=False)' \
    "$READ_PREDICTIONS"
fi

printf 'overall\tok\tnative RiboBA P-site profiles generated for %s/%s samples; emitted %s candidate rows\n' \
  "$successful_profiles" "$attempted_profiles" "$prediction_rows" >> "$PIPELINE_STATUS"
