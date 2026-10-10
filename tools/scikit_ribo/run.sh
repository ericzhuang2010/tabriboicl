#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
DATASET_ROOT="${DATASET_ROOT:-${TABRIBO_DATA_ROOT:-${PROJECT_ROOT}/data}/datasets/p_site_offset_prediction/landmark_riboseq}"
RESULT_ROOT="${RESULT_ROOT:-${PROJECT_ROOT}/results/p_site_offset_prediction/landmark_riboseq}"
NATIVE_OUTPUT_ROOT="${NATIVE_OUTPUT_ROOT:-${RESULT_ROOT}/native_outputs}"

if [[ -f "${PROJECT_ROOT}/scripts/local_env.sh" ]]; then
  source "${PROJECT_ROOT}/scripts/local_env.sh"
fi

SAMPLE_TABLE="${SAMPLE_TABLE:-${DATASET_ROOT}/samples.tsv}"
BAM_DIR="${BAM_DIR:-${DATASET_ROOT}/alignments/genome_bams}"
OUTDIR="${OUTDIR:-${NATIVE_OUTPUT_ROOT}/scikit_ribo}"
INDEX_ROOT="${INDEX_ROOT:-${DATASET_ROOT}/views/scikit_ribo_indexes}"
SCIKIT_RIBO_SRC="${SCIKIT_RIBO_SRC:-${PROJECT_ROOT}/tools/scikit_ribo/src}"
PREDICTIONS="${PREDICTIONS:-${RESULT_ROOT}/within_sample/scikit_ribo/candidate_predictions.parquet}"
READ_PREDICTIONS="${READ_PREDICTIONS:-${RESULT_ROOT}/within_sample/scikit_ribo/read_predictions.parquet}"
TEST="${TEST:-${DATASET_ROOT}/splits/within_sample/test.parquet}"
SPLIT="${SPLIT:-p_site_offset_prediction_landmark_riboseq_within_sample}"
THREADS="${THREADS:-4}"
FORCE="${FORCE:-1}"
SCIKIT_RIBO_MAPQ="${SCIKIT_RIBO_MAPQ:-20}"
SCIKIT_RIBO_MIN_READ_LEN="${SCIKIT_RIBO_MIN_READ_LEN:-10}"
SCIKIT_RIBO_MAX_READ_LEN="${SCIKIT_RIBO_MAX_READ_LEN:-35}"

if [[ ! -d "$SCIKIT_RIBO_SRC/scikit_ribo" ]]; then
  printf 'Missing scikit-ribo source at %s. Run scripts/fetch_tool_sources.sh --yes first.\n' "$SCIKIT_RIBO_SRC" >&2
  exit 1
fi

mkdir -p "$OUTDIR" "$INDEX_ROOT" "logs/baselines" "$(dirname "$PREDICTIONS")"
mkdir -p "$(dirname "$READ_PREDICTIONS")"
STATUS="${OUTDIR}/scikit_ribo_status.tsv"
printf 'sample_id\tstatus\tmessage\n' > "$STATUS"
force_args=()
if [[ "$FORCE" == "1" ]]; then
  force_args=(--force)
fi

for ref_key in dmel scer; do
  ref_index="${INDEX_ROOT}/${ref_key}"
  if [[ "$FORCE" != "1" && -s "${ref_index}/${ref_key}.codons.df" ]]; then
    printf '[%s] skipping scikit-ribo indexes for %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$ref_key"
  else
    python "${PROJECT_ROOT}/tools/scikit_ribo/asite.py" \
      --src-dir "$SCIKIT_RIBO_SRC" \
      build \
      --gtf "${DATASET_ROOT}/references/${ref_key}/annotation.gtf" \
      --fasta "${DATASET_ROOT}/references/${ref_key}/genome.fa" \
      --prefix "$ref_key" \
      --outdir "$ref_index" \
      "${force_args[@]}" \
      > "logs/baselines/${ref_key}.scikit_ribo_build.stdout.log" \
      2> "logs/baselines/${ref_key}.scikit_ribo_build.stderr.log"
  fi
done

tail -n +2 "$SAMPLE_TABLE" | while IFS=$'\t' read -r sample organism ref_key study_key sra_gb rest; do
  [[ -z "$sample" ]] && continue
  bam="${BAM_DIR}/${sample}.Aligned.sortedByCoord.out.bam"
  sample_out="${OUTDIR}/${sample}"
  mkdir -p "$sample_out"

  if [[ ! -s "$bam" ]]; then
    printf '%s\tfailed\tmissing genome BAM\n' "$sample" >> "$STATUS"
    continue
  fi

  if python "${PROJECT_ROOT}/tools/scikit_ribo/asite.py" \
      --src-dir "$SCIKIT_RIBO_SRC" \
      run \
      --bam "$bam" \
      --index-dir "${INDEX_ROOT}/${ref_key}" \
      --prefix "$ref_key" \
      --outdir "$sample_out" \
      --mapq "$SCIKIT_RIBO_MAPQ" \
      --min-read-len "$SCIKIT_RIBO_MIN_READ_LEN" \
      --max-read-len "$SCIKIT_RIBO_MAX_READ_LEN" \
      --threads "$THREADS" \
      --target-test "$TEST" \
      "${force_args[@]}" \
      > "logs/baselines/${sample}.scikit_ribo.stdout.log" \
      2> "logs/baselines/${sample}.scikit_ribo.stderr.log"; then
    printf '%s\tok\tprofile written to %s/riboseq_input.txt\n' "$sample" "$sample_out" >> "$STATUS"
  else
    printf '%s\tfailed\tscikit-ribo failed; see logs/baselines/%s.scikit_ribo.stderr.log\n' "$sample" "$sample" >> "$STATUS"
  fi
done

python "${PROJECT_ROOT}/tools/scikit_ribo/convert.py" \
  --test "$TEST" \
  --baseline-root "$OUTDIR" \
  --output "$PREDICTIONS"

python "${PROJECT_ROOT}/tabriboicl.py" evaluate \
  --test "$TEST" \
  --predictions "$PREDICTIONS" \
  --read-output "$READ_PREDICTIONS" \
  --model scikit_ribo \
  --split "$SPLIT"
