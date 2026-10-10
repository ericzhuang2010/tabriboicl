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
BAM_DIR="${BAM_DIR:-${DATASET_ROOT}/views/transcriptome_bowtie2}"
ANNOT_DIR="${ANNOT_DIR:-${DATASET_ROOT}/views/annotations}"
TRANSCRIPT_FASTA_DIR="${TRANSCRIPT_FASTA_DIR:-${DATASET_ROOT}/views/external_transcripts}"
OUTDIR="${OUTDIR:-${NATIVE_OUTPUT_ROOT}/psite}"
PREDICTIONS="${PREDICTIONS:-${RESULT_ROOT}/within_sample/psite/candidate_predictions.parquet}"
READ_PREDICTIONS="${READ_PREDICTIONS:-${RESULT_ROOT}/within_sample/psite/read_predictions.parquet}"
TEST="${TEST:-${DATASET_ROOT}/splits/within_sample/test.parquet}"
SPLIT="${SPLIT:-p_site_offset_prediction_landmark_riboseq_within_sample}"
THREADS="${THREADS:-4}"
PSITE_FALLBACK_ON_TRAIN_FAILURE="${PSITE_FALLBACK_ON_TRAIN_FAILURE:-0}"
PSITE_FALLBACK_OFFSET="${PSITE_FALLBACK_OFFSET:-12}"
PSITE_FALLBACK_RLEN_MIN="${PSITE_FALLBACK_RLEN_MIN:-25}"
PSITE_FALLBACK_RLEN_MAX="${PSITE_FALLBACK_RLEN_MAX:-40}"

PSITE_TRAIN_MAX_ALIGNMENTS="${PSITE_TRAIN_MAX_ALIGNMENTS:-2000000}"
TABRIBO_PSITE_RANDOM_SEED="${TABRIBO_PSITE_RANDOM_SEED:-20260517}"
export TABRIBO_PSITE_RANDOM_SEED

detect_psite_python() {
  if [[ -x "${TABRIBO_EXTERNAL_ENV_ROOT:-${PROJECT_ROOT}/external_tools/envs}/psite/bin/python" ]]; then
    printf '%s\n' "${TABRIBO_EXTERNAL_ENV_ROOT:-${PROJECT_ROOT}/external_tools/envs}/psite/bin/python"
    return 0
  fi
  command -v python >/dev/null 2>&1 || return 1
  command -v python
}

detect_psite_site_packages() {
  local matches=("${PROJECT_ROOT}"/external_tools/envs/psite/lib/python*/site-packages)
  if [[ ${#matches[@]} -gt 0 && -d "${matches[0]}" ]]; then
    printf '%s\n' "${matches[0]}"
  fi
}

if [[ -z "${PSITE_PYTHON:-}" ]]; then
  PSITE_PYTHON="$(detect_psite_python || true)"
fi
if [[ -z "${PSITE_SITE_PACKAGES:-}" ]]; then
  PSITE_SITE_PACKAGES="$(detect_psite_site_packages || true)"
fi

if [[ -z "${PSITE_PYTHON:-}" ]]; then
  printf 'PSite python is unavailable. Set PSITE_PYTHON or create external_tools/envs/psite.\n' >&2
  exit 1
fi

mkdir -p "$OUTDIR/models" "$OUTDIR/tagged_bams" "${PROJECT_ROOT}/logs/baselines" "$(dirname "$PREDICTIONS")"
mkdir -p "$(dirname "$READ_PREDICTIONS")"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${PROJECT_ROOT}/.matplotlib}"
mkdir -p "$MPLCONFIGDIR"

STATUS="${OUTDIR}/psite_status.tsv"
printf 'sample_id\tstatus\tmessage\n' > "$STATUS"

run_psite() {
  if [[ -n "${PSITE_SITE_PACKAGES:-}" ]]; then
    PYTHONPATH="${PSITE_SITE_PACKAGES}${PYTHONPATH:+:${PYTHONPATH}}" \
      "$PSITE_PYTHON" "${PROJECT_ROOT}/tools/psite/cli.py" "$@"
  else
    "$PSITE_PYTHON" "${PROJECT_ROOT}/tools/psite/cli.py" "$@"
  fi
}

if ! run_psite --help >/dev/null 2>&1; then
  printf 'PSite is not importable. Set PSITE_PYTHON/PSITE_SITE_PACKAGES or run scripts/setup_linux_p_site_offset_prediction.sh.\n' >&2
  exit 1
fi

tail -n +2 "$SAMPLE_TABLE" | while IFS=$'\t' read -r sample organism ref_key study_key sra_gb rest; do
  [[ -z "$sample" ]] && continue
  bam="${BAM_DIR}/${sample}.transcriptome.bam"
  txinfo="${ANNOT_DIR}/${ref_key}.psite_txinfo.tsv"
  transcript_fasta="${TRANSCRIPT_FASTA_DIR}/${ref_key}.fa"
  if [[ ! -s "$transcript_fasta" ]]; then
    transcript_fasta="${DATASET_ROOT}/references/${ref_key}/cdna.fa"
  fi
  model_prefix="${OUTDIR}/models/${sample}"
  tagged_bam="${OUTDIR}/tagged_bams/${sample}.psite.tagged.bam"

  if [[ ! -s "$bam" ]]; then
    printf '%s\tfailed\tmissing transcriptome BAM\n' "$sample" >> "$STATUS"
    continue
  fi

  rm -f "${model_prefix}.dam.pickle" "${model_prefix}.gbt.pickle" "$tagged_bam"

  train_fraction="1.0"
  if [[ "$PSITE_TRAIN_MAX_ALIGNMENTS" -gt 0 ]]; then
    mapped_alignments="$(samtools view -c -F 4 "$bam")"
    if [[ "$mapped_alignments" -gt "$PSITE_TRAIN_MAX_ALIGNMENTS" ]]; then
      train_fraction="$(awk -v cap="$PSITE_TRAIN_MAX_ALIGNMENTS" -v total="$mapped_alignments" 'BEGIN { printf "%.12f", cap / total }')"
    fi
  fi

  train_status="ok"
  if ! run_psite train \
      "$transcript_fasta" \
      "$bam" \
      "$model_prefix" \
      "$txinfo" \
      -f "$train_fraction" \
      > "${PROJECT_ROOT}/logs/baselines/${sample}.psite_train.stdout.log" \
      2> "${PROJECT_ROOT}/logs/baselines/${sample}.psite_train.stderr.log"; then
    train_status="failed_after_model"
    if [[ ! -s "${model_prefix}.gbt.pickle" ]]; then
      if [[ "$PSITE_FALLBACK_ON_TRAIN_FAILURE" == "1" ]]; then
        if run_psite setp \
            "$bam" \
            "$tagged_bam" \
            -l "$PSITE_FALLBACK_RLEN_MIN" \
            -u "$PSITE_FALLBACK_RLEN_MAX" \
            -n "$PSITE_FALLBACK_OFFSET" \
            > "${PROJECT_ROOT}/logs/baselines/${sample}.psite_setp.stdout.log" \
            2> "${PROJECT_ROOT}/logs/baselines/${sample}.psite_setp.stderr.log"; then
          printf '%s\tok\tfallback setp fixed offset %s after training failed; see logs/baselines/%s.psite_train.stderr.log\n' "$sample" "$PSITE_FALLBACK_OFFSET" "$sample" >> "$STATUS"
          continue
        fi
        printf '%s\tfailed\ttraining failed and fallback setp failed; see logs/baselines/%s.psite_train.stderr.log and %s.psite_setp.stderr.log\n' "$sample" "$sample" "$sample" >> "$STATUS"
        continue
      fi
      printf '%s\tfailed\ttraining failed; see logs/baselines/%s.psite_train.stderr.log and %s.log\n' "$sample" "$sample" "$model_prefix" >> "$STATUS"
      continue
    fi
  fi

  if [[ -s "${model_prefix}.gbt.pickle" ]]; then
    if run_psite predict \
        "$transcript_fasta" \
        "$bam" \
        "${model_prefix}.gbt.pickle" \
        "$tagged_bam" \
        > "${PROJECT_ROOT}/logs/baselines/${sample}.psite_predict.stdout.log" \
        2> "${PROJECT_ROOT}/logs/baselines/${sample}.psite_predict.stderr.log"; then
      if [[ "$train_status" == "ok" ]]; then
        printf '%s\tok\ttrained and predicted\n' "$sample" >> "$STATUS"
      else
        printf '%s\tok\ttrained and predicted; train exited after writing GBT model\n' "$sample" >> "$STATUS"
      fi
    else
      rm -f "$tagged_bam"
      printf '%s\tfailed\tprediction failed; see logs/baselines/%s.psite_predict.stderr.log\n' "$sample" "$sample" >> "$STATUS"
    fi
  fi
done

python "${PROJECT_ROOT}/tools/psite/convert.py" \
  --sample-table "$SAMPLE_TABLE" \
  --bam-dir "${OUTDIR}/tagged_bams" \
  --test "$TEST" \
  --output "$PREDICTIONS"

python "${PROJECT_ROOT}/tabriboicl.py" evaluate \
  --test "$TEST" \
  --predictions "$PREDICTIONS" \
  --read-output "$READ_PREDICTIONS" \
  --model psite \
  --split "$SPLIT"
