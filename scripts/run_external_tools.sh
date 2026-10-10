#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
DATASET_ROOT="${DATASET_ROOT:-${TABRIBO_DATA_ROOT:-${PROJECT_ROOT}/data}/datasets/p_site_offset_prediction/landmark_riboseq}"
RESULT_ROOT="${RESULT_ROOT:-${PROJECT_ROOT}/results/p_site_offset_prediction/landmark_riboseq}"
NATIVE_OUTPUT_ROOT="${NATIVE_OUTPUT_ROOT:-${RESULT_ROOT}/native_outputs}"
SAMPLE_TABLE="${SAMPLE_TABLE:-${DATASET_ROOT}/samples.tsv}"
THREADS="${THREADS:-4}"
TRANSCRIPTOME_MAX_READS="${TRANSCRIPTOME_MAX_READS:-all}"
EXTERNAL_TOOLS_STR="${EXTERNAL_TOOLS_STR:-psite ribowaltz ribodeblur riboba scikit_ribo asite_ip}"
FORCE="${FORCE:-1}"
EXTERNAL_STRICT="${EXTERNAL_STRICT:-0}"
RIBOBA_RUN_PIPELINE="${RIBOBA_RUN_PIPELINE:-1}"
RESULT_WITHIN_SAMPLE_DIR="${RESULT_ROOT}/within_sample"

if [[ -f "${PROJECT_ROOT}/scripts/local_env.sh" ]]; then
  source "${PROJECT_ROOT}/scripts/local_env.sh"
fi

mkdir -p "$NATIVE_OUTPUT_ROOT" "$RESULT_WITHIN_SAMPLE_DIR" "${PROJECT_ROOT}/logs/baselines"
for tool in psite ribowaltz ribodeblur riboba scikit_ribo asite_ip; do
  mkdir -p "$RESULT_WITHIN_SAMPLE_DIR/$tool"
  printf 'step,status,message\n' > "$RESULT_WITHIN_SAMPLE_DIR/$tool/run_status.csv"
done

external_tool_enabled() {
  local tool="$1"
  [[ " $EXTERNAL_TOOLS_STR " == *" $tool "* ]]
}

record_external_status() {
  local tool="$1"
  local status="$2"
  local message="$3"
  case "$tool" in
    psite|ribowaltz|ribodeblur|riboba|scikit_ribo|asite_ip)
      local status_dir="$RESULT_WITHIN_SAMPLE_DIR/$tool"
      local status_file="$status_dir/run_status.csv"
      mkdir -p "$status_dir"
      printf '%s,%s,%s\n' "$tool" "$status" "$message" >> "$status_file"
      ;;
  esac
}

run_external_step() {
  local tool="$1"
  local sentinel="$2"
  shift 2

  if [[ "$FORCE" != "1" && -s "$sentinel" ]]; then
    printf '[%s] skipping %s; found %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$tool" "$sentinel"
    record_external_status "$tool" "skipped" "found $sentinel"
    return 0
  fi

  printf '[%s] running %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$tool"
  if "$@"; then
    if [[ "$(basename "$sentinel")" == .* || "$sentinel" == *.ok ]]; then
      mkdir -p "$(dirname "$sentinel")"
      touch "$sentinel"
    fi
    record_external_status "$tool" "ok" "completed"
    return 0
  fi

  record_external_status "$tool" "failed" "command failed; inspect logs/baselines and tool-specific status files"
  if [[ "$EXTERNAL_STRICT" == "1" ]]; then
    return 1
  fi
  return 0
}

printf '[%s] external tool preparation\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
run_external_step \
  "annotations" \
  "${DATASET_ROOT}/views/annotations/.prepared.ok" \
  python "${PROJECT_ROOT}/tools/shared/prepare_annotations.py" \
    --outdir "${DATASET_ROOT}/views/annotations" \
    --transcript-outdir "${DATASET_ROOT}/views/external_transcripts"

if external_tool_enabled psite || external_tool_enabled ribowaltz || external_tool_enabled ribodeblur || external_tool_enabled asite_ip; then
  if [[ "$FORCE" == "1" ]]; then
    rm -f "${DATASET_ROOT}/views/transcriptome_bowtie2"/.prepared_*.ok
  fi
  run_external_step \
    "transcriptome_bowtie2" \
    "${DATASET_ROOT}/views/transcriptome_bowtie2/.prepared_${TRANSCRIPTOME_MAX_READS}.ok" \
    env DATASET_ROOT="$DATASET_ROOT" SAMPLE_TABLE="$SAMPLE_TABLE" OUTDIR="${DATASET_ROOT}/views/transcriptome_bowtie2" THREADS="$THREADS" MAX_READS="$TRANSCRIPTOME_MAX_READS" FORCE="$FORCE" \
      bash "${PROJECT_ROOT}/tools/shared/prepare_transcriptome_bams.sh"
fi

if external_tool_enabled psite; then
  run_external_step \
    "psite" \
    "$RESULT_WITHIN_SAMPLE_DIR/psite/read_predictions.parquet" \
    env DATASET_ROOT="$DATASET_ROOT" RESULT_ROOT="$RESULT_ROOT" NATIVE_OUTPUT_ROOT="$NATIVE_OUTPUT_ROOT" SAMPLE_TABLE="$SAMPLE_TABLE" \
      BAM_DIR="${DATASET_ROOT}/views/transcriptome_bowtie2" ANNOT_DIR="${DATASET_ROOT}/views/annotations" \
      OUTDIR="${NATIVE_OUTPUT_ROOT}/psite" PREDICTIONS="$RESULT_WITHIN_SAMPLE_DIR/psite/candidate_predictions.parquet" \
      READ_PREDICTIONS="$RESULT_WITHIN_SAMPLE_DIR/psite/read_predictions.parquet" THREADS="$THREADS" \
      bash "${PROJECT_ROOT}/tools/psite/run.sh"
else
  record_external_status "psite" "skipped" "not listed in EXTERNAL_TOOLS_STR"
fi

if external_tool_enabled ribowaltz; then
  run_external_step \
    "ribowaltz" \
    "$RESULT_WITHIN_SAMPLE_DIR/ribowaltz/read_predictions.parquet" \
    env DATASET_ROOT="$DATASET_ROOT" RESULT_ROOT="$RESULT_ROOT" NATIVE_OUTPUT_ROOT="$NATIVE_OUTPUT_ROOT" SAMPLE_TABLE="$SAMPLE_TABLE" \
      BAM_DIR="${NATIVE_OUTPUT_ROOT}/ribowaltz/bams" ANNOT_DIR="${DATASET_ROOT}/views/annotations" \
      OUTDIR="${NATIVE_OUTPUT_ROOT}/ribowaltz" PREDICTIONS="$RESULT_WITHIN_SAMPLE_DIR/ribowaltz/candidate_predictions.parquet" \
      READ_PREDICTIONS="$RESULT_WITHIN_SAMPLE_DIR/ribowaltz/read_predictions.parquet" \
      bash "${PROJECT_ROOT}/tools/ribowaltz/run.sh"
else
  record_external_status "ribowaltz" "skipped" "not listed in EXTERNAL_TOOLS_STR"
fi

if external_tool_enabled ribodeblur; then
  run_external_step \
    "ribodeblur" \
    "$RESULT_WITHIN_SAMPLE_DIR/ribodeblur/candidate_predictions.parquet" \
    env DATASET_ROOT="$DATASET_ROOT" RESULT_ROOT="$RESULT_ROOT" NATIVE_OUTPUT_ROOT="$NATIVE_OUTPUT_ROOT" SAMPLE_TABLE="$SAMPLE_TABLE" \
      BAM_DIR="${DATASET_ROOT}/views/transcriptome_bowtie2" ANNOT_DIR="${DATASET_ROOT}/views/annotations" \
      OUTDIR="${NATIVE_OUTPUT_ROOT}/ribodeblur" PREDICTIONS="$RESULT_WITHIN_SAMPLE_DIR/ribodeblur/candidate_predictions.parquet" \
      READ_PREDICTIONS="$RESULT_WITHIN_SAMPLE_DIR/ribodeblur/read_predictions.parquet" \
      TEST="${DATASET_ROOT}/splits/within_sample/test.parquet" SPLIT="p_site_offset_prediction_landmark_riboseq_within_sample" \
      bash "${PROJECT_ROOT}/tools/ribodeblur/run.sh"
else
  record_external_status "ribodeblur" "skipped" "not listed in EXTERNAL_TOOLS_STR"
fi

if external_tool_enabled riboba; then
  run_external_step \
    "riboba" \
    "$RESULT_WITHIN_SAMPLE_DIR/riboba/candidate_predictions.parquet" \
    env DATASET_ROOT="$DATASET_ROOT" RESULT_ROOT="$RESULT_ROOT" NATIVE_OUTPUT_ROOT="$NATIVE_OUTPUT_ROOT" OUTDIR="${NATIVE_OUTPUT_ROOT}/riboba" RUN_PIPELINE="$RIBOBA_RUN_PIPELINE" FORCE="$FORCE" \
      SAMPLE_TABLE="$SAMPLE_TABLE" ANNOT_DIR="${DATASET_ROOT}/views/annotations" RIBOBA_THREADS="$THREADS" \
      PREDICTIONS="$RESULT_WITHIN_SAMPLE_DIR/riboba/candidate_predictions.parquet" READ_PREDICTIONS="$RESULT_WITHIN_SAMPLE_DIR/riboba/read_predictions.parquet" \
      TEST="${DATASET_ROOT}/splits/within_sample/test.parquet" SPLIT="p_site_offset_prediction_landmark_riboseq_within_sample" \
      bash "${PROJECT_ROOT}/tools/riboba/run.sh"
else
  record_external_status "riboba" "skipped" "not listed in EXTERNAL_TOOLS_STR"
fi

if external_tool_enabled scikit_ribo; then
  run_external_step \
    "scikit_ribo" \
    "$RESULT_WITHIN_SAMPLE_DIR/scikit_ribo/candidate_predictions.parquet" \
    env DATASET_ROOT="$DATASET_ROOT" RESULT_ROOT="$RESULT_ROOT" NATIVE_OUTPUT_ROOT="$NATIVE_OUTPUT_ROOT" SAMPLE_TABLE="$SAMPLE_TABLE" \
      BAM_DIR="${DATASET_ROOT}/alignments/genome_bams" OUTDIR="${NATIVE_OUTPUT_ROOT}/scikit_ribo" INDEX_ROOT="${DATASET_ROOT}/views/scikit_ribo_indexes" \
      PREDICTIONS="$RESULT_WITHIN_SAMPLE_DIR/scikit_ribo/candidate_predictions.parquet" READ_PREDICTIONS="$RESULT_WITHIN_SAMPLE_DIR/scikit_ribo/read_predictions.parquet" \
      TEST="${DATASET_ROOT}/splits/within_sample/test.parquet" SPLIT="p_site_offset_prediction_landmark_riboseq_within_sample" \
      THREADS="$THREADS" FORCE="$FORCE" \
      bash "${PROJECT_ROOT}/tools/scikit_ribo/run.sh"
else
  record_external_status "scikit_ribo" "skipped" "not listed in EXTERNAL_TOOLS_STR"
fi

run_asite_ip_with_conversion() {
  DATASET_ROOT="$DATASET_ROOT" \
    RESULT_ROOT="$RESULT_ROOT" \
    NATIVE_OUTPUT_ROOT="$NATIVE_OUTPUT_ROOT" \
    SAMPLE_TABLE="$SAMPLE_TABLE" \
    BAM_DIR="${DATASET_ROOT}/views/transcriptome_bowtie2" \
    ANNOT_DIR="${DATASET_ROOT}/views/annotations" \
    REFERENCE_DIR="${DATASET_ROOT}/references" \
    OUTDIR="${NATIVE_OUTPUT_ROOT}/asite_ip" \
    FORCE="$FORCE" \
    bash "${PROJECT_ROOT}/tools/asite_ip/run.sh" || return $?

  python "${PROJECT_ROOT}/tools/asite_ip/convert.py" \
    --test "${DATASET_ROOT}/splits/within_sample/test.parquet" \
    --baseline-root "${NATIVE_OUTPUT_ROOT}/asite_ip" \
    --output "$RESULT_WITHIN_SAMPLE_DIR/asite_ip/candidate_predictions.parquet" || return $?

  prediction_rows="$(
    python - "$RESULT_WITHIN_SAMPLE_DIR/asite_ip/candidate_predictions.parquet" "$RESULT_WITHIN_SAMPLE_DIR/asite_ip/read_predictions.parquet" <<'PYCODE'
import sys

import pandas as pd

predictions = pd.read_parquet(sys.argv[1])
if predictions.empty:
    pd.DataFrame(
        columns=["read_id", "pred_offset", "pred_frame", "p_pred", "true_offset", "true_frame", "q_true"]
    ).to_parquet(sys.argv[2], index=False)
print(len(predictions))
PYCODE
  )" || return $?
  if [[ "$prediction_rows" == "0" ]]; then
    return 0
  fi

  python "${PROJECT_ROOT}/tabriboicl.py" evaluate \
    --test "${DATASET_ROOT}/splits/within_sample/test.parquet" \
    --predictions "$RESULT_WITHIN_SAMPLE_DIR/asite_ip/candidate_predictions.parquet" \
    --read-output "$RESULT_WITHIN_SAMPLE_DIR/asite_ip/read_predictions.parquet" \
    --model asite_ip \
    --split "p_site_offset_prediction_landmark_riboseq_within_sample" || return $?
}

if external_tool_enabled asite_ip; then
  run_external_step \
    "asite_ip" \
    "$RESULT_WITHIN_SAMPLE_DIR/asite_ip/read_predictions.parquet" \
    run_asite_ip_with_conversion
else
  record_external_status "asite_ip" "skipped" "not listed in EXTERNAL_TOOLS_STR"
fi

run_external_step \
  "external_summary" \
  "$RESULT_WITHIN_SAMPLE_DIR/.external_summary_prepared.ok" \
  python "${PROJECT_ROOT}/tools/shared/summarize_external_status.py" \
    --baseline-root "$NATIVE_OUTPUT_ROOT" \
    --metrics-dir "$RESULT_WITHIN_SAMPLE_DIR"
