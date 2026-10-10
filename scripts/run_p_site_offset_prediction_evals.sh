#!/usr/bin/env bash
set -euo pipefail

TABRIBO_DATA_ROOT="${TABRIBO_DATA_ROOT:-data}"
DATASET_ROOT="${DATASET_ROOT:-${TABRIBO_DATA_ROOT}/datasets/p_site_offset_prediction/landmark_riboseq}"
ROOT="${ROOT:-$DATASET_ROOT}"
RESULT_ROOT="${RESULT_ROOT:-results/p_site_offset_prediction/landmark_riboseq}"
NATIVE_OUTPUT_ROOT="${NATIVE_OUTPUT_ROOT:-${RESULT_ROOT}/native_outputs}"
BASELINE_ROOT="${BASELINE_ROOT:-$NATIVE_OUTPUT_ROOT}"
BASELINE_INPUT_ROOT="${BASELINE_INPUT_ROOT:-$ROOT/views}"
SAMPLE_TABLE="${SAMPLE_TABLE:-$ROOT/samples.tsv}"
EVAL_CONFIG="${EVAL_CONFIG:-configs/p_site_offset_prediction_landmark_riboseq.tsv}"
TEST_READS="${TEST_READS:-5000}"
TABICL_ESTIMATORS="${TABICL_ESTIMATORS:-4}"
TABICL_DEVICE="${TABICL_DEVICE:-auto}"
TABICL_KV_CACHE="${TABICL_KV_CACHE:-kv}"
FEATURE_SET="${FEATURE_SET:-no_boundary}"
FORCE="${FORCE:-1}"
RUN_EXTERNAL_TOOLS="${RUN_EXTERNAL_TOOLS:-1}"
TRANSCRIPTOME_MAX_READS="${TRANSCRIPTOME_MAX_READS:-all}"
EXTERNAL_STRICT="${EXTERNAL_STRICT:-1}"
RIBOBA_RUN_PIPELINE="${RIBOBA_RUN_PIPELINE:-1}"
CLEAN_STALE_RESULT_TABLES="${CLEAN_STALE_RESULT_TABLES:-1}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      printf 'Usage: %s [--config configs/p_site_offset_prediction_landmark_riboseq.tsv]\n' "$0"
      printf 'Run from the repository root after preparing the dataset; see docs/reproduction.md.\n'
      exit 0
      ;;
    --config)
      EVAL_CONFIG="$2"
      shift 2
      ;;
    *)
      if [[ $# -eq 1 ]]; then
        EVAL_CONFIG="$1"
        shift
      else
        printf 'unknown argument: %s\n' "$1" >&2
        exit 2
      fi
      ;;
  esac
done

if [[ -f "$(dirname "$0")/local_env.sh" ]]; then
  source "$(dirname "$0")/local_env.sh"
fi

COMMON_TEST_READS="${COMMON_TEST_READS:-${TEST_READS:-5000}}"
COMMON_TEST_SEED="${COMMON_TEST_SEED:-20260518}"

if [[ ! -s "$EVAL_CONFIG" ]]; then
  printf 'missing eval config: %s\n' "$EVAL_CONFIG" >&2
  exit 2
fi

mkdir -p "$RESULT_ROOT"

config_requests_asite_ip() {
  local experiment enabled algorithms rest algorithm algorithm_words
  while IFS=$'\t' read -r experiment enabled algorithms rest; do
    [[ -n "${experiment:-}" ]] || continue
    [[ "$experiment" == \#* ]] && continue
    [[ "$experiment" == "experiment" ]] && continue
    case "$enabled" in
      1|true|TRUE|True|yes|YES|Yes) ;;
      *) continue ;;
    esac
    algorithm_words="${algorithms//,/ }"
    algorithm_words="${algorithm_words//;/ }"
    for algorithm in $algorithm_words; do
      case "$algorithm" in
        asite_ip|all|external)
          return 0
          ;;
      esac
    done
  done < "$EVAL_CONFIG"
  return 1
}

if [[ "$RUN_EXTERNAL_TOOLS" == "1" ]] && config_requests_asite_ip && [[ -z "${ASITE_IP_SRC:-}" ]]; then
  printf 'ASITE_IP_SRC is required because %s requests asite_ip.\n' "$EVAL_CONFIG" >&2
  printf 'Example:\n' >&2
  printf '  mkdir -p external_tools/src\n' >&2
  printf '  git clone https://github.com/nabeel-bioinfo/Asite_IP_method external_tools/src/Asite_IP_method\n' >&2
  printf '  ASITE_IP_SRC=external_tools/src/Asite_IP_method bash %s --config %s\n' "$0" "$EVAL_CONFIG" >&2
  exit 2
fi


cleanup_stale_result_tables() {
  [[ "$CLEAN_STALE_RESULT_TABLES" == "1" ]] || return 0
  [[ -d "$RESULT_ROOT" ]] || return 0

  find "$RESULT_ROOT" -type f \( \
    -name "metrics.csv" -o \
    -name "metrics.tsv" -o \
    -name "metrics_model_subset.csv" -o \
    -name "model_summary.csv" -o \
    -name "model_summary.tsv" -o \
    -name "offset_baseline_summary.csv" -o \
    -name "offset_baseline_summary.tsv" -o \
    -name "offset_baseline_model_subset_summary.csv" -o \
    -name "offset_baseline_model_subset_summary.tsv" -o \
    -name "summary.csv" -o \
    -name "summary_common_test_subset.csv" -o \
    -name "leave_study_model_metrics.tsv" -o \
    -name "leave_study_offset_baseline_metrics.tsv" -o \
    -name "leave_study_offset_baseline_metrics_model_subset.tsv" \
  \) -delete
}

run_offsets() {
  local train="$1"
  local test="$2"
  local outdir="$3"
  local split="$4"
  local models_str="$5"
  read -r -a offset_models <<< "$models_str"

  python tabriboicl.py baselines \
    --train "$train" \
    --test "$test" \
    --outdir "$outdir" \
    --split "$split" \
    --models "${offset_models[@]}"
}


prepare_common_test_subset() {
  local test="$1"
  local outdir="$2"
  local subset_path="$outdir/common_test_subset_reads.tsv"
  mkdir -p "$outdir"

  python - "$test" "$subset_path" "$COMMON_TEST_READS" "$COMMON_TEST_SEED" <<'PYCODE'
import os
import sys
from pathlib import Path

import pandas as pd

from common.subsets import (
    read_subset_metadata,
    select_read_ids,
    write_read_subset,
)

test_path = Path(sys.argv[1])
subset_path = Path(sys.argv[2])
requested = sys.argv[3]
seed = int(sys.argv[4])

test_reads = pd.read_parquet(test_path, columns=["read_id"])
read_ids = select_read_ids(test_reads, requested, seed)
metadata = read_subset_metadata(
    test_path=test_path,
    test_reads_requested=requested,
    test_reads_selected=len(read_ids),
    selection_seed=seed,
    created_by="scripts/run_p_site_offset_prediction_evals.sh",
    extra={
        "COMMON_TEST_READS": requested,
        "COMMON_TEST_SEED": seed,
        "TEST_READS": os.environ.get("TEST_READS"),
        "test_table_rows": int(len(test_reads)),
        "test_table_unique_reads": int(test_reads["read_id"].astype(str).nunique()),
    },
)
write_read_subset(subset_path, read_ids, metadata)
PYCODE

  printf '%s\n' "$subset_path"
}

list_to_words() {
  local value="$1"
  value="${value//,/ }"
  value="${value//;/ }"
  printf '%s\n' "$value"
}

select_algorithms() {
  local algorithms="$1"
  shift
  local selected=""
  local algorithm
  local allowed
  for algorithm in $(list_to_words "$algorithms"); do
    case "$algorithm" in
      all)
        algorithm="fixed_offset read_length_offset tabicl xgboost lightgbm psite ribowaltz ribodeblur riboba scikit_ribo asite_ip"
        ;;
      offset_baselines)
        algorithm="fixed_offset read_length_offset"
        ;;
      learned)
        algorithm="tabicl xgboost lightgbm"
        ;;
      external)
        algorithm="psite ribowaltz ribodeblur riboba scikit_ribo asite_ip"
        ;;
    esac
    for allowed in "$@"; do
      case " $algorithm " in
        *" $allowed "*) selected="${selected}${selected:+ }${allowed}" ;;
      esac
    done
  done
  printf '%s\n' "$selected"
}

run_configured_experiment() {
  local experiment="$1"
  local train="$2"
  local test="$3"
  local outdir="$4"
  local split_name="$5"
  local algorithms="$6"
  local train_specs="$7"
  local target_reads
  target_reads="$(prepare_common_test_subset "$test" "$outdir")"

  local offset_models
  offset_models="$(select_algorithms "$algorithms" fixed_offset read_length_offset)"
  if [[ -n "$offset_models" ]]; then
    printf '[%s] %s offset baselines: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$experiment" "$offset_models"
    run_offsets "$train" "$test" "$outdir" "$split_name" "$offset_models"
  fi

  local learned_models
  learned_models="$(select_algorithms "$algorithms" tabicl xgboost lightgbm logreg)"
  if [[ -n "$learned_models" ]]; then
    printf '[%s] %s learned models: %s; train_specs=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$experiment" "$learned_models" "$train_specs"
    TRAIN="$train" \
    TEST="$test" \
    OUTDIR="$outdir" \
    COMMON_TEST_READS="$COMMON_TEST_READS" \
    TEST_READS="$COMMON_TEST_READS" \
    TARGET_READS="$target_reads" \
    K_VALUES_STR="$(list_to_words "$train_specs")" \
    MODEL_FAMILIES_STR="$learned_models" \
    SPLIT_NAME="${split_name}_${COMMON_TEST_READS}" \
    TABICL_ESTIMATORS="$TABICL_ESTIMATORS" \
    TABICL_DEVICE="$TABICL_DEVICE" \
    TABICL_KV_CACHE="$TABICL_KV_CACHE" \
    FEATURE_SET="$FEATURE_SET" \
    FORCE="$FORCE" \
    bash tools/shared/run_candidate_models.sh
  fi

  local external_tools
  external_tools="$(select_algorithms "$algorithms" psite ribowaltz ribodeblur riboba scikit_ribo asite_ip)"
  if [[ "$experiment" == "within_sample" && -n "$external_tools" ]]; then
    if [[ "$RUN_EXTERNAL_TOOLS" == "1" ]]; then
      printf '[%s] within-sample external tools: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$external_tools"
      ROOT="$ROOT" \
      DATASET_ROOT="$ROOT" \
      RESULT_ROOT="$RESULT_ROOT" \
      BASELINE_ROOT="$BASELINE_ROOT" \
      NATIVE_OUTPUT_ROOT="$NATIVE_OUTPUT_ROOT" \
      SAMPLE_TABLE="$SAMPLE_TABLE" \
      THREADS="${THREADS:-4}" \
      TRANSCRIPTOME_MAX_READS="$TRANSCRIPTOME_MAX_READS" \
      EXTERNAL_TOOLS_STR="$external_tools" \
      EXTERNAL_STRICT="$EXTERNAL_STRICT" \
      RIBOBA_RUN_PIPELINE="$RIBOBA_RUN_PIPELINE" \
      FORCE="$FORCE" \
        bash scripts/run_external_tools.sh
    else
      printf '[%s] skipping within-sample external tools; RUN_EXTERNAL_TOOLS=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$RUN_EXTERNAL_TOOLS"
    fi
  elif [[ "$experiment" == leave_one_*/* && -n "$external_tools" ]]; then
    printf '[%s] %s external direct-predictor rows will be materialized during summary when artifacts exist: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$experiment" "$external_tools"
  fi
}

dispatch_config_row() {
  local experiment="$1"
  local algorithms="$2"
  local train_specs="$3"

  train_specs="${train_specs:-all}"
  case "$experiment" in
    within_sample)
      run_configured_experiment \
        "$experiment" \
        "$ROOT/splits/within_sample/train.parquet" \
        "$ROOT/splits/within_sample/test.parquet" \
        "$RESULT_ROOT/within_sample" \
        "p_site_offset_prediction_landmark_riboseq_within_sample" \
        "$algorithms" \
        "$train_specs"
      ;;
    leave_one_*/*)
      local split_family="${experiment%%/*}"
      local selector="${experiment#*/}"
      if [[ "$split_family" == "leave_one_study" ]]; then
        split_family="leave_one_study_id"
      fi
      local split_dir
      for split_dir in "$ROOT"/splits/"$split_family"/$selector; do
        [[ -d "$split_dir" ]] || continue
        local heldout
        heldout="$(basename "$split_dir")"
        run_configured_experiment \
          "$split_family/$heldout" \
          "$split_dir/train.parquet" \
          "$split_dir/test.parquet" \
          "$RESULT_ROOT/$split_family/$heldout" \
          "p_site_offset_prediction_landmark_riboseq_${split_family}_${heldout}" \
          "$algorithms" \
          "$train_specs"
      done
      ;;
    *)
      printf 'unsupported experiment in %s: %s\n' "$EVAL_CONFIG" "$experiment" >&2
      exit 2
      ;;
  esac
}

cleanup_stale_result_tables

while IFS=$'\t' read -r experiment enabled algorithms train_specs compare_train_spec rest; do
  [[ -n "${experiment:-}" ]] || continue
  [[ "$experiment" == \#* ]] && continue
  [[ "$experiment" == "experiment" ]] && continue
  case "$enabled" in
    1|true|TRUE|True|yes|YES|Yes) ;;
    *) continue ;;
  esac
  dispatch_config_row "$experiment" "$algorithms" "$train_specs"
done < "$EVAL_CONFIG"

ROOT="$ROOT" RESULT_ROOT="$RESULT_ROOT" BASELINE_ROOT="$BASELINE_ROOT" BASELINE_INPUT_ROOT="$BASELINE_INPUT_ROOT" SAMPLE_TABLE="$SAMPLE_TABLE" EVAL_CONFIG="$EVAL_CONFIG" NATIVE_OUTPUT_ROOT="$NATIVE_OUTPUT_ROOT" COMMON_TEST_READS="$COMMON_TEST_READS" COMMON_TEST_SEED="$COMMON_TEST_SEED" \
  python scripts/summarize_p_site_offset_prediction_evals.py
