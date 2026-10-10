#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
DATASET_ROOT="${DATASET_ROOT:-${TABRIBO_DATA_ROOT:-${PROJECT_ROOT}/data}/datasets/p_site_offset_prediction/landmark_riboseq}"
RESULT_ROOT="${RESULT_ROOT:-${PROJECT_ROOT}/results/p_site_offset_prediction/landmark_riboseq}"
NATIVE_OUTPUT_ROOT="${NATIVE_OUTPUT_ROOT:-${RESULT_ROOT}/native_outputs}"

if [[ -f "${PROJECT_ROOT}/scripts/local_env.sh" ]]; then
  source "${PROJECT_ROOT}/scripts/local_env.sh"
fi

TRAIN="${TRAIN:-${DATASET_ROOT}/splits/within_sample/train.parquet}"
TEST="${TEST:-${DATASET_ROOT}/splits/within_sample/test.parquet}"
OUTDIR="${OUTDIR:-${RESULT_ROOT}/within_sample}"
TEST_READS="${TEST_READS:-${COMMON_TEST_READS:-5000}}"
COMMON_TEST_READS="${COMMON_TEST_READS:-$TEST_READS}"
TARGET_READS="${TARGET_READS:-}"
FEATURE_SET="${FEATURE_SET:-no_boundary}"
TABICL_ESTIMATORS="${TABICL_ESTIMATORS:-4}"
TABICL_DEVICE="${TABICL_DEVICE:-auto}"
TABICL_KV_CACHE="${TABICL_KV_CACHE:-kv}"
SPLIT_NAME="${SPLIT_NAME:-p_site_offset_prediction_landmark_riboseq_within_sample_${TEST_READS}}"
K_VALUES_STR="${K_VALUES_STR:-all}"
MODEL_FAMILIES_STR="${MODEL_FAMILIES_STR:-tabicl xgboost lightgbm}"
FORCE="${FORCE:-1}"
read -r -a K_VALUES <<< "$K_VALUES_STR"
read -r -a MODEL_FAMILIES <<< "$MODEL_FAMILIES_STR"

uses_all_test_reads() {
  case "$TEST_READS" in
    all|full|none|0) return 0 ;;
    *) return 1 ;;
  esac
}

mkdir -p "$OUTDIR"
export OUTDIR

run_eval() {
  local model_name="$1"
  local predictions="$2"
  local reads="$3"

  python "${PROJECT_ROOT}/tabriboicl.py" evaluate \
    --test "$TEST" \
    --predictions "$predictions" \
    --read-output "$reads" \
    --model "$model_name" \
    --split "$SPLIT_NAME"
}

model_enabled() {
  local target="$1"
  local model
  for model in "${MODEL_FAMILIES[@]}"; do
    [[ "$model" == "$target" ]] && return 0
  done
  return 1
}

predictions_match_target() {
  local predictions="$1"
  [[ -n "$TARGET_READS" ]] || return 0
  python - "$predictions" "$TARGET_READS" <<'PYCODE'
import sys
from pathlib import Path

import pandas as pd

from common.subsets import (
    assert_prediction_read_set_matches_target,
    load_read_subset,
)

predictions_path = Path(sys.argv[1])
target_path = Path(sys.argv[2])
predictions = pd.read_parquet(predictions_path, columns=["read_id"])
target_read_ids = load_read_subset(target_path)
assert_prediction_read_set_matches_target(predictions, target_read_ids)
PYCODE
}

should_reuse_predictions() {
  local model_name="$1"
  local train_spec="$2"
  local predictions="$3"
  local reads="$4"

  if [[ "$FORCE" != "1" && -s "$predictions" && -s "$reads" ]]; then
    if predictions_match_target "$predictions"; then
      printf '[%s] skipping %s train_reads=%s; predictions exist\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$model_name" "$train_spec"
      return 0
    fi
    printf '[%s] rerunning %s train_reads=%s; cached predictions do not match TARGET_READS\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$model_name" "$train_spec"
  fi
  return 1
}

run_tabicl_train() {
  local output="$1"
  local top_output="$2"
  local train_spec="$3"
  local args=(
    --train "$TRAIN"
    --test "$TEST"
    --output "$output"
    --top-output "$top_output"
    --feature-set "$FEATURE_SET"
    --n-estimators "$TABICL_ESTIMATORS"
    --device "$TABICL_DEVICE"
    --kv-cache "$TABICL_KV_CACHE"
  )
  if [[ "$train_spec" != "all" ]]; then
    args+=(--max-train-reads "$train_spec")
  fi
  if [[ -n "$TARGET_READS" ]]; then
    args+=(--target-read-ids "$TARGET_READS")
  elif ! uses_all_test_reads; then
    args+=(--max-test-reads "$TEST_READS")
  fi
  python "${PROJECT_ROOT}/tabriboicl.py" predict --model tabicl "${args[@]}"
}

run_gbdt_train() {
  local model="$1"
  local output="$2"
  local top_output="$3"
  local train_spec="$4"

  local args=(
    --train "$TRAIN"
    --test "$TEST"
    --output "$output"
    --top-output "$top_output"
    --model "$model"
    --feature-set "$FEATURE_SET"
  )
  if [[ "$train_spec" != "all" ]]; then
    args+=(--max-train-reads "$train_spec")
  fi
  if [[ -n "$TARGET_READS" ]]; then
    args+=(--target-read-ids "$TARGET_READS")
  elif ! uses_all_test_reads; then
    args+=(--max-test-reads "$TEST_READS")
  fi
  python "${PROJECT_ROOT}/tabriboicl.py" predict "${args[@]}"
}

for k in "${K_VALUES[@]}"; do
  case "$k" in
    all|full|none)
      train_spec="all"
      model_suffix="all"
      ;;
    *)
      train_spec="$k"
      model_suffix="k${k}"
      ;;
  esac

  printf '[%s] learned models train_reads=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$train_spec"

  if model_enabled tabicl; then
    tabicl_dir="$OUTDIR/tabicl_${model_suffix}"
    mkdir -p "$tabicl_dir"
    if should_reuse_predictions "tabicl" "$train_spec" "$tabicl_dir/candidate_predictions.parquet" "$tabicl_dir/read_predictions.parquet"; then
      :
    else
      run_tabicl_train \
        "$tabicl_dir/candidate_predictions.parquet" \
        "$tabicl_dir/top_predictions.parquet" \
        "$train_spec"
      run_eval "tabicl_${model_suffix}" \
        "$tabicl_dir/candidate_predictions.parquet" \
        "$tabicl_dir/read_predictions.parquet"
    fi
  fi

  for model in "${MODEL_FAMILIES[@]}"; do
    [[ "$model" == "tabicl" ]] && continue
    case "$model" in
      xgboost|lightgbm|logreg) ;;
      *)
        printf 'unsupported learned model family in MODEL_FAMILIES_STR: %s\n' "$model" >&2
        exit 2
        ;;
    esac
    model_dir="$OUTDIR/${model}_${model_suffix}"
    mkdir -p "$model_dir"
    if should_reuse_predictions "$model" "$train_spec" "$model_dir/candidate_predictions.parquet" "$model_dir/read_predictions.parquet"; then
      :
    else
      run_gbdt_train \
        "$model" \
        "$model_dir/candidate_predictions.parquet" \
        "$model_dir/top_predictions.parquet" \
        "$train_spec"
      run_eval "${model}_${model_suffix}" \
        "$model_dir/candidate_predictions.parquet" \
        "$model_dir/read_predictions.parquet"
    fi
  done
done
