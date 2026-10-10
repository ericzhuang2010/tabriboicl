#!/usr/bin/env bash
# Source this file to use the local TabRiboICL environment from this repo.

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ -x "${PROJECT_ROOT}/.conda-envs/tabribo-local/bin/python" ]]; then
  export PATH="${PROJECT_ROOT}/.conda-envs/tabribo-local/bin:${PATH}"
fi
export MPLCONFIGDIR="${MPLCONFIGDIR:-${PROJECT_ROOT}/.cache/matplotlib}"
export NCBI_HOME="${NCBI_HOME:-${PROJECT_ROOT}/.ncbi-local}"
export HF_HOME="${HF_HOME:-${PROJECT_ROOT}/.hf-cache}"
export TABRIBO_DATA_ROOT="${TABRIBO_DATA_ROOT:-${PROJECT_ROOT}/data}"
export TABRIBO_RESULT_ROOT="${TABRIBO_RESULT_ROOT:-${PROJECT_ROOT}/results}"
export TABRIBO_LOG_ROOT="${TABRIBO_LOG_ROOT:-${PROJECT_ROOT}/logs}"
export TABRIBO_EXTERNAL_TOOLS_ROOT="${TABRIBO_EXTERNAL_TOOLS_ROOT:-${PROJECT_ROOT}/external_tools}"
export TABRIBO_EXTERNAL_ENV_ROOT="${TABRIBO_EXTERNAL_ENV_ROOT:-${TABRIBO_EXTERNAL_TOOLS_ROOT}/envs}"
export TABRIBO_EXTERNAL_RLIB_ROOT="${TABRIBO_EXTERNAL_RLIB_ROOT:-${TABRIBO_EXTERNAL_TOOLS_ROOT}/Rlib}"
mkdir -p "$MPLCONFIGDIR"
mkdir -p "$NCBI_HOME"
mkdir -p "$HF_HOME"

# Runtime-specific OpenMP overrides, if needed, are set by the caller.
