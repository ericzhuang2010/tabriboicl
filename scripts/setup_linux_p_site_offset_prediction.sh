#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DOWNLOAD_REFS=0
INSECURE_SSL=0

usage() {
  printf 'Usage: %s [--download-refs] [--insecure-ssl]\n' "$0"
  printf 'Prepare a fresh Linux machine to rerun scripts/run_p_site_offset_prediction_evals.sh.\n'
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --download-refs)
      DOWNLOAD_REFS=1
      shift
      ;;
    --insecure-ssl)
      INSECURE_SSL=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      usage
      exit 2
      ;;
  esac
done

if [[ "$(uname -s)" != "Linux" ]]; then
  printf 'This bootstrap script is only for Linux hosts.\n' >&2
  exit 2
fi

cd "$PROJECT_ROOT"

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || MISSING_CMDS+=("$1")
}

MISSING_CMDS=()
for cmd in bash curl tar git python3; do
  require_cmd "$cmd"
done
require_cmd Rscript

if [[ ${#MISSING_CMDS[@]} -gt 0 ]]; then
  printf 'Missing required system commands: %s\n' "${MISSING_CMDS[*]}" >&2
  printf 'Install the Linux prerequisites first; see docs/installation.md.\n' >&2
  exit 1
fi

install_micromamba() {
  local micromamba="${PROJECT_ROOT}/.tools/bin/micromamba"
  local arch
  local subdir
  local tmp_dir

  if [[ -x "$micromamba" ]] && "$micromamba" --version >/dev/null 2>&1; then
    printf 'using existing micromamba at %s\n' "$micromamba"
    return 0
  fi

  arch="$(uname -m)"
  case "$arch" in
    x86_64|amd64)
      subdir="linux-64"
      ;;
    aarch64|arm64)
      subdir="linux-aarch64"
      ;;
    *)
      printf 'Unsupported Linux architecture for micromamba bootstrap: %s\n' "$arch" >&2
      exit 2
      ;;
  esac

  mkdir -p "${PROJECT_ROOT}/.tools/bin"
  tmp_dir="$(mktemp -d)"
  curl -Ls "https://micro.mamba.pm/api/micromamba/${subdir}/latest" \
    | tar -xvj -C "$tmp_dir" bin/micromamba >/dev/null
  mv "$tmp_dir/bin/micromamba" "$micromamba"
  chmod +x "$micromamba"
  rm -rf "$tmp_dir"
  printf 'installed micromamba at %s\n' "$micromamba"
}

setup_local_env() {
  if [[ "$INSECURE_SSL" == "1" ]]; then
    bash scripts/setup_local_env.sh --insecure-ssl
  else
    bash scripts/setup_local_env.sh
  fi

  if [[ ! -x ".conda-envs/tabribo-local/bin/bowtie" ]]; then
    printf 'installing missing bowtie into .conda-envs/tabribo-local\n'
    ./.tools/bin/micromamba install -y \
      -r ".micromamba-root" \
      -p ".conda-envs/tabribo-local" \
      -c conda-forge \
      -c bioconda \
      bowtie
  fi

  (
    source scripts/local_env.sh
    python -c 'import pandas, pyarrow, pysam, sklearn, pybedtools, gffutils, seaborn'
  )
}

setup_baseline_code() {
  bash scripts/fetch_tool_sources.sh --yes
}

setup_psite_env() {
  local local_python="${PROJECT_ROOT}/.conda-envs/tabribo-local/bin/python"
  local psite_python="${TABRIBO_EXTERNAL_ENV_ROOT:-${PROJECT_ROOT}/external_tools/envs}/psite/bin/python"

  if [[ -d "${TABRIBO_EXTERNAL_ENV_ROOT:-${PROJECT_ROOT}/external_tools/envs}/psite" ]]; then
    if ! "$psite_python" --version >/dev/null 2>&1; then
      printf 'Existing external_tools/envs/psite is not runnable on this Linux host.\n' >&2
      printf 'Remove external_tools/envs/psite and rerun this script.\n' >&2
      exit 1
    fi
  else
    "$local_python" -m venv "${TABRIBO_EXTERNAL_ENV_ROOT:-${PROJECT_ROOT}/external_tools/envs}/psite"
  fi

  "$psite_python" -m pip install --upgrade pip setuptools wheel
  "$psite_python" -m pip install --upgrade psite
  "$psite_python" -c 'import psite'
}

setup_r_baselines() {
  export R_LIBS_USER="${R_LIBS_USER:-${TABRIBO_EXTERNAL_RLIB_ROOT:-${PROJECT_ROOT}/external_tools/Rlib}}"
  export RIBOBA_LOCAL_SRC="${PROJECT_ROOT}/tools/riboba/src"
  local riboba_patch="${PROJECT_ROOT}/tools/riboba/patches/riboba_tabriboicl.patch"
  mkdir -p "$R_LIBS_USER"

  if git -C "$RIBOBA_LOCAL_SRC" apply --reverse --check "$riboba_patch"; then
    printf 'RiboBA compatibility patch is already applied.\n'
  elif git -C "$RIBOBA_LOCAL_SRC" apply --check "$riboba_patch"; then
    git -C "$RIBOBA_LOCAL_SRC" apply "$riboba_patch"
    printf 'Applied RiboBA compatibility patch.\n'
  else
    printf 'RiboBA source does not match the expected compatibility patch.\n' >&2
    exit 1
  fi

  Rscript - <<'EOF'
lib_dir <- Sys.getenv("R_LIBS_USER")
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_dir, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

ensure_pkg("pak")
ensure_pkg("devtools")
ensure_pkg("data.table")
ensure_pkg("BiocManager")

if (!requireNamespace("riboWaltz", quietly = TRUE)) {
  pak::pkg_install("LabTranslationalArchitectomics/riboWaltz")
}
if (!requireNamespace("RiboBA", quietly = TRUE)) {
  pak::pkg_install("Bai-JunYu/RiboBA")
}

local_riboba <- Sys.getenv("RIBOBA_LOCAL_SRC")
if (!dir.exists(local_riboba)) {
  stop("Local patched RiboBA source is missing: ", local_riboba)
}
install.packages(
  local_riboba,
  repos = NULL,
  type = "source",
  dependencies = FALSE
)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(riboWaltz))
suppressPackageStartupMessages(library(RiboBA))
EOF
}

download_refs() {
  bash scripts/download_references.sh dmel --yes
  bash scripts/download_references.sh scer --yes
}

install_micromamba
setup_local_env
setup_baseline_code
setup_psite_env
setup_r_baselines

if [[ "$DOWNLOAD_REFS" == "1" ]]; then
  download_refs
fi

printf '\nLinux bootstrap complete.\n'
printf 'If you copied data/references/, data/intermediates/alignments/genome_bams/, and data/datasets/p_site_offset_prediction/ from the source machine,\n'
printf 'you can now rerun after checking out Asite IP externally:\n'
printf '  mkdir -p external_tools/src\n'
printf '  git clone https://github.com/nabeel-bioinfo/Asite_IP_method external_tools/src/Asite_IP_method\n'
printf '  ASITE_IP_SRC=external_tools/src/Asite_IP_method TABICL_KV_CACHE=false bash scripts/run_p_site_offset_prediction_evals.sh\n'
