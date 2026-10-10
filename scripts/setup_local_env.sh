#!/usr/bin/env bash
set -euo pipefail

SSL_MODE="verified"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREFIX="${PROJECT_ROOT}/.conda-envs/tabribo-local"
ROOT_PREFIX="${PROJECT_ROOT}/.micromamba-root"
ENV_FILE="${PROJECT_ROOT}/environment.yml"

usage() {
  printf 'Usage: %s [--insecure-ssl]\n' "$0"
  printf 'Create or reuse the local environment, then install Python dependencies from pyproject.toml.\n'
  printf 'TLS verification is enabled; --insecure-ssl applies only to micromamba.\n'
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --insecure-ssl)
      SSL_MODE="insecure"
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

resolve_ca_bundle() {
  if [[ -n "${SSL_CERT_FILE:-}" && -r "${SSL_CERT_FILE}" ]]; then
    printf '%s\n' "${SSL_CERT_FILE}"
    return 0
  fi

  local candidate
  for candidate in \
    /etc/ssl/cert.pem \
    /etc/ssl/certs/ca-certificates.crt \
    /etc/pki/tls/certs/ca-bundle.crt \
    /etc/ssl/ca-bundle.pem
  do
    if [[ -r "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done

  return 1
}

if [[ -d "$PREFIX" ]]; then
  printf 'Using existing environment at %s\n' "$PREFIX"
else
  if [[ -x "${PROJECT_ROOT}/.tools/bin/micromamba" ]]; then
    MICROMAMBA="${PROJECT_ROOT}/.tools/bin/micromamba"
  elif command -v micromamba >/dev/null 2>&1; then
    MICROMAMBA="$(command -v micromamba)"
  else
    printf 'Missing micromamba. Put it on PATH or run scripts/setup_linux_p_site_offset_prediction.sh on Linux.\n' >&2
    exit 1
  fi

  if [[ "$SSL_MODE" == "insecure" ]]; then
    printf 'Using micromamba with --ssl-verify false because --insecure-ssl was passed.\n' >&2
    "$MICROMAMBA" create -y \
      -r "$ROOT_PREFIX" \
      -p "$PREFIX" \
      -f "$ENV_FILE" \
      --ssl-verify false
  elif CERT_BUNDLE="$(resolve_ca_bundle)"; then
    printf 'Using micromamba with verified TLS via %s.\n' "$CERT_BUNDLE" >&2
    SSL_CERT_FILE="$CERT_BUNDLE" CURL_CA_BUNDLE="$CERT_BUNDLE" \
      "$MICROMAMBA" create -y \
        -r "$ROOT_PREFIX" \
        -p "$PREFIX" \
        -f "$ENV_FILE" \
        --ssl-verify "$CERT_BUNDLE"
  else
    printf 'Using micromamba with verified TLS via its default trust store.\n' >&2
    "$MICROMAMBA" create -y -r "$ROOT_PREFIX" -p "$PREFIX" -f "$ENV_FILE"
  fi
fi

if ! "$PREFIX/bin/python" --version >/dev/null 2>&1; then
  printf 'Environment at %s has no runnable Python; repair it or move it aside and rerun setup.\n' "$PREFIX" >&2
  exit 1
fi

"$PREFIX/bin/python" -m pip install -e "${PROJECT_ROOT}[predictor,comparators,external]"

printf 'Environment ready at %s\n' "$PREFIX"
