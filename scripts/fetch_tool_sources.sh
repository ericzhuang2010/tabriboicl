#!/usr/bin/env bash
set -euo pipefail

DEST_ROOT="tools"
YES=0

usage() {
  printf "Usage: %s --yes [--dest-root tools]\n" "$0"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dest-root)
      DEST_ROOT="$2"
      shift 2
      ;;
    --yes)
      YES=1
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

if [[ "$YES" -ne 1 ]]; then
  printf "Refusing external code download without explicit --yes.\n" >&2
  exit 3
fi

mkdir -p \
  "$DEST_ROOT/ribodeblur" \
  "$DEST_ROOT/scikit_ribo" \
  "$DEST_ROOT/riboba"

clone_or_skip() {
  local url="$1"
  local target="$2"
  local branch="${3:-}"

  if [[ -d "$target/.git" ]]; then
    printf "already present: %s\n" "$target"
    return 0
  fi

  if [[ -n "$branch" ]]; then
    git clone --depth 1 --branch "$branch" "$url" "$target"
  else
    git clone --depth 1 "$url" "$target"
  fi
}

clone_or_skip https://github.com/Kingsford-Group/ribodeblur.git "$DEST_ROOT/ribodeblur/src"
clone_or_skip https://github.com/schatzlab/scikit-ribo.git "$DEST_ROOT/scikit_ribo/src"
clone_or_skip https://github.com/Bai-JunYu/RiboBA.git "$DEST_ROOT/riboba/src"

printf "Tool source repositories are under %s/<tool>/src\n" "$DEST_ROOT"
