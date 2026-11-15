#!/usr/bin/env bash
set -euo pipefail

RUN_LIST=""
OUT_DIR="data/raw"
THREADS=8
YES=0
MAX_RUNS=""

usage() {
  printf 'Usage: %s --run-list data/metadata/stage1_b1_core.runs.txt --out data/raw --threads 8 --yes [--max-runs N]\n' "$0"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-list)
      RUN_LIST="$2"
      shift 2
      ;;
    --out)
      OUT_DIR="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --max-runs)
      MAX_RUNS="$2"
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

if [[ -z "$RUN_LIST" ]]; then
  usage
  exit 2
fi

if [[ "$YES" -ne 1 ]]; then
  printf 'Refusing bulk download without explicit approval. Re-run with --yes only after approval.\n' >&2
  exit 3
fi

command -v prefetch >/dev/null
command -v fasterq-dump >/dev/null
command -v pigz >/dev/null

mkdir -p "$OUT_DIR/sra" "$OUT_DIR/fastq" logs

count=0
while IFS= read -r srr; do
  [[ -z "$srr" ]] && continue
  [[ "$srr" =~ ^# ]] && continue
  count=$((count + 1))
  if [[ -n "$MAX_RUNS" && "$count" -gt "$MAX_RUNS" ]]; then
    break
  fi

  printf '[%s] prefetch %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$srr"
  prefetch "$srr" --output-directory "$OUT_DIR/sra"

  printf '[%s] fasterq-dump %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$srr"
  fasterq-dump "$OUT_DIR/sra/$srr" --split-files --threads "$THREADS" --outdir "$OUT_DIR/fastq"

  printf '[%s] pigz %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$srr"
  pigz -p "$THREADS" "$OUT_DIR"/fastq/"$srr"*.fastq
done < "$RUN_LIST"
