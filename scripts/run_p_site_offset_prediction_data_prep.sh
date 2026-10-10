#!/usr/bin/env bash
set -euo pipefail

TABRIBO_DATA_ROOT="${TABRIBO_DATA_ROOT:-data}"
RAW_ROOT="${RAW_ROOT:-${TABRIBO_DATA_ROOT}/raw}"
TRIMMED_ROOT="${TRIMMED_ROOT:-${TABRIBO_DATA_ROOT}/intermediates/trimmed_fastq}"
ALIGN_ROOT="${ALIGN_ROOT:-${TABRIBO_DATA_ROOT}/intermediates/alignments/genome_bams}"
SAMPLE_TABLE="${SAMPLE_TABLE:-${TABRIBO_DATA_ROOT}/datasets/p_site_offset_prediction/landmark_riboseq/samples.tsv}"
THREADS=4
MIN_FREE_GB=30
MAX_READS=""

if [[ -f "$(dirname "$0")/local_env.sh" ]]; then
  source "$(dirname "$0")/local_env.sh"
fi

usage() {
  printf 'Usage: %s [--sample-table ${TABRIBO_DATA_ROOT}/datasets/p_site_offset_prediction/landmark_riboseq/samples.tsv] [--threads 4] [--max-reads N]\n' "$0"
  printf 'Default ALIGNER is bowtie2 for the local run. Set ALIGNER=star to test STAR.\n'
  printf 'Default cutadapt bounds are MIN_LEN=20 and MAX_LEN=120 for heterogeneous B1 reads.\n'
  printf 'Sample tables may include optional adapter, min_len, max_len, and trim5 columns.\n'
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-table)
      SAMPLE_TABLE="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --min-free-gb)
      MIN_FREE_GB="$2"
      shift 2
      ;;
    --max-reads)
      MAX_READS="$2"
      shift 2
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

command -v prefetch >/dev/null
command -v fasterq-dump >/dev/null
command -v cutadapt >/dev/null
command -v bowtie2 >/dev/null
command -v STAR >/dev/null
command -v samtools >/dev/null
command -v pigz >/dev/null

if [[ -n "$MAX_READS" ]]; then
  command -v seqtk >/dev/null
fi

check_storage() {
  local path="$TABRIBO_DATA_ROOT"
  [[ -e "$path" ]] || path="."
  python3 - "$MIN_FREE_GB" "$path" <<PY
import shutil, sys
min_free_gb = float(sys.argv[1])
usage = shutil.disk_usage(sys.argv[2])
free_gb = usage.free / 1024**3
total_gb = usage.total / 1024**3
used_gb = usage.used / 1024**3
print(f"disk total={total_gb:.1f} GiB used={used_gb:.1f} GiB free={free_gb:.1f} GiB")
if free_gb < min_free_gb:
    print(f"free space is below required threshold: {min_free_gb:.1f} GiB", file=sys.stderr)
    raise SystemExit(1)
PY
}

mkdir -p "${RAW_ROOT}/sra" "${RAW_ROOT}/fastq" "$TRIMMED_ROOT" "$ALIGN_ROOT" logs results/qc

tail -n +2 "$SAMPLE_TABLE" | while IFS=$'\t' read -r sample organism ref_key study_key sra_gb adapter min_len max_len trim5 _extra; do
  [[ -z "$sample" ]] && continue
  adapter="${adapter:-}"
  min_len="${min_len:-}"
  max_len="${max_len:-}"
  trim5="${trim5:-}"

  printf '[%s] starting %s (%s, %s)\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$sample" "$study_key" "$ref_key"
  check_storage

  if [[ -s "results/qc/${sample}.flagstat.txt" ]]; then
    printf '[%s] skipping %s; QC output already exists\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$sample"
    continue
  fi

  if [[ ! -f "${RAW_ROOT}/fastq/${sample}.fastq.gz" ]]; then
    if [[ ! -f "${RAW_ROOT}/fastq/${sample}.fastq" && ! -f "${RAW_ROOT}/fastq/${sample}_1.fastq" ]]; then
      prefetch "$sample" --output-directory ${RAW_ROOT}/sra
      fasterq-dump "${RAW_ROOT}/sra/${sample}" --split-files --threads "$THREADS" --outdir ${RAW_ROOT}/fastq
    fi

    if [[ -f "${RAW_ROOT}/fastq/${sample}.fastq" ]]; then
      fastq_in="${RAW_ROOT}/fastq/${sample}.fastq"
    elif [[ -f "${RAW_ROOT}/fastq/${sample}_1.fastq" ]]; then
      fastq_in="${RAW_ROOT}/fastq/${sample}_1.fastq"
    else
      printf 'Could not find FASTQ output for %s\n' "$sample" >&2
      exit 1
    fi

    if [[ -n "$MAX_READS" ]]; then
      seqtk sample -s100 "$fastq_in" "$MAX_READS" \
        | pigz -p "$THREADS" > "${RAW_ROOT}/fastq/${sample}.fastq.gz"
      rm -f "${RAW_ROOT}/fastq/${sample}.fastq" "${RAW_ROOT}/fastq/${sample}_1.fastq" "${RAW_ROOT}/fastq/${sample}_2.fastq"
    else
      if [[ "$fastq_in" != "${RAW_ROOT}/fastq/${sample}.fastq" ]]; then
        mv "$fastq_in" "${RAW_ROOT}/fastq/${sample}.fastq"
      fi
      rm -f "${RAW_ROOT}/fastq/${sample}_2.fastq"
      pigz -p "$THREADS" "${RAW_ROOT}/fastq/${sample}.fastq"
    fi
  fi

  preprocess_env=(THREADS="$THREADS")
  if [[ -n "$adapter" && "$adapter" != "." ]]; then
    preprocess_env+=(ADAPTER="$adapter")
  fi
  if [[ -n "$min_len" && "$min_len" != "." ]]; then
    preprocess_env+=(MIN_LEN="$min_len")
  fi
  if [[ -n "$max_len" && "$max_len" != "." ]]; then
    preprocess_env+=(MAX_LEN="$max_len")
  fi
  if [[ -n "$trim5" && "$trim5" != "." ]]; then
    preprocess_env+=(TRIM_5P="$trim5")
  fi

  env "${preprocess_env[@]}" scripts/preprocess_sample.sh "$sample" "$ref_key"

  samtools flagstat "${ALIGN_ROOT}/${sample}.Aligned.sortedByCoord.out.bam" \
    > "results/qc/${sample}.flagstat.txt"

  check_storage
  printf '[%s] finished %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$sample"
done
