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
REFERENCE_DIR="${REFERENCE_DIR:-${DATASET_ROOT}/references}"
TRANSCRIPT_FASTA_DIR="${TRANSCRIPT_FASTA_DIR:-${DATASET_ROOT}/views/external_transcripts}"
OUTDIR="${OUTDIR:-${NATIVE_OUTPUT_ROOT}/asite_ip}"
LOG_DIR="${LOG_DIR:-${PROJECT_ROOT}/logs/baselines}"
ASITE_IP_SRC="${ASITE_IP_SRC:-}"
ASITE_IP_PYTHON="${ASITE_IP_PYTHON:-$(command -v python || command -v python3 || true)}"
ASITE_IP_FRAG_MIN="${ASITE_IP_FRAG_MIN:-20}"
ASITE_IP_FRAG_MAX="${ASITE_IP_FRAG_MAX:-75}"
ASITE_IP_MIN_AVG_READS="${ASITE_IP_MIN_AVG_READS:-1}"
ASITE_IP_UNIQUE_THRESHOLD="${ASITE_IP_UNIQUE_THRESHOLD:-70}"
ASITE_IP_SECOND_THRESHOLD="${ASITE_IP_SECOND_THRESHOLD:-5}"
ASITE_IP_FORCE="${ASITE_IP_FORCE:-1}"
ASITE_IP_KEEP_WORKDIR="${ASITE_IP_KEEP_WORKDIR:-0}"

if [[ -n "$ASITE_IP_PYTHON" && ! -x "$ASITE_IP_PYTHON" ]]; then
  RESOLVED_ASITE_IP_PYTHON="$(command -v "$ASITE_IP_PYTHON" || true)"
  if [[ -n "$RESOLVED_ASITE_IP_PYTHON" ]]; then
    ASITE_IP_PYTHON="$RESOLVED_ASITE_IP_PYTHON"
  fi
fi

if [[ -z "$ASITE_IP_SRC" ]]; then
  printf 'ASITE_IP_SRC is required. Set it to an external Asite_IP_method checkout; upstream has no visible license, so this repo does not vendor it.\n' >&2
  exit 1
fi

PREPROCESS_SRC="${ASITE_IP_SRC}/src/asite_ip_preprocess.py"
RUN_SRC="${ASITE_IP_SRC}/src/asite_ip_run.py"
if [[ ! -s "$PREPROCESS_SRC" || ! -s "$RUN_SRC" ]]; then
  printf 'Missing Asite IP upstream scripts under %s/src. Expected asite_ip_preprocess.py and asite_ip_run.py.\n' "$ASITE_IP_SRC" >&2
  exit 1
fi

if [[ -z "$ASITE_IP_PYTHON" || ! -x "$ASITE_IP_PYTHON" ]]; then
  printf 'No usable Python found. Set ASITE_IP_PYTHON or source scripts/local_env.sh.\n' >&2
  exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
  printf 'samtools is required for BAM-to-SAM conversion. Source scripts/local_env.sh or put samtools on PATH.\n' >&2
  exit 1
fi

if ! "$ASITE_IP_PYTHON" -c 'import Bio, numpy, scipy' >/dev/null 2>&1; then
  printf 'ASITE_IP_PYTHON must import Bio, numpy, and scipy.\n' >&2
  exit 1
fi

mkdir -p "$OUTDIR" "$LOG_DIR"
absolute_file() {
  local path="$1"
  case "$path" in
    /*) printf '%s
' "$path" ;;
    *) printf '%s/%s
' "$(cd "$(dirname "$path")" && pwd)" "$(basename "$path")" ;;
  esac
}
absolute_dir() {
  local path="$1"
  case "$path" in
    /*) printf '%s
' "$path" ;;
    *) (cd "$path" && pwd) ;;
  esac
}
SAMPLE_TABLE="$(absolute_file "$SAMPLE_TABLE")"
BAM_DIR="$(absolute_dir "$BAM_DIR")"
ANNOT_DIR="$(absolute_dir "$ANNOT_DIR")"
REFERENCE_DIR="$(absolute_dir "$REFERENCE_DIR")"
TRANSCRIPT_FASTA_DIR="$(absolute_dir "$TRANSCRIPT_FASTA_DIR")"
OUTDIR="$(absolute_dir "$OUTDIR")"
LOG_DIR="$(absolute_dir "$LOG_DIR")"
STATUS="${OUTDIR}/asite_ip_status.tsv"
printf 'sample_id	status	message
' > "$STATUS"

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/tabribo_asite_ip.XXXXXX")"
cleanup() {
  if [[ "$ASITE_IP_KEEP_WORKDIR" == "1" ]]; then
    printf 'Keeping Asite IP temporary workdir: %s\n' "$WORKDIR" >&2
  else
    rm -rf "$WORKDIR"
  fi
}
trap cleanup EXIT

COMPAT_DIR="${WORKDIR}/py3"
mkdir -p "$COMPAT_DIR"
if ! "$ASITE_IP_PYTHON" -m lib2to3 -w -n -o "$COMPAT_DIR" -W "$PREPROCESS_SRC" "$RUN_SRC" \
    > "${LOG_DIR}/asite_ip_lib2to3.stdout.log" \
    2> "${LOG_DIR}/asite_ip_lib2to3.stderr.log"; then
  printf 'Failed to convert upstream Asite IP scripts to Python 3; see %s/asite_ip_lib2to3.stderr.log\n' "$LOG_DIR" >&2
  exit 1
fi

PREPROCESS_PY="${COMPAT_DIR}/asite_ip_preprocess.py"
RUN_PY="${COMPAT_DIR}/asite_ip_run.py"
if ! "$ASITE_IP_PYTHON" "${PROJECT_ROOT}/tools/asite_ip/patch_preprocess.py" "$PREPROCESS_PY"; then
  printf 'Failed to apply the bounded-memory patch to converted Asite IP preprocessing.\n' >&2
  exit 1
fi
if ! "$ASITE_IP_PYTHON" - "$RUN_PY" <<'PYCODE'
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
text = path.read_text()
old = "elif perc_dict[sorted_perc[0]] >= offset_threshold and no_of_genes >= 10:"
new = "elif no_of_genes >= 10 and perc_dict[sorted_perc[0]] >= offset_threshold:"
if old not in text:
    raise SystemExit("expected Python 3 compatibility target not found in converted asite_ip_run.py")
path.write_text(text.replace(old, new, 1))
PYCODE
then
  printf 'Failed to patch converted Asite IP scripts for Python 3 compatibility.\n' >&2
  exit 1
fi
if ! "$ASITE_IP_PYTHON" "${PROJECT_ROOT}/tools/asite_ip/patch_run.py" "$RUN_PY"; then
  printf 'Failed to apply the compact-count patch to converted Asite IP runtime.\n' >&2
  exit 1
fi
if ! "$ASITE_IP_PYTHON" -m py_compile "$PREPROCESS_PY" "$RUN_PY" \
    > "${LOG_DIR}/asite_ip_py_compile.stdout.log" \
    2> "${LOG_DIR}/asite_ip_py_compile.stderr.log"; then
  printf 'Converted Asite IP scripts failed Python compilation; see %s/asite_ip_py_compile.stderr.log\n' "$LOG_DIR" >&2
  exit 1
fi

has_unique_offsets() {
  local result_file="$1"
  awk '
    BEGIN { in_table = 0; found = 0 }
    /^Frag size\tFrame_0\tFrame_1\tFrame_2$/ { in_table = 1; next }
    in_table && NF == 0 { exit }
    in_table && $1 ~ /^[0-9]+$/ {
      for (i = 2; i <= 4; i++) {
        if ($i ~ /^[0-9]+$/) {
          found = 1
        }
      }
    }
    END { exit(found ? 0 : 1) }
  ' "$result_file"
}

tail -n +2 "$SAMPLE_TABLE" | while IFS=$'\t' read -r sample organism ref_key study_key sra_gb rest; do
  [[ -z "$sample" ]] && continue

  bam="${BAM_DIR}/${sample}.transcriptome.bam"
  annotation="${ANNOT_DIR}/${ref_key}.asite_ip_start_index.tsv"
  transcripts="${TRANSCRIPT_FASTA_DIR}/${ref_key}.fa"
  if [[ ! -s "$transcripts" ]]; then
    transcripts="${REFERENCE_DIR}/${ref_key}/cdna.fa"
  fi
  sample_out="${OUTDIR}/${sample}"
  output_dir="${sample_out}/output"
  result_file="${output_dir}/Results_IP_algorithm.tab"
  multi_file="${output_dir}/Multiple_mapped_gene_read_counts_${ASITE_IP_FRAG_MIN}_${ASITE_IP_FRAG_MAX}.tab"
  sample_work="${WORKDIR}/samples/${sample}"
  sam_file="${sample_work}/${sample}.sam"

  if [[ "$ASITE_IP_FORCE" != "1" && -s "$result_file" ]]; then
    if has_unique_offsets "$result_file"; then
      printf '%s\tok\texisting native result kept at %s\n' "$sample" "$result_file" >> "$STATUS"
    else
      printf '%s\tempty_or_ambiguous\texisting native result has no unique numeric offsets\n' "$sample" >> "$STATUS"
    fi
    continue
  fi

  if [[ ! -s "$bam" ]]; then
    printf '%s\tfailed\tmissing transcriptome BAM\n' "$sample" >> "$STATUS"
    continue
  fi
  if [[ ! -s "$annotation" ]]; then
    printf '%s\tfailed\tmissing Asite IP annotation view for %s\n' "$sample" "$ref_key" >> "$STATUS"
    continue
  fi
  if [[ ! -s "$transcripts" ]]; then
    printf '%s\tfailed\tmissing transcript FASTA for %s\n' "$sample" "$ref_key" >> "$STATUS"
    continue
  fi

  rm -rf "$sample_out"
  mkdir -p "$sample_out" "$sample_work"

  if ! samtools view -h "$bam" \
      | awk 'BEGIN { OFS = "\t" } /^@/ { print; next } { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11 }' \
      > "$sam_file"; then
    rm -f "$sam_file"
    printf '%s\tfailed\tsamtools view failed for transcriptome BAM\n' "$sample" >> "$STATUS"
    continue
  fi

  if ! (
      cd "$sample_out"
      "$ASITE_IP_PYTHON" "$PREPROCESS_PY" \
        -i "$sam_file" \
        -a transcriptome \
        -g "$annotation" \
        -e "$transcripts" \
        -m "$ASITE_IP_FRAG_MIN" \
        -x "$ASITE_IP_FRAG_MAX" \
        -o "$sample_out"
    ) > "${LOG_DIR}/${sample}.asite_ip_preprocess.stdout.log" \
      2> "${LOG_DIR}/${sample}.asite_ip_preprocess.stderr.log"; then
    rm -f "$sam_file"
    printf '%s\tfailed\tpreprocess failed; see logs/baselines/%s.asite_ip_preprocess.stderr.log\n' "$sample" "$sample" >> "$STATUS"
    continue
  fi
  rm -f "$sam_file"

  if [[ ! -e "$multi_file" ]]; then
    printf '%s\tfailed\tpreprocess did not write %s\n' "$sample" "$multi_file" >> "$STATUS"
    continue
  fi

  if ! (
      cd "$sample_out"
      "$ASITE_IP_PYTHON" "$RUN_PY" \
        -i "${output_dir}/" \
        -j "$multi_file" \
        -m "$ASITE_IP_FRAG_MIN" \
        -x "$ASITE_IP_FRAG_MAX" \
        -o "$sample_out" \
        -n "$ASITE_IP_MIN_AVG_READS" \
        -t "$ASITE_IP_UNIQUE_THRESHOLD" \
        -s "$ASITE_IP_SECOND_THRESHOLD"
    ) > "${LOG_DIR}/${sample}.asite_ip_run.stdout.log" \
      2> "${LOG_DIR}/${sample}.asite_ip_run.stderr.log"; then
    printf '%s\tfailed\tnative IP run failed; see logs/baselines/%s.asite_ip_run.stderr.log\n' "$sample" "$sample" >> "$STATUS"
    continue
  fi

  if [[ ! -s "$result_file" ]]; then
    printf '%s\tfailed\tnative IP run did not write Results_IP_algorithm.tab\n' "$sample" >> "$STATUS"
  elif has_unique_offsets "$result_file"; then
    printf '%s\tok\tnative result written to %s\n' "$sample" "$result_file" >> "$STATUS"
  else
    printf '%s\tempty_or_ambiguous\tnative result has no unique numeric offsets\n' "$sample" >> "$STATUS"
  fi
done
