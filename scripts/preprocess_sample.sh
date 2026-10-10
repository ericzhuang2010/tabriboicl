#!/usr/bin/env bash
set -euo pipefail

SAMPLE="${1:-}"
ORGANISM="${2:-dmel}"
ADAPTER="${ADAPTER:-CTGTAGGCACCATCAAT}"
MIN_LEN="${MIN_LEN:-20}"
MAX_LEN="${MAX_LEN:-120}"
TRIM_5P="${TRIM_5P:-0}"
THREADS="${THREADS:-8}"
ALIGNER="${ALIGNER:-bowtie2}"
STAR_QUANT_MODE="${STAR_QUANT_MODE:-none}"

if [[ -f "$(dirname "$0")/local_env.sh" ]]; then
  source "$(dirname "$0")/local_env.sh"
fi

TABRIBO_DATA_ROOT="${TABRIBO_DATA_ROOT:-${PROJECT_ROOT}/data}"
RAW_ROOT="${RAW_ROOT:-${TABRIBO_DATA_ROOT}/raw}"
TRIMMED_ROOT="${TRIMMED_ROOT:-${TABRIBO_DATA_ROOT}/intermediates/trimmed_fastq}"
ALIGN_ROOT="${ALIGN_ROOT:-${TABRIBO_DATA_ROOT}/intermediates/alignments/genome_bams}"

if [[ -z "$SAMPLE" ]]; then
  printf 'Usage: %s SRRXXXXXXX [dmel]\n' "$0" >&2
  exit 2
fi

case "$ORGANISM" in
  dmel)
    REF_DIR="${TABRIBO_DATA_ROOT}/references/dmel"
    STAR_SA_INDEX_NBASES=12
    ;;
  scer)
    REF_DIR="${TABRIBO_DATA_ROOT}/references/scer"
    STAR_SA_INDEX_NBASES=10
    ;;
  *)
    printf 'Unsupported organism target: %s\n' "$ORGANISM" >&2
    exit 2
    ;;
esac

mkdir -p "$TRIMMED_ROOT" "$ALIGN_ROOT" logs "$REF_DIR"

CUTADAPT_ARGS=(-m "$MIN_LEN" -M "$MAX_LEN")
if [[ -n "$TRIM_5P" && "$TRIM_5P" != "0" && "$TRIM_5P" != "." ]]; then
  CUTADAPT_ARGS+=(-u "$TRIM_5P")
fi

cutadapt "${CUTADAPT_ARGS[@]}" -a "$ADAPTER" \
  -o "${TRIMMED_ROOT}/${SAMPLE}.trimmed.fastq.gz" \
  "${RAW_ROOT}/fastq/${SAMPLE}.fastq.gz" \
  > "logs/${SAMPLE}.cutadapt.log" 2>&1

if [[ ! -f "${REF_DIR}/ncrna.1.bt2" ]]; then
  bowtie2-build "${REF_DIR}/ncrna.fa" "${REF_DIR}/ncrna"
fi

bowtie2 -x "${REF_DIR}/ncrna" \
  -U "${TRIMMED_ROOT}/${SAMPLE}.trimmed.fastq.gz" \
  --un-gz "${TRIMMED_ROOT}/${SAMPLE}.no_ncRNA.fastq.gz" \
  -S "${ALIGN_ROOT}/${SAMPLE}.contaminants.sam" \
  2> "logs/${SAMPLE}.contaminants.log"

gzip -cd "${TRIMMED_ROOT}/${SAMPLE}.no_ncRNA.fastq.gz" > "${TRIMMED_ROOT}/${SAMPLE}.no_ncRNA.fastq"

if [[ "$ALIGNER" == "star" && ! -d "${REF_DIR}/star" ]]; then
  mkdir -p "${REF_DIR}/star"
  STAR --runThreadN "$THREADS" \
    --runMode genomeGenerate \
    --genomeDir "${REF_DIR}/star" \
    --genomeFastaFiles "${REF_DIR}/genome.fa" \
    --sjdbGTFfile "${REF_DIR}/annotation.gtf" \
    --sjdbOverhang 29 \
    --genomeSAindexNbases "$STAR_SA_INDEX_NBASES"
fi

case "$ALIGNER" in
  star)
    rm -rf "${ALIGN_ROOT}/${SAMPLE}._STARtmp"
    STAR_ARGS=(
      --runThreadN "$THREADS"
      --genomeDir "${REF_DIR}/star"
      --readFilesIn "${TRIMMED_ROOT}/${SAMPLE}.no_ncRNA.fastq"
      --outSAMtype BAM SortedByCoordinate
      --outSAMattributes All
      --outFilterMultimapNmax 1
      --outFileNamePrefix "${ALIGN_ROOT}/${SAMPLE}."
    )

    if [[ "$STAR_QUANT_MODE" != "none" ]]; then
      STAR_ARGS+=(--quantMode "$STAR_QUANT_MODE")
    fi

    STAR "${STAR_ARGS[@]}"
    ;;
  bowtie2)
    if [[ ! -f "${REF_DIR}/genome_bowtie2.1.bt2" ]]; then
      bowtie2-build "${REF_DIR}/genome.fa" "${REF_DIR}/genome_bowtie2"
    fi
    bowtie2 -x "${REF_DIR}/genome_bowtie2" \
      -U "${TRIMMED_ROOT}/${SAMPLE}.no_ncRNA.fastq" \
      -S "${ALIGN_ROOT}/${SAMPLE}.genome.sam" \
      -p "$THREADS" \
      2> "logs/${SAMPLE}.genome_bowtie2.log"
    samtools view -@ "$THREADS" -bS "${ALIGN_ROOT}/${SAMPLE}.genome.sam" \
      | samtools sort -@ "$THREADS" -o "${ALIGN_ROOT}/${SAMPLE}.Aligned.sortedByCoord.out.bam"
    rm -f "${ALIGN_ROOT}/${SAMPLE}.genome.sam"
    ;;
  *)
    printf 'Unsupported ALIGNER: %s\n' "$ALIGNER" >&2
    exit 2
    ;;
esac

samtools index "${ALIGN_ROOT}/${SAMPLE}.Aligned.sortedByCoord.out.bam"
