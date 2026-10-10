#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
DATASET_ROOT="${DATASET_ROOT:-${TABRIBO_DATA_ROOT:-${PROJECT_ROOT}/data}/datasets/p_site_offset_prediction/landmark_riboseq}"
SAMPLE_TABLE="${SAMPLE_TABLE:-${DATASET_ROOT}/samples.tsv}"
OUTDIR="${OUTDIR:-${DATASET_ROOT}/views/transcriptome_bowtie2}"
TRANSCRIPT_FASTA_DIR="${TRANSCRIPT_FASTA_DIR:-${DATASET_ROOT}/views/external_transcripts}"
THREADS="${THREADS:-4}"
MAX_READS="${MAX_READS:-250000}"
case "$MAX_READS" in
  all|full|none) MAX_READS=0 ;;
esac
FORCE="${FORCE:-1}"

if [[ -f "${PROJECT_ROOT}/scripts/local_env.sh" ]]; then
  source "${PROJECT_ROOT}/scripts/local_env.sh"
fi

mkdir -p "$OUTDIR" "logs/baselines"

for ref_key in dmel scer; do
  ref_fa="${TRANSCRIPT_FASTA_DIR}/${ref_key}.fa"
  idx_prefix="${TRANSCRIPT_FASTA_DIR}/${ref_key}.bowtie2"
  if [[ ! -s "$ref_fa" ]]; then
    printf 'warning: %s is missing; using unflanked reference cDNA\n' "$ref_fa" >&2
    ref_fa="${DATASET_ROOT}/references/${ref_key}/cdna.fa"
    idx_prefix="${DATASET_ROOT}/references/${ref_key}/cdna_bowtie2"
  fi
  if [[ ! -f "${idx_prefix}.1.bt2" ]]; then
    bowtie2-build "$ref_fa" "$idx_prefix" > "logs/baselines/${ref_key}.cdna_bowtie2_build.log" 2>&1
  fi
done

tail -n +2 "$SAMPLE_TABLE" | while IFS=$'\t' read -r sample organism ref_key study_key sra_gb rest; do
  [[ -z "$sample" ]] && continue
  genome_bam="${DATASET_ROOT}/alignments/genome_bams/${sample}.Aligned.sortedByCoord.out.bam"
  out_bam="${OUTDIR}/${sample}.transcriptome.bam"
  tmp_fastq="${OUTDIR}/${sample}.mapped.fastq"
  tmp_sam="${OUTDIR}/${sample}.transcriptome.sam"

  if [[ "$FORCE" != "1" && -s "$out_bam" && -s "${out_bam}.bai" ]]; then
    printf '[%s] skipping %s transcriptome BAM\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$sample"
    continue
  fi
  if [[ ! -s "$genome_bam" ]]; then
    printf 'missing genome BAM for %s: %s\n' "$sample" "$genome_bam" >&2
    exit 1
  fi

  printf '[%s] preparing transcriptome BAM for %s (%s)\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$sample" "$ref_key"
  python "${PROJECT_ROOT}/tools/shared/bam_to_fastq_subset.py" \
    --bam "$genome_bam" \
    --output "$tmp_fastq" \
    --max-reads "$MAX_READS" \
    > "logs/baselines/${sample}.bam_to_fastq.log" 2>&1

  idx_prefix="${TRANSCRIPT_FASTA_DIR}/${ref_key}.bowtie2"
  if [[ ! -f "${idx_prefix}.1.bt2" ]]; then
    idx_prefix="${DATASET_ROOT}/references/${ref_key}/cdna_bowtie2"
  fi

  bowtie2 -x "$idx_prefix" \
    -U "$tmp_fastq" \
    --no-unal \
    -p "$THREADS" \
    -S "$tmp_sam" \
    2> "logs/baselines/${sample}.transcriptome_bowtie2.log"

  samtools view -@ "$THREADS" -bS "$tmp_sam" \
    | samtools sort -@ "$THREADS" -o "$out_bam"
  samtools index "$out_bam"
  rm -f "$tmp_fastq" "$tmp_sam"
done
