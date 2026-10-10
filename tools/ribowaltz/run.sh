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
BAM_DIR="${BAM_DIR:-${NATIVE_OUTPUT_ROOT}/ribowaltz/bams}"
ANNOT_DIR="${ANNOT_DIR:-${DATASET_ROOT}/views/annotations}"
OUTDIR="${OUTDIR:-${NATIVE_OUTPUT_ROOT}/ribowaltz}"
PREDICTIONS="${PREDICTIONS:-${RESULT_ROOT}/within_sample/ribowaltz/candidate_predictions.parquet}"
READ_PREDICTIONS="${READ_PREDICTIONS:-${RESULT_ROOT}/within_sample/ribowaltz/read_predictions.parquet}"
TEST="${TEST:-${DATASET_ROOT}/splits/within_sample/test.parquet}"
SPLIT="${SPLIT:-p_site_offset_prediction_landmark_riboseq_within_sample}"

R_LIBS_USER="${R_LIBS_USER:-${TABRIBO_EXTERNAL_RLIB_ROOT:-${PROJECT_ROOT}/external_tools/Rlib}}"
export R_LIBS_USER

THREADS="${THREADS:-4}"
RIBOWALTZ_MAX_ALIGNMENTS="${RIBOWALTZ_MAX_ALIGNMENTS:-2000000}"
RIBOWALTZ_SAMPLE_SEED="${RIBOWALTZ_SAMPLE_SEED:-20260517}"

mkdir -p "$BAM_DIR" "$OUTDIR" "$(dirname "$PREDICTIONS")"
mkdir -p "$(dirname "$READ_PREDICTIONS")"

tail -n +2 "$SAMPLE_TABLE" | while IFS=$'\t' read -r sample organism ref_key study_key sra_gb rest; do
  [[ -z "$sample" ]] && continue
  source_bam="${DATASET_ROOT}/views/transcriptome_bowtie2/${sample}.transcriptome.bam"
  source_bai="${source_bam}.bai"
  target_bam="${BAM_DIR}/${sample}.bam"
  target_bai="${target_bam}.bai"
  mapped_alignments="$(samtools view -c -F 4 "$source_bam")"

  rm -f "$target_bam" "$target_bai"
  if [[ "$RIBOWALTZ_MAX_ALIGNMENTS" -gt 0 && "$mapped_alignments" -gt "$RIBOWALTZ_MAX_ALIGNMENTS" ]]; then
    sample_fraction="$(awk -v cap="$RIBOWALTZ_MAX_ALIGNMENTS" -v total="$mapped_alignments" 'BEGIN { printf "%.12f", cap / total }')"
    sample_digits="${sample_fraction#0.}"
    sample_spec="${RIBOWALTZ_SAMPLE_SEED}.${sample_digits}"
    tmp_bam="${target_bam}.tmp.bam"
    tmp_bai="${tmp_bam}.bai"
    rm -f "$tmp_bam" "$tmp_bai"
    samtools view \
      -@ "$THREADS" \
      -b \
      -F 4 \
      -s "$sample_spec" \
      -o "$tmp_bam" \
      "$source_bam"
    samtools index -@ "$THREADS" "$tmp_bam"
    mv "$tmp_bam" "$target_bam"
    mv "$tmp_bai" "$target_bai"
  else
    ln -s "$(cd "$(dirname "$source_bam")" && pwd)/$(basename "$source_bam")" "$target_bam"
    ln -s "$(cd "$(dirname "$source_bai")" && pwd)/$(basename "$source_bai")" "$target_bai"
  fi
done

Rscript "${PROJECT_ROOT}/tools/ribowaltz/run.R" "$SAMPLE_TABLE" "$BAM_DIR" "$ANNOT_DIR" "$OUTDIR"

python "${PROJECT_ROOT}/tools/ribowaltz/convert.py" \
  --test "$TEST" \
  --offset-dir "${OUTDIR}/offsets" \
  --output "$PREDICTIONS"

python "${PROJECT_ROOT}/tabriboicl.py" evaluate \
  --test "$TEST" \
  --predictions "$PREDICTIONS" \
  --read-output "$READ_PREDICTIONS" \
  --model ribowaltz \
  --split "$SPLIT"
