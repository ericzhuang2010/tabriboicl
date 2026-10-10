#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
DATASET_ROOT="${DATASET_ROOT:-${TABRIBO_DATA_ROOT:-${PROJECT_ROOT}/data}/datasets/p_site_offset_prediction/landmark_riboseq}"
RESULT_ROOT="${RESULT_ROOT:-${PROJECT_ROOT}/results/p_site_offset_prediction/landmark_riboseq}"
NATIVE_OUTPUT_ROOT="${NATIVE_OUTPUT_ROOT:-${RESULT_ROOT}/native_outputs}"

if [[ -f "${PROJECT_ROOT}/scripts/local_env.sh" ]]; then
  source "${PROJECT_ROOT}/scripts/local_env.sh"
fi

RIBODEBLUR_SRC="${RIBODEBLUR_SRC:-${PROJECT_ROOT}/tools/ribodeblur/src}"
SAMPLE_TABLE="${SAMPLE_TABLE:-${DATASET_ROOT}/samples.tsv}"
BAM_DIR="${BAM_DIR:-${DATASET_ROOT}/views/transcriptome_bowtie2}"
ANNOT_DIR="${ANNOT_DIR:-${DATASET_ROOT}/views/annotations}"
OUTDIR="${OUTDIR:-${NATIVE_OUTPUT_ROOT}/ribodeblur}"
PREDICTIONS="${PREDICTIONS:-${RESULT_ROOT}/within_sample/ribodeblur/candidate_predictions.parquet}"
READ_PREDICTIONS="${READ_PREDICTIONS:-${RESULT_ROOT}/within_sample/ribodeblur/read_predictions.parquet}"
TEST="${TEST:-${DATASET_ROOT}/splits/within_sample/test.parquet}"
SPLIT="${SPLIT:-p_site_offset_prediction_landmark_riboseq_within_sample}"
RIBODEBLUR_COVER_RATIO="${RIBODEBLUR_COVER_RATIO:-0.01}"
RIBODEBLUR_MIN_HIGH_COVER_READ_LENGTHS="${RIBODEBLUR_MIN_HIGH_COVER_READ_LENGTHS:-1}"
RIBODEBLUR_LOWEST_FRAME_CNT="${RIBODEBLUR_LOWEST_FRAME_CNT:-1}"
THREADS="${THREADS:-4}"
RIBODEBLUR_MAX_PROFILE_ALIGNMENTS="${RIBODEBLUR_MAX_PROFILE_ALIGNMENTS:-2000000}"
RIBODEBLUR_SAMPLE_SEED="${RIBODEBLUR_SAMPLE_SEED:-20260517}"

mkdir -p "$OUTDIR" "logs/baselines" "$(dirname "$PREDICTIONS")"
mkdir -p "$(dirname "$READ_PREDICTIONS")"
STATUS="${OUTDIR}/ribodeblur_status.tsv"
printf 'sample_id\tstatus\tmessage\n' > "$STATUS"

tail -n +2 "$SAMPLE_TABLE" | while IFS=$'\t' read -r sample organism ref_key study_key sra_gb rest; do
  [[ -z "$sample" ]] && continue
  bam="${BAM_DIR}/${sample}.transcriptome.bam"
  ref="${ANNOT_DIR}/${ref_key}.ribodeblur.fa"
  filtered_bam="${OUTDIR}/${sample}/${sample}.cds_transcripts.bam"
  prefix="${OUTDIR}/${sample}/${sample}"
  mkdir -p "${OUTDIR}/${sample}"

  if [[ ! -s "$bam" ]]; then
    printf '%s\tfailed\tmissing transcriptome BAM\n' "$sample" >> "$STATUS"
    continue
  fi

  python "${PROJECT_ROOT}/tools/ribodeblur/filter_bam_to_fasta_refs.py" \
    --bam "$bam" \
    --fasta "$ref" \
    --output "$filtered_bam" \
    > "logs/baselines/${sample}.ribodeblur_filter.log" 2>&1
  samtools index "$filtered_bam"

  profile_bam="$filtered_bam"
  mapped_alignments="$(samtools view -c -F 4 "$filtered_bam")"
  if [[ "$RIBODEBLUR_MAX_PROFILE_ALIGNMENTS" -gt 0 && "$mapped_alignments" -gt "$RIBODEBLUR_MAX_PROFILE_ALIGNMENTS" ]]; then
    sample_fraction="$(awk -v cap="$RIBODEBLUR_MAX_PROFILE_ALIGNMENTS" -v total="$mapped_alignments" 'BEGIN { printf "%.12f", cap / total }')"
    sample_digits="${sample_fraction#0.}"
    sample_spec="${RIBODEBLUR_SAMPLE_SEED}.${sample_digits}"
    profile_bam="${OUTDIR}/${sample}/${sample}.profile_input.bam"
    tmp_profile_bam="${profile_bam}.tmp.bam"
    tmp_profile_bai="${tmp_profile_bam}.bai"
    rm -f "$profile_bam" "${profile_bam}.bai" "$tmp_profile_bam" "$tmp_profile_bai"
    samtools view \
      -@ "$THREADS" \
      -b \
      -F 4 \
      -s "$sample_spec" \
      -o "$tmp_profile_bam" \
      "$filtered_bam"
    samtools index -@ "$THREADS" "$tmp_profile_bam"
    mv "$tmp_profile_bam" "$profile_bam"
    mv "$tmp_profile_bai" "${profile_bam}.bai"
  else
    sample_fraction="1.0"
  fi

  dominant_rlen="$(samtools view -F 2324 "$profile_bam" | awk '{print length($10)}' | sort | uniq -c | sort -nr | awk 'NR == 1 {print $2}')"
  max_rlen="$(samtools view -F 2324 "$profile_bam" | awk '{print length($10)}' | sort -nr | awk 'NR == 1 {print $1}')"
  if [[ -z "$dominant_rlen" || -z "$max_rlen" ]]; then
    printf '%s\tfailed\tno forward primary transcriptome alignments after filtering\n' "$sample" >> "$STATUS"
    continue
  fi
  rlen_max="$max_rlen"
  if [[ "$rlen_max" -lt 35 ]]; then
    rlen_max=35
  fi

  if RIBODEBLUR_ANCHOR_READ_LENGTH="$dominant_rlen" \
      RIBODEBLUR_RLEN_MAX="$rlen_max" \
      RIBODEBLUR_COVER_RATIO="$RIBODEBLUR_COVER_RATIO" \
      RIBODEBLUR_MIN_HIGH_COVER_READ_LENGTHS="$RIBODEBLUR_MIN_HIGH_COVER_READ_LENGTHS" \
      RIBODEBLUR_LOWEST_FRAME_CNT="$RIBODEBLUR_LOWEST_FRAME_CNT" \
      python "${RIBODEBLUR_SRC}/deblur_pipeline.py" \
      -r "$ref" \
      -b "$profile_bam" \
      -o "$prefix" \
      -f \
      > "logs/baselines/${sample}.ribodeblur.stdout.log" \
      2> "logs/baselines/${sample}.ribodeblur.stderr.log"; then
    printf '%s\tok\tprofile written to %s.profile; anchor_read_length=%s; cover_ratio=%s; profile_fraction=%s\n' "$sample" "$prefix" "$dominant_rlen" "$RIBODEBLUR_COVER_RATIO" "$sample_fraction" >> "$STATUS"
  else
    printf '%s\tfailed\tdeblur failed; see logs/baselines/%s.ribodeblur.stderr.log\n' "$sample" "$sample" >> "$STATUS"
  fi
done

python "${PROJECT_ROOT}/tools/shared/convert_profile_predictions.py" \
  --method ribodeblur \
  --test "$TEST" \
  --sample-table "$SAMPLE_TABLE" \
  --baseline-root "$OUTDIR" \
  --annotation-dir "$ANNOT_DIR" \
  --output "$PREDICTIONS"

python "${PROJECT_ROOT}/tabriboicl.py" evaluate \
  --test "$TEST" \
  --predictions "$PREDICTIONS" \
  --read-output "$READ_PREDICTIONS" \
  --model ribodeblur \
  --split "$SPLIT"
