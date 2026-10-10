#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(riboWaltz)
})

args <- commandArgs(trailingOnly = TRUE)
sample_table <- ifelse(length(args) >= 1, args[[1]], "data/datasets/p_site_offset_prediction/landmark_riboseq/samples.tsv")
bam_dir <- ifelse(length(args) >= 2, args[[2]], "results/p_site_offset_prediction/landmark_riboseq/native_outputs/ribowaltz/bams")
annot_dir <- ifelse(length(args) >= 3, args[[3]], "data/datasets/p_site_offset_prediction/landmark_riboseq/views/annotations")
outdir <- ifelse(length(args) >= 4, args[[4]], "results/p_site_offset_prediction/landmark_riboseq/native_outputs/ribowaltz")

dir.create(file.path(outdir, "offsets"), recursive = TRUE, showWarnings = FALSE)
dir.create("logs/baselines", recursive = TRUE, showWarnings = FALSE)

samples <- fread(sample_table)
status <- data.table(sample_id = character(), status = character(), message = character())

for (idx in seq_len(nrow(samples))) {
  sample_id <- samples$run[[idx]]
  ref_key <- samples$ref_key[[idx]]
  annotation_path <- file.path(annot_dir, paste0(ref_key, ".ribowaltz_annotation.tsv"))
  bam_path <- file.path(bam_dir, paste0(sample_id, ".bam"))
  if (!file.exists(bam_path)) {
    status <- rbind(status, data.table(sample_id = sample_id, status = "failed", message = "missing BAM symlink"))
    next
  }

  message(sprintf("[%s] riboWaltz %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), sample_id))
  result <- tryCatch({
    annotation <- fread(annotation_path)
    name_samples <- setNames(sample_id, sample_id)
    reads <- bamtolist(
      bamfolder = bam_dir,
      annotation = annotation,
      transcript_align = TRUE,
      name_samples = name_samples
    )
    offsets <- psite(reads, flanking = 6, plot = FALSE)
    fwrite(offsets, file.path(outdir, "offsets", paste0(sample_id, ".ribowaltz_offsets.tsv")), sep = "\t")
    status <- rbind(status, data.table(sample_id = sample_id, status = "ok", message = "offsets computed"))
    status
  }, error = function(e) {
    status <- rbind(status, data.table(sample_id = sample_id, status = "failed", message = conditionMessage(e)))
    status
  })
  status <- result
}

fwrite(status, file.path(outdir, "ribowaltz_status.tsv"), sep = "\t")
