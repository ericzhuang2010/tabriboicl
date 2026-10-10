configure_learn_bias_parameters <- function() {
  suppressPackageStartupMessages(library(RiboBA))
  ns <- asNamespace("RiboBA")
  original <- get("learn_bias_parameters", envir = ns)
  if (isTRUE(attr(original, "tabriboicl_configured"))) {
    return(invisible())
  }

  source_text <- paste(deparse(original, width.cutoff = 500), collapse = "\n")
  old_text <- paste(
    "        if (number_min < 1000) {",
    "            number_min <- 1000",
    "        }",
    sep = "\n"
  )
  new_text <- paste(
    "        min_group_count <- suppressWarnings(as.integer(Sys.getenv(\"RIBOBA_OFFSET_MIN_GROUP_COUNT\", \"1000\")))",
    "        if (is.na(min_group_count)) {",
    "            min_group_count <- 1000L",
    "        }",
    "        if (number_min < min_group_count) {",
    "            number_min <- min_group_count",
    "        }",
    sep = "\n"
  )
  configured_text <- sub(old_text, new_text, source_text, fixed = TRUE)
  if (identical(configured_text, source_text)) {
    stop("Unable to configure RiboBA learn_bias_parameters read-length threshold")
  }

  old_offset_text <- "            offset_p <- as.integer(as.numeric(names((which.max(freq_v)))))"
  new_offset_text <- paste(
    "            if ((all(is.na(freq_v)) || sum(freq_v, na.rm = TRUE) == 0L) && identical(Sys.getenv(\"RIBOBA_ALLOW_SYNTHETIC_OFFSETS\", \"0\"), \"1\")) {",
    "                frame_counts <- table(df_read$frame_idx)",
    "                if (length(frame_counts) == 0L) {",
    "                  stop(\"Cannot infer fallback RiboBA offset: no frame counts\")",
    "                }",
    "                target_frame <- as.integer(names(which.max(frame_counts)))",
    "                target_mod <- (1L - target_frame) %% 3L",
    "                candidates <- limited_range[(limited_range + 1L) %% 3L == target_mod]",
    "                if (length(candidates) == 0L) {",
    "                  candidates <- limited_range",
    "                }",
    "                offset_p <- candidates[which.min(abs(candidates - round(df_read$qwidth[1] / 3)))]",
    "            }",
    "            else {",
    "                offset_p <- as.integer(as.numeric(names((which.max(freq_v)))))",
    "            }",
    sep = "\n"
  )
  offset_configured_text <- sub(old_offset_text, new_offset_text, configured_text, fixed = TRUE)
  if (identical(offset_configured_text, configured_text)) {
    stop("Unable to configure RiboBA empirical-offset fallback")
  }
  configured_text <- offset_configured_text

  old_trunc_text <- "    trunc_cds <- IRanges::IRanges(start = tx_info$mix_tx_pos$utr5_p3 + 16L, end = tx_info$mix_tx_pos$utr3_p5 - 16L)"
  new_trunc_text <- paste(
    "    trunc_margin <- suppressWarnings(as.integer(Sys.getenv(\"RIBOBA_TRUNC_MARGIN\", \"16\")))",
    "    if (is.na(trunc_margin)) {",
    "        trunc_margin <- 16L",
    "    }",
    "    trunc_start <- tx_info$mix_tx_pos$utr5_p3 + trunc_margin",
    "    trunc_end <- tx_info$mix_tx_pos$utr3_p5 - trunc_margin",
    "    trunc_end <- pmax(trunc_start, trunc_end)",
    "    trunc_cds <- IRanges::IRanges(start = trunc_start, end = trunc_end)",
    sep = "\n"
  )
  trunc_configured_text <- sub(old_trunc_text, new_trunc_text, configured_text, fixed = TRUE)
  if (identical(trunc_configured_text, configured_text)) {
    stop("Unable to configure RiboBA CDS trimming margin")
  }
  configured_text <- trunc_configured_text

  old_empty_candidate_text <- paste(
    "        candi_cut3 <- candi_cut3[idx_read, ]",
    "        shrink_p <- as.factor(candi_psite@x)",
    sep = "\n"
  )
  new_empty_candidate_text <- paste(
    "        candi_cut3 <- candi_cut3[idx_read, ]",
    "        if (nrow(df_rpf) == 0L || length(candi_psite@x) == 0L) {",
    "            return(rep(0, length(pos_psite)))",
    "        }",
    "        shrink_p <- as.factor(candi_psite@x)",
    sep = "\n"
  )
  empty_candidate_configured_text <- sub(
    old_empty_candidate_text,
    new_empty_candidate_text,
    configured_text,
    fixed = TRUE
  )
  if (identical(empty_candidate_configured_text, configured_text)) {
    stop("Unable to configure RiboBA empty-candidate P-site handling")
  }
  configured_text <- empty_candidate_configured_text

  configured <- eval(parse(text = configured_text))
  environment(configured) <- environment(original)
  attr(configured, "tabriboicl_configured") <- TRUE

  if (bindingIsLocked("learn_bias_parameters", ns)) {
    unlockBinding("learn_bias_parameters", ns)
  }
  assign("learn_bias_parameters", configured, envir = ns)
  lockBinding("learn_bias_parameters", ns)
  invisible()
}

patch_learn_bias_parameters <- function() {
  suppressPackageStartupMessages(library(RiboBA))
  ns <- asNamespace("RiboBA")
  original <- get("learn_bias_parameters", envir = ns)
  if (isTRUE(attr(original, "tabriboicl_wrapped"))) {
    return(invisible())
  }

  use_no_add <- function(lst_bam) {
    if (is.list(lst_bam) && "no_add" %in% names(lst_bam)) {
      return(list(no_add = lst_bam$no_add))
    }
    lst_bam
  }

  replacement <- function(tx_info, lst_bam, rnase_bias, RNase, ligat_par,
                          init_num = 20, maxfeval = 5000,
                          plot_lst_rpf_dist = FALSE, plot_file = NULL) {
    run_original <- function(reads) {
      withCallingHandlers(
        original(
          tx_info = tx_info,
          lst_bam = reads,
          rnase_bias = rnase_bias,
          RNase = RNase,
          ligat_par = ligat_par,
          init_num = init_num,
          maxfeval = maxfeval,
          plot_lst_rpf_dist = plot_lst_rpf_dist,
          plot_file = plot_file
        ),
        error = function(err) {
          if (identical(Sys.getenv("RIBOBA_DEBUG_ERROR_CALLS"), "1")) {
            calls <- sys.calls()
            call_text <- vapply(
              calls,
              function(call) paste(deparse(call, width.cutoff = 160), collapse = " "),
              character(1)
            )
            message("RIBOBA error call stack:\n", paste(call_text, collapse = "\n"))
            for (frame in sys.frames()) {
              if (
                exists("prob_mp", envir = frame, inherits = FALSE) &&
                  exists("prod_hd", envir = frame, inherits = FALSE)
              ) {
                prob_mp <- get("prob_mp", envir = frame, inherits = FALSE)
                prod_hd <- get("prod_hd", envir = frame, inherits = FALSE)
                message(
                  "RIBOBA failing dimensions: prob_mp=",
                  paste(dim(prob_mp), collapse = "x"),
                  "; prod_hd_length=",
                  length(prod_hd)
                )
              }
            }
          }
        }
      )
    }

    if (identical(Sys.getenv("RIBOBA_FORCE_NO_ADD5"), "1")) {
      message("RIBOBA_FORCE_NO_ADD5=1; using no_add reads for bias learning.")
      return(run_original(use_no_add(lst_bam)))
    }

    result <- tryCatch(
      run_original(lst_bam),
      error = function(err) err
    )
    if (!inherits(result, "error")) {
      return(result)
    }

    can_retry <- identical(Sys.getenv("RIBOBA_RETRY_NO_ADD5", "1"), "1") &&
      is.list(lst_bam) &&
      "no_add" %in% names(lst_bam)
    if (!can_retry) {
      stop(result)
    }

    message(
      "learn_bias_parameters failed: ",
      conditionMessage(result),
      "; retrying without add5 reads."
    )
    run_original(use_no_add(lst_bam))
  }
  attr(replacement, "tabriboicl_wrapped") <- TRUE

  if (bindingIsLocked("learn_bias_parameters", ns)) {
    unlockBinding("learn_bias_parameters", ns)
  }
  assign("learn_bias_parameters", replacement, envir = ns)
  lockBinding("learn_bias_parameters", ns)
  invisible()
}

require_env <- function(name) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) {
    stop(name, " is required")
  }
  value
}

get_prepare_bam_inputs <- function(bam_dir) {
  bam_cds <- sort(list.files(
    bam_dir,
    pattern = "_cds_txsorted\\.bam$",
    full.names = TRUE
  ))
  bam_ncrna <- sort(list.files(
    bam_dir,
    pattern = "_lncrna_txsorted\\.bam$",
    full.names = TRUE
  ))
  bai_cds <- paste0(bam_cds, ".bai")
  bai_ncrna <- paste0(bam_ncrna, ".bai")
  if (length(bam_cds) == 0L || length(bam_ncrna) == 0L) {
    stop("No CDS/lncRNA BAM files found in bam_dir: ", bam_dir)
  }
  if (length(bam_cds) != length(bam_ncrna)) {
    stop("CDS and lncRNA BAM counts do not match in bam_dir: ", bam_dir)
  }
  if (!all(file.exists(c(bai_cds, bai_ncrna)))) {
    stop("Missing BAM index (.bai) files in bam_dir: ", bam_dir)
  }
  cds_sample <- sub("_cds_txsorted\\.bam$", "", basename(bam_cds))
  ncrna_sample <- sub("_lncrna_txsorted\\.bam$", "", basename(bam_ncrna))
  if (!identical(cds_sample, ncrna_sample)) {
    stop("Sample names do not match between CDS and lncRNA BAM files in bam_dir: ", bam_dir)
  }

  requested_sample <- Sys.getenv("RIBOBA_SAMPLE")
  if (nzchar(requested_sample)) {
    selected <- which(cds_sample == requested_sample)
    if (length(selected) != 1L) {
      stop(
        "RIBOBA_SAMPLE must match exactly one prepared sample; got ",
        requested_sample
      )
    }
    bam_cds <- bam_cds[selected]
    bai_cds <- bai_cds[selected]
    bam_ncrna <- bam_ncrna[selected]
    bai_ncrna <- bai_ncrna[selected]
    cds_sample <- cds_sample[selected]
  }

  list(
    bam_cds = bam_cds,
    bai_cds = bai_cds,
    bam_ncrna = bam_ncrna,
    bai_ncrna = bai_ncrna,
    sample_names = cds_sample
  )
}

validate_rnase <- function(rnase, n_sample) {
  allowed <- c("rnase-i", "mnase", "p1")
  rnase <- tolower(rnase)
  if (!all(rnase %in% allowed)) {
    stop("RNase must be one of: ", paste(allowed, collapse = ", "))
  }
  if (length(rnase) == 1L) {
    return(rep(rnase, n_sample))
  }
  if (length(rnase) != n_sample) {
    stop("RNase length must be 1 or equal to number of samples: ", n_sample)
  }
  rnase
}

resolve_rnase <- function(bam_dir) {
  bam_inputs <- get_prepare_bam_inputs(bam_dir)
  table_path <- Sys.getenv("RIBOBA_RNASE_TABLE")
  if (!nzchar(table_path)) {
    return(validate_rnase(require_env("RIBOBA_RNASE"), length(bam_inputs$sample_names)))
  }
  if (!file.exists(table_path)) {
    stop("RIBOBA_RNASE_TABLE does not exist: ", table_path)
  }

  rnase_table <- utils::read.delim(
    table_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c("sample_id", "rnase")
  missing <- setdiff(required, names(rnase_table))
  if (length(missing) > 0L) {
    stop(
      "RIBOBA_RNASE_TABLE is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
  if (anyDuplicated(rnase_table$sample_id)) {
    stop("RIBOBA_RNASE_TABLE contains duplicate sample_id values")
  }

  matched <- match(bam_inputs$sample_names, rnase_table$sample_id)
  if (anyNA(matched)) {
    stop(
      "RIBOBA_RNASE_TABLE has no entry for: ",
      paste(bam_inputs$sample_names[is.na(matched)], collapse = ", ")
    )
  }
  validate_rnase(rnase_table$rnase[matched], length(bam_inputs$sample_names))
}

export_tx_map <- function(transcript_file, output_path) {
  if (!nzchar(output_path)) {
    return(invisible(NULL))
  }

  load(transcript_file)
  if (!exists("tx_info_lst", inherits = FALSE)) {
    stop("tx_information.RData does not contain tx_info_lst")
  }
  tx_lens <- tx_info_lst$tx_info$tx_lens
  mix_tx_pos <- tx_info_lst$tx_info$mix_tx_pos
  if (is.null(tx_lens) || is.null(mix_tx_pos)) {
    stop("tx_info_lst$tx_info must contain tx_lens and mix_tx_pos")
  }
  if (nrow(tx_lens) != nrow(mix_tx_pos)) {
    stop("tx_lens and mix_tx_pos row counts do not match")
  }

  tx_map <- data.frame(
    tx_id = tx_lens$tx_id,
    tx_name = tx_lens$tx_name,
    gene_id = tx_lens$gene_id,
    tx_len = tx_lens$tx_len,
    cds_len = tx_lens$cds_len,
    utr5_len = tx_lens$utr5_len,
    utr3_len = tx_lens$utr3_len,
    mix_start = mix_tx_pos$utr5_p5,
    mix_cds_start = mix_tx_pos$utr5_p3 + 1L,
    mix_cds_end = mix_tx_pos$utr3_p5 - 1L,
    mix_end = mix_tx_pos$utr3_p3
  )
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    tx_map,
    file = output_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  invisible(tx_map)
}

run_native_profiles <- function(transcript_file, bam_dir, rnase,
                                 downstream_dir, export_downstream) {
  ns <- asNamespace("RiboBA")
  bam_inputs <- get_prepare_bam_inputs(bam_dir)
  rnase <- validate_rnase(rnase, length(bam_inputs$sample_names))

  message("Running native RiboBA generative profiles for ", length(bam_inputs$sample_names), " sample(s).")
  load(transcript_file)
  reads_bam <- get("read_bam_alignments", ns)(
    bam_path = bam_inputs$bam_cds,
    bai_path = bam_inputs$bai_cds,
    tx_info = tx_info_lst$tx_info,
    only_uni_map = TRUE,
    add_ratio = 0.1
  )

  result <- vector("list", length(reads_bam))
  names(result) <- bam_inputs$sample_names
  status_rows <- vector("list", length(reads_bam))
  for (i in seq_along(reads_bam)) {
    sample_name <- bam_inputs$sample_names[i]
    rnase_i <- rnase[i]
    rnase_bias <- identical(rnase_i, "mnase")
    sample_dir <- file.path(downstream_dir, sample_name)
    output_files <- file.path(
      sample_dir,
      paste0(
        sample_name,
        c(
          ".riboba_profile.tsv",
          ".riboba_psite_profile.tsv",
          ".riboba_corrected_profile.tsv"
        )
      )
    )
    if (isTRUE(export_downstream)) {
      dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
      unlink(output_files, force = TRUE)
    }

    message(
      sprintf("Native profile for sample %s (RNase=%s) ... ", sample_name, rnase_i),
      appendLF = FALSE
    )
    sample_result <- tryCatch(
      {
        par_rpf <- get("learn_bias_parameters", ns)(
          tx_info = tx_info_lst$tx_info,
          lst_bam = reads_bam[[i]],
          rnase_bias = rnase_bias,
          RNase = rnase_i,
          ligat_par = TRUE
        )
        profile <- par_rpf$lst_rpf$rpf_overlap
        corrected_profile <- par_rpf$rpf_corrected
        if (is.null(par_rpf$ribo_num)) {
          stop("Installed RiboBA does not expose learned ribo_num P-site occupancies")
        }
        positive_psites <- which(
          is.finite(par_rpf$ribo_num) & par_rpf$ribo_num > 0
        )
        psite_profile <- data.frame(
          pos = positive_psites,
          weight = as.numeric(par_rpf$ribo_num[positive_psites])
        )
        if (nrow(psite_profile) == 0L) {
          stop("Native RiboBA learned an empty P-site occupancy profile")
        }
        list(
          value = list(
            mode = "native_profile",
            bias_info = par_rpf$lig_par,
            par_0 = par_rpf$par_0,
            profile = profile,
            corrected_profile = corrected_profile,
            psite_profile = psite_profile
          ),
          profile = profile,
          corrected_profile = corrected_profile,
          psite_profile = psite_profile
        )
      },
      error = function(err) err
    )

    if (inherits(sample_result, "error")) {
      message("Failed: ", conditionMessage(sample_result))
      status_rows[[i]] <- data.frame(
        sample_id = sample_name,
        status = "failed",
        message = conditionMessage(sample_result)
      )
      next
    }

    result[[i]] <- sample_result$value
    if (isTRUE(export_downstream)) {
      utils::write.table(
        sample_result$profile,
        file = output_files[1],
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
      utils::write.table(
        sample_result$psite_profile,
        file = output_files[2],
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
      if (!is.null(sample_result$corrected_profile)) {
        utils::write.table(
          sample_result$corrected_profile,
          file = output_files[3],
          sep = "\t",
          row.names = FALSE,
          quote = FALSE
        )
      }
    }
    status_rows[[i]] <- data.frame(
      sample_id = sample_name,
      status = "ok",
      message = paste0("native RiboBA P-site profile; RNase=", rnase_i)
    )
    message("Done")
  }

  profile_status_path <- Sys.getenv("RIBOBA_PROFILE_STATUS")
  profile_status <- do.call(rbind, status_rows)
  if (nzchar(profile_status_path)) {
    dir.create(
      dirname(profile_status_path),
      recursive = TRUE,
      showWarnings = FALSE
    )
    utils::write.table(
      profile_status,
      file = profile_status_path,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }

  successful <- !vapply(result, is.null, logical(1))
  result <- result[successful]
  if (length(result) == 0L) {
    stop("Native RiboBA profile learning failed for every sample")
  }
  attr(result, "mode") <- "native_profile"
  result
}

configure_learn_bias_parameters()
patch_learn_bias_parameters()

tx_info <- require_env("RIBOBA_TX_INFO")
bam_dir <- require_env("RIBOBA_BAM_DIR")
rnase <- resolve_rnase(bam_dir)
result_path <- require_env("RIBOBA_RESULT")
tx_map_path <- Sys.getenv("RIBOBA_TX_MAP", file.path(dirname(result_path), "riboba_tx_map.tsv"))
downstream_dir <- Sys.getenv("RIBOBA_DOWNSTREAM_DIR", file.path(dirname(result_path), "downstream"))
export_downstream <- !identical(Sys.getenv("RIBOBA_EXPORT_DOWNSTREAM", "1"), "0")
min_lnc_reads <- as.integer(Sys.getenv("RIBOBA_MIN_LNC_READS", "100"))
if (is.na(min_lnc_reads)) {
  stop("RIBOBA_MIN_LNC_READS must be an integer")
}

allow_profile_fallback <- identical(Sys.getenv("RIBOBA_ALLOW_PROFILE_FALLBACK", "1"), "1")
profile_only <- !identical(Sys.getenv("RIBOBA_PROFILE_ONLY", "1"), "0")

if (export_downstream) {
  dir.create(downstream_dir, recursive = TRUE, showWarnings = FALSE)
}
export_tx_map(tx_info, tx_map_path)
if (profile_only) {
  res <- run_native_profiles(
    transcript_file = tx_info,
    bam_dir = bam_dir,
    rnase = rnase,
    downstream_dir = downstream_dir,
    export_downstream = export_downstream
  )
  saveRDS(res, file = result_path)
  quit(save = "no", status = 0L)
}


res <- tryCatch(
  RiboBA::run_riboba_pipeline(
    transcript.file = tx_info,
    bam_dir = bam_dir,
    RNase = rnase,
    min_lnc_reads_for_prediction = min_lnc_reads,
    export_downstream = export_downstream,
    downstream_output_dir = downstream_dir
  ),
  error = function(err) err
)
if (inherits(res, "error")) {
  if (!allow_profile_fallback) {
    stop(res)
  }
  full_pipeline_error <- conditionMessage(res)
  message("RiboBA full pipeline failed: ", full_pipeline_error)
  message("Falling back to RiboBA bias/profile output.")
  res <- run_native_profiles(
    transcript_file = tx_info,
    bam_dir = bam_dir,
    rnase = rnase,
    downstream_dir = downstream_dir,
    export_downstream = export_downstream
  )
  attr(res, "full_pipeline_error") <- full_pipeline_error
}
saveRDS(res, file = result_path)
