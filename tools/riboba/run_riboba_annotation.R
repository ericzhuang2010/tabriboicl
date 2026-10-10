patch_devtools_install <- function() {
  suppressPackageStartupMessages(library(devtools))
  ns <- asNamespace("devtools")
  if (bindingIsLocked("install", ns)) {
    unlockBinding("install", ns)
  }
  assign(
    "install",
    function(pkg, ..., upgrade = FALSE) {
      pkg_name <- basename(normalizePath(pkg, mustWork = FALSE))
      if (!requireNamespace(pkg_name, quietly = TRUE)) {
        utils::install.packages(
          pkg,
          repos = NULL,
          type = "source",
          dependencies = FALSE,
          quiet = TRUE
        )
      }
      invisible(pkg)
    },
    envir = ns
  )
  lockBinding("install", ns)
}

annot_dir <- Sys.getenv("RIBOBA_ANNOT_DIR")
if (!nzchar(annot_dir)) {
  stop("RIBOBA_ANNOT_DIR is required")
}
if (identical(Sys.getenv("RIBOBA_CLEAR_ANNOT_DIR"), "1")) {
  unlink(annot_dir, recursive = TRUE, force = TRUE)
}

patch_devtools_install()
suppressPackageStartupMessages(library(RiboBA))

build_annotation_index(
  output_dir = annot_dir,
  scientific_name = Sys.getenv("RIBOBA_SCIENTIFIC_NAME"),
  annotation_version = Sys.getenv("RIBOBA_ANNOTATION_VERSION"),
  genome_file = Sys.getenv("RIBOBA_GENOME_FILE"),
  gff_file = Sys.getenv("RIBOBA_GFF_FILE")
)
