# =============================================================================
# catcheR_10XcatchQC.R
# =============================================================================

#' Exploratory QC step for 10X CATCHER data
#'
#' Runs the exploratory analysis on already-extracted cell barcode files
#' (produced by the slicing step of [catcheR_10Xcatch()]). This step:
#' \itemize{
#'   \item Loads the per-cell UCI, UMI, barcode, and reference tables.
#'   \item Filters reads by the expected reference sequence.
#'   \item Deduplicates UMIs.
#'   \item Plots UMI x UCI distributions and estimates a UMI threshold
#'     (either from bimodal valley detection or noise level).
#'   \item Writes `UMI_threshold.txt` used by [catcheR_filtercatch()].
#' }
#'
#' Use this function when you want to inspect the UMI distribution and set
#' a manual threshold before running [catcheR_filtercatch()], rather than
#' running the full pipeline with [catcheR_10Xcatch()].
#'
#' @param folder    Path to the working folder.
#' @param reference Reference sequence to select reads containing the UCI.
#'   Default `"GGCGCGTTCATCTGGGGGAGCCG"`.
#' @param mode      `"bimodal"` (default) or `"noise"`. See
#'   [catcheR_10Xcatch()] for details.
#' @param samples   Number of samples. Default `1`.
#' @param x         X-axis limit for zoomed UMI x UCI plot. Default `100`.
#' @param y         Y-axis limit for zoomed UMI x UCI plot. Default `400`.
#' @param rscript   Path to Rscript executable. Default `"Rscript"`.
#'
#' @return Invisibly returns the path to `folder`.
#'
#' @section Output files written to `folder/Results_<i>/`:
#' \describe{
#'   \item{`complete_table_fin_<i>.csv/.rds`}{Full annotated table.}
#'   \item{`UMI_threshold.txt`}{Estimated UMI threshold for this sample.}
#'   \item{`barcode_distribution.pdf/.jpg`}{UMIs per shRNA barcode.}
#'   \item{`gene_distribution.pdf/.jpg`}{UMIs per gene.}
#'   \item{`UMIxUCI.pdf/.jpg`}{Full UMI x UCI distribution histogram.}
#'   \item{`UMIxUCI_<x>_<y>.pdf/.jpg`}{Zoomed UMI x UCI histogram.}
#'   \item{`UCIxcell.pdf/.jpg`}{UCIs per cell histogram.}
#'   \item{`percentage_of_UMIxUCI_dist.pdf/.jpg`}{UMI percentage distribution.}
#'   \item{`2D_percentage_of_UMIxUCI_UMI_count.pdf/.jpg`}{2D scatter.}
#' }
#'
#' @examples
#' \dontrun{
#' catcheR_10XcatchQC(
#'   folder    = "/path/to/experiment/",
#'   reference = "GGCGCGTTCATCTGGGGGAGCCG",
#'   mode      = "noise",
#'   samples   = 4
#' )
#' }
#'
#' @seealso [catcheR_10Xcatch()], [catcheR_filtercatch()]
#' @export
catcheR_10XcatchQC <- function(folder,
                                output_dir = "Output/",
                                sample_name = NULL,
                                reference = "GGCGCGTTCATCTGGGGGAGCCG",
                                mode      = c("bimodal", "noise"),
                                samples   = 1,
                                x         = 100,
                                y         = 400,
                                rscript   = "Rscript") {

  mode <- match.arg(mode)
  ptm  <- proc.time()

  if (!dir.exists(folder))
    stop("Folder not found: ", folder, call. = FALSE)
  folder <- normalizePath(folder, mustWork = TRUE)

  r_script <- system.file("scripts",
                           "barcode_silencing_explorative_analysis.R",
                           package = "catcheR")
  if (!nzchar(r_script))
    stop("Cannot find barcode_silencing_explorative_analysis.R.", call. = FALSE)

  dirs <- .setup_output_dirs(output_dir, "10XcatchQC", sample_name)

  message("=== catcheR::catcheR_10XcatchQC() ===")
  message(sprintf("  Folder: %s | Mode: %s | Samples: %d", folder, mode, samples))

  for (i in seq_len(samples)) {
    message(sprintf("  Processing sample %d / %d ...", i, samples))
    ret <- system2(rscript,
      args   = c(r_script, paste0(folder,"/"), reference, mode,
                 as.character(i), as.character(x), as.character(y)),
      stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(ret,"status")) && attr(ret,"status") != 0) {
      message(paste(ret, collapse="\n"))
      stop(sprintf("barcode_silencing_explorative_analysis.R failed on sample %d.", i),
           call. = FALSE)
    }
  }

  # Collect outputs
  rds_files <- list.files(folder, pattern="complete_table_fin_.*\\.rds",
                           full.names=TRUE, recursive=TRUE)
  thr_files <- list.files(folder, pattern="UMI_threshold\\.txt",
                           full.names=TRUE, recursive=TRUE)
  plot_files <- list.files(folder,
    pattern="barcode_distribution|gene_distribution|UMIxUCI|UCIxcell|percentage_of_UMI",
    full.names=TRUE, recursive=TRUE)
  plot_files <- plot_files[grepl("\\.(pdf|jpg)$", plot_files)]

  for (f in plot_files)
    file.copy(f, file.path(dirs$plots, basename(f)), overwrite=TRUE)

  # Read thresholds
  thresholds <- sapply(thr_files, function(f)
    as.numeric(readLines(f, warn=FALSE)[1]))

  # Load complete tables as R objects
  tables <- lapply(rds_files, readRDS)
  names(tables) <- paste0("sample_", seq_along(tables))

  result <- list(
    complete_tables = tables,
    thresholds      = thresholds,
    paths = list(
      output     = dirs$root,
      plots      = dirs$plots,
      R_objects  = dirs$R_objects,
      to_scratch = dirs$to_scratch
    )
  )

  .save_r_objects(result, dirs, "10XcatchQC")
  # to_scratch: complete tables and thresholds needed by catcheR_filtercatch
  .to_scratch(c(rds_files, thr_files), dirs)

  elapsed <- proc.time() - ptm
  message(sprintf("  Done in %.1f s. Output: %s", elapsed["elapsed"], dirs$root))
  invisible(result)
}



# =============================================================================
# catcheR_filtercatch
# =============================================================================
catcheR_filtercatch <- function(folder,
                                 expression.matrix,
                                 output_dir  = "Output/",
                                 sample_name = NULL,
                                 UMI.count   = NULL,
                                 percentage  = 15,
                                 ratio       = 5,
                                 samples     = 1,
                                 rscript     = "Rscript") {

  ptm <- proc.time()

  if (!dir.exists(folder))
    stop("Folder not found: ", folder, call. = FALSE)
  folder <- normalizePath(folder, mustWork = TRUE)

  r_script <- system.file("scripts", "barcode_silencing_cell_filtering.R",
                           package = "catcheR")
  if (!nzchar(r_script))
    stop("Cannot find barcode_silencing_cell_filtering.R.", call. = FALSE)

  dirs <- .setup_output_dirs(output_dir, "filtercatch", sample_name)

  message("=== catcheR::catcheR_filtercatch() ===")
  message(sprintf("  Folder: %s | Matrix: %s", folder, expression.matrix))

  for (i in seq_len(samples)) {
    umi_val <- UMI.count
    if (is.null(umi_val)) {
      thr_file <- file.path(folder, paste0("Results_",i), "UMI_threshold.txt")
      if (!file.exists(thr_file))
        stop("UMI_threshold.txt not found for sample ", i,
             ". Run catcheR_10XcatchQC() first.", call. = FALSE)
      umi_val <- as.numeric(readLines(thr_file, warn=FALSE)[1])
      message(sprintf("  Sample %d: auto UMI threshold = %g", i, umi_val))
    }
    ret <- system2(rscript,
      args   = c(r_script, paste0(folder,"/"),
                 as.character(percentage), as.character(umi_val),
                 as.character(ratio), expression.matrix, as.character(i)),
      stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(ret,"status")) && attr(ret,"status") != 0) {
      message(paste(ret, collapse="\n"))
      stop(sprintf("barcode_silencing_cell_filtering.R failed on sample %d.", i),
           call. = FALSE)
    }
  }

  # Collect outputs
  rds_files  <- list.files(folder, pattern="silencing_matrix_\\d+\\.rds",
                            full.names=TRUE, recursive=TRUE)
  csv_files  <- list.files(folder,
    pattern="new_cell_names|table_true_UCI|cells_with_at_least",
    full.names=TRUE, recursive=TRUE)
  csv_files  <- csv_files[grepl("\\.csv$", csv_files)]
  plot_files <- list.files(folder,
    pattern="true_UCIxcell|cellsxclone|2D_percentage.*thresh|ValidCell",
    full.names=TRUE, recursive=TRUE)
  plot_files <- plot_files[grepl("\\.(pdf|jpg)$", plot_files)]

  for (f in c(csv_files))
    file.copy(f, file.path(dirs$csv, basename(f)), overwrite=TRUE)
  for (f in plot_files)
    file.copy(f, file.path(dirs$plots, basename(f)), overwrite=TRUE)

  # Load log stats
  log_files <- list.files(folder, pattern="log_part3\\.txt",
                           full.names=TRUE, recursive=TRUE)
  log_text  <- unlist(lapply(log_files, readLines))
  if (length(log_text) > 0)
    .save_stats(log_text, dirs, "filtering_stats")

  matrices <- lapply(rds_files, readRDS)
  names(matrices) <- paste0("sample_", seq_along(matrices))

  result <- list(
    matrices = matrices,
    paths = list(
      output     = dirs$root,
      csv        = dirs$csv,
      plots      = dirs$plots,
      R_objects  = dirs$R_objects,
      to_scratch = dirs$to_scratch
    )
  )

  .save_r_objects(result, dirs, "filtercatch")
  .to_scratch(rds_files, dirs)

  elapsed <- proc.time() - ptm
  message(sprintf("  Done in %.1f s. Output: %s", elapsed["elapsed"], dirs$root))
  invisible(result)
}



# =============================================================================
# catcheR_nocatch
# =============================================================================
catcheR_nocatch <- function(folder,
                             expression.matrix,
                             output_dir     = "Output/",
                             sample_name    = NULL,
                             threshold      = NULL,
                             samples        = 1,
                             ref            = "TACGCGTTCATCTGGGGGAGCCG",
                             merge.samples  = TRUE,
                             rscript        = "Rscript") {

  ptm <- proc.time()

  if (!dir.exists(folder))
    stop("Folder not found: ", folder, call. = FALSE)
  folder <- normalizePath(folder, mustWork = TRUE)

  r_empty <- system.file("scripts", "barcode_silencing_empty_selection.R",
                          package = "catcheR")
  r_merge <- system.file("scripts", "barcode_silencing_all_samples.R",
                          package = "catcheR")
  if (!nzchar(r_empty))
    stop("Cannot find barcode_silencing_empty_selection.R.", call. = FALSE)

  dirs <- .setup_output_dirs(output_dir, "nocatch", sample_name)

  message("=== catcheR::catcheR_nocatch() ===")
  message(sprintf("  Folder: %s | Ref: %s", folder, ref))

  for (i in seq_len(samples)) {
    thr_val <- threshold
    if (is.null(thr_val)) {
      thr_file <- file.path(folder, paste0("Results_",i), "UMI_threshold.txt")
      if (!file.exists(thr_file))
        stop("UMI_threshold.txt not found for sample ", i, call. = FALSE)
      thr_val <- as.numeric(readLines(thr_file, warn=FALSE)[1])
    }
    ret <- system2(rscript,
      args   = c(r_empty, folder, expression.matrix,
                 as.character(thr_val), as.character(i), ref),
      stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(ret,"status")) && attr(ret,"status") != 0) {
      message(paste(ret, collapse="\n"))
      stop(sprintf("barcode_silencing_empty_selection.R failed on sample %d.", i),
           call. = FALSE)
    }
  }

  if (merge.samples && nzchar(r_merge)) {
    message("  Merging samples ...")
    ret_m <- system2(rscript,
      args   = c(r_merge, paste0(folder,"/"), as.character(samples), "TRUE"),
      stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(ret_m,"status")) && attr(ret_m,"status") != 0)
      warning("barcode_silencing_all_samples.R failed.", call. = FALSE)
  }

  # Collect outputs
  complete_rds <- list.files(folder,
    pattern="silencing_matrix_complete.*\\.rds",
    full.names=TRUE, recursive=TRUE)
  all_csv <- list.files(folder,
    pattern="silencing_matrix_complete_all_samples\\.csv",
    full.names=TRUE)
  for (f in all_csv)
    file.copy(f, file.path(dirs$csv, basename(f)), overwrite=TRUE)

  matrices <- lapply(complete_rds, readRDS)
  names(matrices) <- paste0("sample_", seq_along(matrices))

  result <- list(
    matrices = matrices,
    paths = list(
      output     = dirs$root,
      csv        = dirs$csv,
      R_objects  = dirs$R_objects,
      to_scratch = dirs$to_scratch
    )
  )

  .save_r_objects(result, dirs, "nocatch")
  .to_scratch(c(complete_rds, all_csv), dirs)

  elapsed <- proc.time() - ptm
  message(sprintf("  Done in %.1f s. Output: %s", elapsed["elapsed"], dirs$root))
  invisible(result)
}

