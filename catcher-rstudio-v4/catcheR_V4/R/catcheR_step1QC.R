#' Run step 1 QC — intermediate plasmid library analysis
#'
#' Analyses a plasmid library where **both barcode and shRNA sequence** are
#' present in the read (intermediate form, before final cloning). This is the
#' step that detects barcode–shRNA swaps and confirms correct associations
#' between UCIs, barcodes, and expected shRNA sequences.
#'
#' Internally this function:
#' 1. Runs `plasmid_inter.sh` (bundled in the package) which extracts UMI,
#'    barcode, UCI, and shRNA sequence from the FASTQ using bash `cut`.
#' 2. Calls `plasmid_inter2.R` (bundled) which loads the extracted files,
#'    matches barcodes to the reference CSV, applies DIs and ratio filters,
#'    detects swaps, and produces QC plots.
#'
#' **Required files in `folder`:**
#' - The FASTQ file (named in `fastq.read1`)
#' - `rc_barcodes_genes.csv` — two-column CSV (barcode, shRNA name)
#' - `expected_shRNA_names.csv` — two-column CSV (actual_name, shRNA sequence)
#' - Optionally `colors.csv` — two-column CSV (shRNA name, hex colour)
#' - Optionally a clones-of-interest TXT file (one barcode_UCI per line)
#'
#' @param folder         Path to the working folder containing input files.
#' @param fastq.read1    Filename of the read 1 FASTQ (or FASTQ.gz).
#' @param DIs            Minimum number of reads for the dominant shRNA
#'   assignment to a UCI-BC to be considered reliable. Default `100`.
#' @param ratio          Minimum ratio of reads for the top vs second shRNA
#'   assignment. Default `10`.
#' @param plot.threshold Minimum number of reads per clone for the output
#'   distribution plot. Default `2000`.
#' @param clones         Filename of a TXT file listing clones of interest
#'   (barcode_UCI format, one per line). Default `NULL`.
#' @param bash           Path to bash executable. Default `"bash"`.
#' @param rscript        Path to Rscript executable. Default `"Rscript"`.
#'
#' @return Invisibly returns the path to `folder`. All output files are
#'   written there.
#'
#' @section Output files written to `folder`:
#' \describe{
#'   \item{`clone_distribution_filter_<threshold>.pdf/.jpg`}{Read distribution
#'     for clones above `plot.threshold`.}
#'   \item{`gene_distribution.pdf/.jpg`}{UMIs per gene.}
#'   \item{`name_distribution.pdf/.jpg`}{UMIs per shRNA name.}
#'   \item{`shRNA_distribution.pdf/.jpg`}{UMIs per actual shRNA.}
#'   \item{`swaps_distribution.pdf/.jpg`}{Distribution of detected barcode
#'     swaps.}
#'   \item{`reliable_clones_<DIs>_<ratio>.csv`}{All reliable UCI-BC
#'     assignments.}
#'   \item{`reliable_clones_swaps_<DIs>_<ratio>.csv`}{Swapped assignments.}
#'   \item{`reliable_clones_confirmations_<DIs>_<ratio>.csv`}{Confirmed
#'     correct assignments.}
#'   \item{`shRNA_distribution/`}{Per-shRNA swap distribution plots.}
#' }
#'
#' @examples
#' \dontrun{
#' catcheR_step1QC(
#'   folder        = "/path/to/step1/",
#'   fastq.read1   = "my_library_R1.fastq",
#'   DIs           = 100,
#'   ratio         = 10,
#'   plot.threshold = 2000,
#'   clones        = "clones_of_interest.txt"
#' )
#' }
#'
#' @seealso [catcheR_step2QC_extraction()], [catcheR_sortcatch()]
#' @export
catcheR_step1QC <- function(folder,
                             fastq.read1,
                             output_dir     = "Output/",
                             sample_name    = NULL,
                             DIs            = 100,
                             ratio          = 10,
                             plot.threshold = 2000,
                             clones         = NULL,
                             bash           = "bash",
                             rscript        = "Rscript") {

  ptm <- proc.time()

  if (!dir.exists(folder))
    stop("Folder not found: ", folder, call. = FALSE)
  folder <- normalizePath(folder, mustWork = TRUE)

  for (req in c("rc_barcodes_genes.csv", "expected_shRNA_names.csv")) {
    if (!file.exists(file.path(folder, req)))
      stop("Required file missing: ", req, call. = FALSE)
  }

  if (!is.null(clones) && !file.exists(file.path(folder, clones)))
    stop("Clones file not found: ", file.path(folder, clones), call. = FALSE)

  sh_script <- system.file("scripts", "plasmid_inter.sh",  package = "catcheR")
  r_script  <- system.file("scripts", "plasmid_inter2.R",  package = "catcheR")
  if (!nzchar(sh_script)) stop("Cannot find plasmid_inter.sh.", call. = FALSE)
  if (!nzchar(r_script))  stop("Cannot find plasmid_inter2.R.", call. = FALSE)

  # Patch Docker paths
  sh_patched <- .patch_sh_script(sh_script, folder)

  dirs        <- .setup_output_dirs(output_dir, "step1QC", sample_name)
  clones_arg  <- if (is.null(clones)) "NULL" else clones

  message("=============================================================")
  message("  catcheR::catcheR_step1QC()")
  message("=============================================================")
  message(sprintf("  Folder    : %s", folder))
  message(sprintf("  FASTQ     : %s", fastq.read1))
  message(sprintf("  DIs       : %d | Ratio: %d | Threshold: %d",
                  DIs, ratio, plot.threshold))

  # Run shell script — writes intermediate txt files to folder
  message("\n  Running plasmid_inter.sh ...")
  ret_sh <- system2(bash,
    args   = c(sh_patched, folder, fastq.read1,
               as.character(plot.threshold), as.character(DIs),
               as.character(ratio), clones_arg),
    stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(ret_sh,"status")) && attr(ret_sh,"status") != 0) {
    message(paste(ret_sh, collapse="\n"))
    stop("plasmid_inter.sh failed.", call. = FALSE)
  }

  # Run R script — reads txt files, writes CSVs and plots to folder
  message("  Running plasmid_inter2.R ...")
  ret_r <- system2(rscript,
    args   = c(r_script, folder,
               as.character(plot.threshold), as.character(DIs),
               as.character(ratio), clones_arg),
    stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(ret_r,"status")) && attr(ret_r,"status") != 0) {
    message(paste(ret_r, collapse="\n"))
    stop("plasmid_inter2.R failed.", call. = FALSE)
  }

  # ── Collect outputs and copy into structured Output/ ──────────────────────
  # CSVs produced by plasmid_inter2.R
  csv_files <- list.files(folder,
    pattern = paste0("reliable_clones.*\\.csv|table_of_int"),
    full.names = TRUE)
  for (f in csv_files)
    file.copy(f, file.path(dirs$csv, basename(f)), overwrite = TRUE)

  # Plots
  plot_files <- list.files(folder,
    pattern = "clone_distribution|gene_distribution|name_distribution|shRNA_distribution|swaps_distribution",
    full.names = TRUE)
  plot_files <- plot_files[grepl("\\.(pdf|jpg)$", plot_files)]
  for (f in plot_files)
    file.copy(f, file.path(dirs$plots, basename(f)), overwrite = TRUE)

  # Load main outputs as R objects for return
  swaps_file <- list.files(folder,
    pattern = paste0("reliable_clones_swaps_", DIs, "_", ratio, "\\.csv"),
    full.names = TRUE)
  correct_file <- list.files(folder,
    pattern = paste0("reliable_clones_confirmations_", DIs, "_", ratio, "\\.csv"),
    full.names = TRUE)
  reliable_file <- list.files(folder,
    pattern = paste0("reliable_clones_", DIs, "_", ratio, "\\.csv"),
    full.names = TRUE)

  load_if_exists <- function(f) {
    if (length(f) > 0 && file.exists(f[1]))
      utils::read.csv(f[1], stringsAsFactors = FALSE)
    else NULL
  }

  swaps    <- load_if_exists(swaps_file)
  correct  <- load_if_exists(correct_file)
  reliable <- load_if_exists(reliable_file)

  if (!is.null(swaps))    .save_csv(swaps,    dirs, paste0("swaps_DIs",DIs,"_ratio",ratio))
  if (!is.null(correct))  .save_csv(correct,  dirs, paste0("correct_DIs",DIs,"_ratio",ratio))
  if (!is.null(reliable)) .save_csv(reliable, dirs, paste0("reliable_DIs",DIs,"_ratio",ratio))

  # ── Assemble result ────────────────────────────────────────────────────────
  result <- list(
    swaps    = swaps,
    correct  = correct,
    reliable = reliable,
    paths = list(
      input      = folder,
      output     = dirs$root,
      csv        = dirs$csv,
      plots      = dirs$plots,
      R_objects  = dirs$R_objects,
      to_scratch = dirs$to_scratch
    )
  )

  .save_r_objects(result, dirs, "step1QC")

  # to_scratch: swap table is needed by catcheR_sortcatch()
  if (!is.null(swaps_file) && length(swaps_file) > 0)
    .to_scratch(swaps_file, dirs)

  elapsed <- proc.time() - ptm
  message(sprintf("\n  Done in %.1f seconds.", elapsed["elapsed"]))
  message(sprintf("  Output in: %s", dirs$root))
  message("=============================================================")

  invisible(result)
}
