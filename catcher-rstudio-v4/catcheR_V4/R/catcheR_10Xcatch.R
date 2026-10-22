#' Run the full 10X CATCHER pipeline
#'
#' Complete data analysis pipeline for iPS2-seq / CATCHER 10X experiments.
#' Processes paired-end FASTQ files from a 10X Genomics experiment together
#' with a gene expression matrix, extracts UCI barcodes from read 2, assigns
#' them to cells, filters valid cells, identifies empty cells, and produces
#' an annotated expression matrix where each cell name encodes its shRNA
#' perturbation identity.
#'
#' This function runs the full pipeline in one call by executing
#' `barcode_silencing_slicing.sh` (bundled in the package), which internally
#' calls the explorative analysis, cell filtering, empty selection, and
#' multi-sample merging steps in sequence.
#'
#' For a step-by-step approach with custom thresholds between steps, use
#' [catcheR_10XcatchQC()], [catcheR_filtercatch()], and [catcheR_nocatch()]
#' individually.
#'
#' @param folder           Path to the working folder containing input files.
#' @param fastq.read1      Filename of read 1 FASTQ (cell barcode + UMI).
#'   Can be `.fastq` or `.fastq.gz`.
#' @param fastq.read2      Filename of read 2 FASTQ (shRNA barcode + UCI).
#'   Can be `.fastq` or `.fastq.gz`.
#' @param expression.matrix Filename of the gene expression matrix CSV
#'   (cells as columns, genes as rows, as output by Cell Ranger).
#' @param reference        Reverse complement of the sequence flanking the
#'   UCI on read 2. Default `"GGCGCGTTCATCTGGGGGAGCCG"`.
#' @param UCI.length       Length of the UCI sequence in nucleotides.
#'   Default `6`.
#' @param threads          Number of parallel threads for cell barcode
#'   matching. Default `2`.
#' @param percentage       Minimum percentage of UMIs supporting a UCI over
#'   total cell UMIs to consider the UCI valid. Default `15`.
#' @param mode             Method for UMI threshold estimation.
#'   `"bimodal"` (default) sets the threshold at the valley of the
#'   UMI x UCI distribution. `"noise"` sets it at
#'   1.35 × (UCIs supported by a single UMI).
#' @param ratio            Minimum ratio of UMIs for the top UCI vs second
#'   UCI to accept a cell as single-integration. Default `5`.
#' @param samples          Number of samples in the aggregated matrix
#'   (e.g. from Cell Ranger `aggr`). Default `1`.
#' @param x                X-axis limit for the zoomed UMI x UCI plot.
#'   Default `100`.
#' @param y                Y-axis limit for the zoomed UMI x UCI plot.
#'   Default `400`.
#' @param bash             Path to bash executable. Default `"bash"`.
#'
#' @return Invisibly returns the path to `folder`. All output files are
#'   written there.
#'
#' @section Required files in `folder`:
#' \describe{
#'   \item{`fastq.read1`, `fastq.read2`}{Paired FASTQ files.}
#'   \item{`expression.matrix`}{Cell Ranger gene expression CSV.}
#'   \item{`rc_barcodes_genes.csv`}{Two-column CSV: barcode, shRNA name.}
#'   \item{`grep.sh`}{Helper shell script — bundled in the package.}
#' }
#'
#' @section Output files written to `folder`:
#' Per sample `i` in `Results_i/`:
#' \describe{
#'   \item{`complete_table_fin_<i>.csv/.rds`}{Full annotated cell table.}
#'   \item{`silencing_matrix_<i>.csv/.rds`}{Valid cell expression matrix.}
#'   \item{`silencing_matrix_empty_<i>.csv/.rds`}{Empty cell matrix.}
#'   \item{`silencing_matrix_complete_<i>.csv/.rds`}{Valid + empty combined.}
#'   \item{`UMI_threshold.txt`}{Estimated UMI threshold.}
#'   \item{`new_cell_names.csv`}{Cell name mapping table.}
#'   \item{QC plots}{`barcode_distribution`, `gene_distribution`,
#'     `UMIxUCI`, `UCIxcell`, `2D_percentage_of_UMIxUCI`, etc.}
#' }
#' Combined across all samples:
#' \describe{
#'   \item{`silencing_matrix_all_samples.csv/.rds`}{Merged valid matrix.}
#'   \item{`silencing_matrix_complete_all_samples.csv/.rds`}{Merged complete.}
#' }
#'
#' @examples
#' \dontrun{
#' catcheR_10Xcatch(
#'   folder            = "/path/to/experiment/",
#'   fastq.read1       = "sample_R1.fastq.gz",
#'   fastq.read2       = "sample_R2.fastq.gz",
#'   expression.matrix = "matrix.csv",
#'   reference         = "GGCGCGTTCATCTGGGGGAGCCG",
#'   UCI.length        = 6,
#'   threads           = 12,
#'   percentage        = 15,
#'   mode              = "noise",
#'   ratio             = 5,
#'   samples           = 4
#' )
#' }
#'
#' @seealso [catcheR_10XcatchQC()], [catcheR_filtercatch()],
#'   [catcheR_nocatch()], [catcheR_sortcatch()]
#' @export
catcheR_10Xcatch <- function(folder,
                              fastq.read1,
                              fastq.read2,
                              expression.matrix,
                              reference  = "GGCGCGTTCATCTGGGGGAGCCG",
                              UCI.length = 6,
                              threads    = 2,
                              percentage = 15,
                              mode       = c("bimodal", "noise"),
                              ratio      = 5,
                              samples    = 1,
                              x          = 100,
                              y          = 400,
                              bash       = "bash") {

  mode <- match.arg(mode)
  ptm  <- proc.time()

  # ── Validate ─────────────────────────────────────────────────────────────────
  if (!dir.exists(folder))
    stop("Folder not found: ", folder, call. = FALSE)
  folder <- normalizePath(folder, mustWork = TRUE)

  for (f in c(fastq.read1, fastq.read2, expression.matrix,
              "rc_barcodes_genes.csv")) {
    if (!file.exists(file.path(folder, f)))
      stop("Required file missing from folder: ", f, call. = FALSE)
  }

  # ── Locate bundled scripts ────────────────────────────────────────────────────
  sh  <- system.file("scripts", "barcode_silencing_slicing.sh",
                     package = "catcheR")
  gsh <- system.file("scripts", "grep.sh", package = "catcheR")

  if (!nzchar(sh))
    stop("Cannot find barcode_silencing_slicing.sh inside the catcheR package.",
         call. = FALSE)

  Sys.chmod(sh,  mode = "0755")
  if (nzchar(gsh)) Sys.chmod(gsh, mode = "0755")

  # ── Patch the shell script to use bundled R scripts ──────────────────────────
  sh_patched <- .patch_sh_script(sh, folder)

  message("=============================================================")
  message("  catcheR::catcheR_10Xcatch()")
  message("=============================================================")
  message(sprintf("  Folder     : %s", folder))
  message(sprintf("  Read 1     : %s", fastq.read1))
  message(sprintf("  Read 2     : %s", fastq.read2))
  message(sprintf("  Matrix     : %s", expression.matrix))
  message(sprintf("  Reference  : %s", reference))
  message(sprintf("  UCI length : %d", UCI.length))
  message(sprintf("  Threads    : %d", threads))
  message(sprintf("  Percentage : %d", percentage))
  message(sprintf("  Mode       : %s", mode))
  message(sprintf("  Ratio      : %d", ratio))
  message(sprintf("  Samples    : %d", samples))

  ret <- system2(
    bash,
    args = c(sh_patched, folder, fastq.read1, fastq.read2,
             expression.matrix, reference,
             as.character(UCI.length), as.character(threads),
             as.character(percentage), mode,
             as.character(samples), as.character(ratio),
             as.character(x), as.character(y)),
    stdout = TRUE, stderr = TRUE
  )

  if (!is.null(attr(ret, "status")) && attr(ret, "status") != 0) {
    message(paste(ret, collapse = "\n"))
    stop("barcode_silencing_slicing.sh failed. See messages above.",
         call. = FALSE)
  }

  elapsed <- proc.time() - ptm
  message(sprintf("\n  Done in %.1f seconds.", elapsed["elapsed"]))
  message("=============================================================")

  invisible(folder)
}
