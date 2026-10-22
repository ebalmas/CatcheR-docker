#' Run the CATCHER plasmid extraction pipeline (Python)
#'
#' Calls `plasmid_final_corrected.py` — the corrected Python script that
#' extracts UMI, barcode (BC), and UCI from every read in a FASTQ file,
#' matches barcodes against a reference CSV, counts reads per clone, and
#' writes all output files to a results folder.
#'
#' **Read structure (1-based positions, corrected from original catcheR):**
#' * pos 1–13 : UMI (13 nt)
#' * pos 14–29: RED primer (16 nt)
#' * pos 30–33: GCGT overhang (4 nt) — old code started here (wrong)
#' * pos 34–41: Barcode (8 nt) — correct start
#' * pos 42–47: UCI (6 nt) — correct
#'
#' The Python script is bundled inside the package at
#' `inst/python/plasmid_final_corrected.py` and does not need to be installed
#' separately.
#'
#' @param fastq        Path to the FASTQ file (`.fastq` or `.fastq.gz`).
#' @param barcodes     Path to the barcode reference CSV.
#'   Two-column, no header: `barcode,shRNA_name`.
#' @param output_dir   Directory for output files. Created if it does not
#'   exist. Default `"results/"`.
#' @param DIs          Minimum read count per clone to be considered genuine
#'   (Distinct Integrations threshold). Default `300`.
#' @param clones       Optional path to a file of clones of interest
#'   (barcode_UCI format, one per line). Default `NULL`.
#' @param python       Path to the Python 3 executable. Default `"python3"`.
#' @param bc_start     Barcode start position (1-based). Default `34`.
#' @param bc_end       Barcode end position (1-based). Default `41`.
#' @param uci_start    UCI start position (1-based). Default `42`.
#' @param uci_end      UCI end position (1-based). Default `47`.
#' @param umi_len      UMI length in nt. Default `13`.
#'
#' @return Invisibly returns the path to the output directory.
#'
#' @section Output files written to `output_dir`:
#' \describe{
#'   \item{`distribution_all_clones.csv`}{One row per unique clone (BC_UCI).
#'     Primary input for all downstream R functions.}
#'   \item{`percentages.csv`}{Clones above DIs with read percentage.}
#'   \item{`counts_per_barcode.tsv`}{Reads and unique UCIs per shRNA barcode.}
#'   \item{`counts_per_gene.tsv`}{Reads and unique UCIs per gene.}
#'   \item{`complete_table.csv`}{All reads with UMI, BC, UCI, name, gene.}
#'   \item{`final_BC.txt`, `final_UCI.txt`, `final_UMI.txt`}{Flat files,
#'     one identifier per read.}
#'   \item{`qc_summary.txt`}{Mapping statistics.}
#' }
#'
#' @examples
#' \dontrun{
#' catcheR_step2QC_extraction(
#'   fastq      = "EB003CS1_KS7_1_finalQC_S1_L001_R1_001.fastq",
#'   barcodes   = "rc_barcodes_genes.csv",
#'   output_dir = "results/",
#'   DIs        = 300
#' )
#' }
#'
#' @export
catcheR_step2QC_extraction <- function(fastq,
                            barcodes,
                            output_dir   = "results/",
                            DIs          = 300,
                            clones       = NULL,
                            python       = "python3",
                            bc_start     = 34,
                            bc_end       = 41,
                            uci_start    = 42,
                            uci_end      = 47,
                            umi_len      = 13) {

  # Locate bundled Python script
  py_script <- system.file("python", "plasmid_final_corrected.py",
                            package = "catcheR")
  if (!nzchar(py_script)) {
    stop("Cannot find plasmid_final_corrected.py inside the catcheR package.",
         call. = FALSE)
  }

  # Validate inputs
  if (!file.exists(fastq))    stop("FASTQ file not found: ",    fastq,    call.=FALSE)
  if (!file.exists(barcodes)) stop("Barcodes file not found: ", barcodes, call.=FALSE)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Build command
  args <- c(
    py_script,
    "--fastq",     fastq,
    "--barcodes",  barcodes,
    "--output",    output_dir,
    "--DIs",       as.character(DIs),
    "--bc_start",  as.character(bc_start),
    "--bc_end",    as.character(bc_end),
    "--uci_start", as.character(uci_start),
    "--uci_end",   as.character(uci_end),
    "--umi_len",   as.character(umi_len)
  )
  if (!is.null(clones)) args <- c(args, "--clones", clones)

  message("=============================================================")
  message("  catcheR_step2QC_extraction()")
  message("=============================================================")
  message("  FASTQ    : ", fastq)
  message("  Barcodes : ", barcodes)
  message("  Output   : ", output_dir)
  message("  DIs      : ", DIs)
  message("  Running Python script...")

  ret <- system2(python, args = args)
  if (ret != 0) {
    stop("Python script exited with code ", ret,
         ".\nCheck that Python 3 is installed and the FASTQ path is correct.",
         call. = FALSE)
  }

  message("  Extraction complete. Results in: ", output_dir)
  invisible(output_dir)
}
