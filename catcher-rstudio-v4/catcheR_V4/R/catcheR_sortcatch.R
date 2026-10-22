#' Correct barcode swaps in the annotated expression matrix
#'
#' Applies barcode swap corrections to a CATCHER expression matrix produced
#' by [catcheR_filtercatch()] or [catcheR_nocatch()]. The matrix cell names
#' encode the barcode identity (`cellID_UMIxUCI_barcode_gene_UCI`). If a swap
#' has occurred — i.e. a UCI was initially associated with the wrong barcode —
#' this function uses a swap table (`reliable_clones_swaps_<DIs>_<ratio>.csv`
#' from [catcheR_step1QC()]) to rename affected cells to their correct
#' barcode/gene assignment.
#'
#' @param folder  Path to the working folder containing the matrix and
#'   barcodes file.
#' @param matrix  Filename of the gene expression matrix CSV (cells as
#'   columns). Typically `"silencing_matrix_all_samples.csv"` or
#'   `"silencing_matrix_complete_all_samples.csv"`.
#' @param swaps   Filename of the swap table CSV, produced by
#'   [catcheR_step1QC()]. Typically
#'   `"reliable_clones_swaps_<DIs>_<ratio>.csv"`.
#'
#' @return A data frame of the corrected expression matrix (invisibly).
#'   The corrected matrix is also written to
#'   `folder/silencing_matrix_updated.csv`.
#'
#' @section How swap correction works:
#' Each cell column in the expression matrix is named
#' `cellID_UMIxUCI_barcode_gene_UCI`. The `barcode_gene_UCI` part is the
#' "clone" identifier. The swap table from [catcheR_step1QC()] maps
#' incorrect (swapped) clone identifiers to their correct counterparts.
#' Cells whose clone is in the swap table are renamed; all others are
#' unchanged.
#'
#' @examples
#' \dontrun{
#' corrected <- catcheR_sortcatch(
#'   folder = "/path/to/experiment/",
#'   matrix = "silencing_matrix_all_samples.csv",
#'   swaps  = "reliable_clones_swaps_100_10.csv"
#' )
#'
#' # Check how many cells were corrected
#' cat("Swaps corrected:", sum(colnames(corrected) != original_names))
#' }
#'
#' @seealso [catcheR_step1QC()], [catcheR_nocatch()]
#' @export
catcheR_sortcatch <- function(folder,
                               matrix,
                               swaps,
                               output_dir  = "Output/",
                               sample_name = NULL) {

  ptm <- proc.time()

  if (!dir.exists(folder))
    stop("Folder not found: ", folder, call. = FALSE)
  folder <- normalizePath(folder, mustWork = TRUE)

  matrix_path <- file.path(folder, matrix)
  swaps_path  <- file.path(folder, swaps)
  bc_path     <- file.path(folder, "rc_barcodes_genes.csv")

  if (!file.exists(matrix_path))
    stop("Matrix file not found: ", matrix_path, call. = FALSE)
  if (!file.exists(swaps_path))
    stop("Swaps file not found: ", swaps_path,
         "\nRun catcheR_step1QC() first.", call. = FALSE)
  if (!file.exists(bc_path))
    stop("rc_barcodes_genes.csv not found in: ", folder, call. = FALSE)

  dirs <- .setup_output_dirs(output_dir, "sortcatch", sample_name)

  message("=== catcheR::catcheR_sortcatch() ===")
  message(sprintf("  Folder : %s", folder))
  message(sprintf("  Matrix : %s", matrix))
  message(sprintf("  Swaps  : %s", swaps))

  # Load
  mat            <- utils::read.csv(matrix_path, row.names=1, header=TRUE)
  swaps_tbl      <- utils::read.table(swaps_path, sep=",", header=TRUE, row.names=1)
  barcodes_genes <- utils::read.csv(bc_path, header=FALSE,
                                     col.names=c("barcode","name"))

  # Build association: swapped clone -> correct clone
  swaps_tbl <- dplyr::left_join(swaps_tbl, barcodes_genes, by=c("name"="name"))
  swaps_tbl$gene <- sub("\\..*","", swaps_tbl$name)
  swaps_tbl <- dplyr::mutate(swaps_tbl,
    clone = paste(barcode, gene, UCI, sep="_"))

  names(barcodes_genes) <- c("actual_barcode","actual_name")
  barcodes_genes <- barcodes_genes[!duplicated(barcodes_genes$actual_name),]
  swaps_tbl <- dplyr::left_join(swaps_tbl, barcodes_genes,
                                 by=c("actual_name"="actual_name"))
  swaps_tbl$actual_gene <- sub("\\..*","", swaps_tbl$actual_name)
  swaps_tbl <- dplyr::mutate(swaps_tbl,
    actual_clone = paste(actual_barcode, actual_gene, UCI, sep="_"))
  association <- swaps_tbl[, c("clone","actual_clone")]

  # Parse column names
  names_df <- as.data.frame(colnames(mat))
  names(names_df) <- "name"
  names_df <- tidyr::separate(names_df, name,
    into=c("cellID","UMIxUCI","barcode","gene","UCI"),
    sep="_", remove=FALSE)
  names_df <- dplyr::mutate(names_df, clone=paste(barcode,gene,UCI,sep="_"))
  names_df <- dplyr::left_join(names_df, association, by="clone")

  n_swaps <- sum(!is.na(names_df$actual_clone))
  message(sprintf("  Swaps corrected: %d", n_swaps))

  names_df$actual_clone <- ifelse(is.na(names_df$actual_clone),
                                   names_df$clone, names_df$actual_clone)
  names_df <- dplyr::mutate(names_df,
    new_names = paste(cellID, UMIxUCI, actual_clone, sep="_"))

  colnames(mat) <- names_df$new_names

  # Save
  .save_csv(mat, dirs, "silencing_matrix_updated")

  # Stats
  .save_stats(
    c(sprintf("Swaps corrected : %d", n_swaps),
      sprintf("Total cells     : %d", ncol(mat)),
      sprintf("Matrix dims     : %d genes x %d cells", nrow(mat), ncol(mat))),
    dirs, "sortcatch_stats"
  )

  result <- list(
    matrix      = mat,
    n_swaps     = n_swaps,
    name_map    = names_df[, c("name","new_names","clone","actual_clone")],
    paths = list(
      output     = dirs$root,
      csv        = dirs$csv,
      R_objects  = dirs$R_objects,
      to_scratch = dirs$to_scratch
    )
  )

  .save_r_objects(result, dirs, "sortcatch")
  # to_scratch: updated matrix is input for downstream single-cell analysis
  .to_scratch(file.path(dirs$csv, "silencing_matrix_updated.csv"), dirs)

  elapsed <- proc.time() - ptm
  message(sprintf("  Done in %.1f seconds. Output in: %s",
                  elapsed["elapsed"], dirs$root))
  invisible(result)
}
