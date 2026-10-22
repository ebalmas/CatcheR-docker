## Part of catcheR. See R/catcheR_pipcatch.R for the entry point and
## inst/rmd/catcheR_pipcatch_report.Rmd for the report template.

# ── 2. catcheR_pip_filtercatch ──────────────────────────────────────────────

#' Re-filter a pipcatch result with refined thresholds
#'
#' @param jitter geometry for the 2D panels; see [catcheR_pipcatch()].
#' Mirrors `catcheR_filtercatch()`. Also writes the threshold-sweep table so the
#' cost of each choice is explicit, and `new_cell_names.csv` in catcheR's
#' `cellID_UMIxUCI_barcode_gene_UCI` format.
#' @export
catcheR_pip_filtercatch <- function(x,
                                    UMI.count   = 0,
                                    percentage  = 15,
                                    ratio       = 2,
                                    output_dir  = NULL,
                                    subdir      = TRUE,
                                    to.scratch  = NULL,
                                    sweep       = c(0, 1),
                                    jitter      = "auto",
                                    plots       = TRUE,
                                    verbose     = TRUE) {
  if (is.character(x)) x <- readRDS(x)
  stopifnot(inherits(x, "pipcatch"))
  if (is.null(output_dir)) output_dir <- x$output_dir
  # Same working-step logic as the other chemistries; PIPseq has no cell floor.
  if (isTRUE(subdir))
    output_dir <- file.path(output_dir,
                            sprintf("filtercatch_p%s_u%s_r%s", percentage, UMI.count, ratio))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  op <- function(f) file.path(output_dir, f)

  val <- .cat_validate(x$tab, UMI.count, percentage, ratio)
  ctf <- val$pairs %>% left_join(x$cu %>% select(cellID, UCI, barcode, shRNA, set, gene),
                                by = c("cellID", "UCI"))
  utils::write.csv(ctf, op(sprintf("complete_table_fin_%s.csv", x$sample)), row.names = FALSE)

  valid <- val$cells %>% filter(ValidCell) %>% pull(cellID)
  utils::write.csv(data.frame(cellID = val$cells$cellID[val$cells$True_UCIxcellID > 0]),
                   op("cells_with_at_least_1_UCI.csv"), row.names = FALSE)
  new_names <- ctf %>% filter(ValidCell, TorF) %>%
    transmute(cellID, Name = paste(cellID, UMIxUCI, barcode, gene, UCI, sep = "_")) %>% distinct()
  utils::write.csv(new_names, op("new_cell_names.csv"), row.names = FALSE)
  # Files the next step consumes go to scratch, renamed with chemistry and
  # sample. Opt-in: nothing is copied unless a path is given.
  if (!is.null(to.scratch))
    .cat_to_scratch(c("cells_with_at_least_1_UCI.csv" = "cells_with_at_least_1_UCI",
                      "new_cell_names.csv"            = "new_cell_names",
                      setNames("complete_table_fin",
                               sprintf("complete_table_fin_%s.csv", x$sample))),
                    from = output_dir, to = to.scratch,
                    chem = "pipcatch", sample = x$sample, verbose = verbose)
  # The input manifest is written by the EXTRACTION step, so it sits one level
  # up from this filter's own folder - copy it from there, not from here.
  if (!is.null(to.scratch))
    .cat_to_scratch(stats::setNames("input_files",
                                    sprintf("%s_%s_input_files.csv", x$sample, "pipcatch")),
                    from = x$output_dir, to = to.scratch,
                    chem = "pipcatch", sample = x$sample, verbose = verbose)

  thr_grid <- sort(unique(c(UMI.count, sweep, x$thresholds$bimodal)))
  thr_grid <- thr_grid[!is.na(thr_grid)]
  filter_returns <- bind_rows(lapply(thr_grid, function(t) {
    r <- .cat_validate(x$tab, t, percentage, ratio)
    data.frame(threshold = t, keeps_pairs_with_min_reads = t + 1L,
                   TorF_pairs = sum(r$pairs$TorF),
                   cells_with_true_UCI = sum(r$cells$True_UCIxcellID > 0),
                   ValidCell = sum(r$cells$ValidCell))
  }))
  utils::write.csv(filter_returns, op("pipcatch_filter_returns.csv"), row.names = FALSE)

  if (plots) {
    pp <- modifyList(x$params, list(UMI.count = UMI.count, percentage = percentage, ratio = ratio))
    # the 6 catcheR cell-filtering plots, recomputed at the refined thresholds
    .cat_plots_stage2(ctf, val$cells, percentage, UMI.count, output_dir, jitter)
    .cat_plots_thresholds(x$tab, x$thresholds, ctf, pp, output_dir,
                          file = "pipcatch_filtercatch_thresholds.png")
  }
  out <- list(sample = x$sample, output_dir = output_dir,
              params = list(UMI.count = UMI.count, percentage = percentage, ratio = ratio),
              complete_table_fin = ctf, cells = val$cells,
              valid_cells = valid, new_names = new_names, filter_returns = filter_returns)
  class(out) <- c("pip_filtercatch", "list")
  saveRDS(out, op(sprintf("pip_filtercatch_%s.rds", x$sample)))
  if (verbose) {
    message(sprintf("UMI.count=%s percentage=%s ratio=%s", UMI.count, percentage, ratio))
    message(sprintf("  cells with >=1 true UCI : %s", format(sum(val$cells$True_UCIxcellID > 0), big.mark = ",")))
    message(sprintf("  ValidCell               : %s", format(length(valid), big.mark = ",")))
    print(as.data.frame(filter_returns))
  }
  invisible(out)
}

# ── 3. catcheR_pip_nocatch ──────────────────────────────────────────────────

#' Identify empty (no-insertion) cells
#'
#' Mirrors `catcheR_nocatch()`. A cell is empty when it carries reads with the
#' no-insertion reference above `threshold` AND is absent from the set of cells
#' with at least one UCI. Note the empty reference is an argument here, not a
#' string buried in a shell script.
#' @export
catcheR_pip_nocatch <- function(x,
                                filtered   = NULL,
                                threshold  = 0,
                                output_dir = NULL,
                                verbose    = TRUE) {
  if (is.character(x)) x <- readRDS(x)
  stopifnot(inherits(x, "pipcatch"))
  if (is.character(filtered)) filtered <- readRDS(filtered)
  if (is.null(output_dir)) output_dir <- x$output_dir
  op <- function(f) file.path(output_dir, f)

  emp <- x$empty_reads
  if (!nrow(emp)) {
    warning("No reads carry the empty reference ", x$params$empty.reference,
            " - empty cells cannot be called at this depth.", call. = FALSE)
    return(invisible(NULL))
  }
  per_cell_empty <- emp %>% count(cellID, name = "reads_empty") %>% filter(reads_empty > threshold)
  with_uci <- if (!is.null(filtered)) filtered$cells$cellID[filtered$cells$True_UCIxcellID > 0]
              else unique(x$assignments$cellID)
  empty_cells <- per_cell_empty %>% filter(!cellID %in% with_uci) %>%
    mutate(Name = paste(cellID, "?", "empty", "empty", "empty", sep = "_"))
  utils::write.csv(empty_cells, op(sprintf("complete_table_empty_%s.csv", x$sample)), row.names = FALSE)

  if (verbose) {
    message(sprintf("Reads with the empty reference      : %s", format(nrow(emp), big.mark = ",")))
    message(sprintf("Cells above threshold %s             : %s", threshold,
                    format(nrow(per_cell_empty), big.mark = ",")))
    message(sprintf("Empty cells (and no UCI)            : %s", format(nrow(empty_cells), big.mark = ",")))
    if (nrow(empty_cells) < 50)
      message("  -> too few to serve as a negative-control population; at this depth a cell\n",
              "     with no call is far more likely to have lost its construct read to\n",
              "     enrichment-primer mispriming than to be genuinely untransduced.")
  }
  out <- list(sample = x$sample, output_dir = output_dir, threshold = threshold,
              empty_reads = emp, empty_cells = empty_cells)
  class(out) <- c("pip_nocatch", "list")
  # the prefix follows the class of `x`, so a 10x or Parse object run through
  # this shared implementation does not get a file called "pip_*"
  pfx <- if (inherits(x, "tenxcatch")) "tenxcatch"
         else if (inherits(x, "parsecatch")) "parsecatch" else "pipcatch"
  saveRDS(out, op(sprintf("%s_nocatch_%s.rds", pfx, x$sample)))
  invisible(out)
}

# ── 4. catcheR_pip_seurat ───────────────────────────────────────────────────

#' Join pipcatch calls to the PIPseeker expression matrix
#'
#' PIPseeker writes `barcodes.tsv.gz` as the **28 nt concatenated 4-tier
#' barcode** - the same string `catcheR_pipcatch()` extracts - so no barcode
#' translation is needed and the join is a plain match. (Verify once per
#' chemistry with `nchar(readLines(gzfile(barcodes)))[1]`.)
#'
#' @param mtx,features,cells  the three PIPseeker files
#' @param protein.coding optional character vector of gene symbols to keep
#' @section Which cells end up in the object:
#' The expression matrix and catcheR do not cover the same cells. The matrix
#' holds every cell the chemistry's own pipeline called; catcheR holds only
#' those whose cassette was sequenced deeply enough to survive the filter. The
#' returned object keeps **all** matrix cells by default and marks each one in
#' `catcheR_status`:
#' \itemize{
#'   \item `valid` - passed the filter (`ValidCell`), one dominant clone;
#'   \item `ambiguous` - a call, but more than one UCI clears the thresholds;
#'   \item `called_not_valid` - a call that failed the ratio test;
#'   \item `empty` - empty-vector cell, if `empty` was supplied;
#'   \item `not_detected` - no cassette read reached this cell.
#' }
#' Pass `subset.valid = TRUE` to keep only `valid` (plus `empty`) cells, which
#' is the right starting point for a silencing comparison; the function then
#' reports how many cells were dropped and in which status.
#'
#' @param mtx either a Cell Ranger `filtered_feature_bc_matrix.h5` (in which
#'   case `features` and `cells` are ignored) or the `matrix.mtx.gz` of a
#'   MatrixMarket triple.
#' @param features,cells the other two files of the triple; leave `NULL` with
#'   a `.h5`.
#' @param gem.well required only for a `cellranger aggr` matrix, which pools
#'   several GEM wells as `-1`/`-2`/...: give the suffix of the library this
#'   run corresponds to (e.g. `1`). The matrix is then restricted to that well.
#'   A matrix with a single suffix needs nothing.
#' @param rename  if TRUE, rename cells to catcheR's
#'        `cellID_UMIxUCI_barcode_gene_UCI` (and `cellID_?_empty_empty_empty`)
#' @param subset.valid keep only ValidCell (and, if supplied, empty) cells
#' @return a Seurat object if Seurat is installed, else the counts matrix with
#'         the metadata attached as an attribute. The metadata CSV is always
#'         written.
#' @export
catcheR_pip_seurat <- function(x,
                               filtered       = NULL,
                               empty          = NULL,
                               mtx, features = NULL, cells = NULL,
                               project        = "pipcatch",
                               sample_id      = NULL,
                               protein.coding = NULL,
                               feature.column = 2,
                               gem.well       = NULL,
                               rename         = FALSE,
                               subset.valid   = FALSE,
                               write.matrix.csv = FALSE,
                               output_dir     = NULL,
                               verbose        = TRUE) {
  if (is.character(x)) x <- readRDS(x)
  stopifnot(inherits(x, "pipcatch"))
  if (is.character(filtered)) filtered <- readRDS(filtered)
  if (is.character(empty))    empty    <- readRDS(empty)
  if (is.null(output_dir)) output_dir <- x$output_dir
  if (is.null(sample_id))  sample_id  <- x$sample
  # output files follow the class of `x`, not the file this code lives in
  pfx <- if (inherits(x, "tenxcatch")) "tenxcatch"
         else if (inherits(x, "parsecatch")) "parsecatch" else "pipcatch"
  op <- function(f) file.path(output_dir, f)

  # A Cell Ranger .h5 carries the barcodes, the features and the matrix in one
  # file, so `features` and `cells` are not needed when `mtx` points at one.
  is_h5 <- grepl("\\.h5$", mtx, ignore.case = TRUE)
  counts <- NULL
  if (is_h5) {
    if (!requireNamespace("Seurat", quietly = TRUE))
      stop("Reading a .h5 needs Seurat (Read10X_h5), which also needs hdf5r. ",
           "Pass the matrix.mtx.gz / features.tsv.gz / barcodes.tsv.gz triple instead.",
           call. = FALSE)
    counts <- Seurat::Read10X_h5(mtx)
    # a multimodal .h5 (e.g. with an Antibody Capture assay) returns a list
    if (is.list(counts)) {
      if (verbose)
        message("  multimodal .h5: using assay '", names(counts)[1], "' of ",
                paste(names(counts), collapse = ", "))
      counts <- counts[[1]]
    }
    gex <- colnames(counts)
  } else {
    gex <- readLines(gzfile(cells))
  }
  # Cell Ranger appends a GEM-well suffix ("-1") that the extracted cell ID
  # does not carry; left in place it makes EVERY barcode fail to join and the
  # function reports 0 calls on a perfectly good matrix. Harmless elsewhere -
  # PIPseeker and split-pipe IDs have no suffix to strip.
  # `cellranger count` and each per-sample output of `cellranger multi` carry a
  # single suffix ("-1"), which is safe to drop. `cellranger aggr` pools several
  # GEM wells as -1/-2/..., where the SAME 16 nt barcode can appear once per
  # well: dropping the suffix there would merge cells from different wells.
  # Since catcheR's cell IDs come from one enrichment library, i.e. one well,
  # the resolution is to say which well that is and keep only those cells.
  sfx   <- ifelse(grepl("-[0-9]+$", gex), sub("^.*-", "", gex), NA_character_)
  wells <- sort(unique(stats::na.omit(sfx)))
  keep_cells <- rep(TRUE, length(gex))
  if (length(wells) > 1L) {
    if (is.null(gem.well))
      stop("This matrix pools ", length(wells), " GEM wells (suffixes -",
           paste(wells, collapse = "/-"), "), as `cellranger aggr` produces. ",
           "catcheR's cell IDs come from one enrichment library, so pass ",
           "`gem.well=` with the suffix of the library this run corresponds to ",
           "(e.g. gem.well = 1). To annotate a whole aggregated object, run this ",
           "once per well and merge the objects afterwards - attaching the ",
           "metadata BEFORE the merge, since merge() re-prefixes the barcodes.",
           call. = FALSE)
    gem.well <- as.character(gem.well)
    if (!gem.well %in% wells)
      stop("gem.well = ", gem.well, " is not among the suffixes present (-",
           paste(wells, collapse = "/-"), ").", call. = FALSE)
    keep_cells <- !is.na(sfx) & sfx == gem.well
    if (verbose)
      message(sprintf("  aggregated matrix: keeping GEM well -%s (%s of %s cells; %s belong to other wells)",
                      gem.well, format(sum(keep_cells), big.mark = ","),
                      format(length(gex), big.mark = ","),
                      format(sum(!keep_cells), big.mark = ",")))
    gex <- gex[keep_cells]
  } else if (length(wells) == 1L && verbose) {
    message(sprintf("  single GEM well (-%s): dropping the suffix from %s barcodes",
                    wells, format(sum(!is.na(sfx)), big.mark = ",")))
  }
  gex <- sub("-[0-9]+$", "", gex)
  # the expected width is chemistry-specific: 16 nt for 10x, 28 nt for the
  # PIPseq 4-tier concatenation, and a well-triplet string for Parse
  expect <- if (inherits(x, "tenxcatch")) 16L else if (inherits(x, "parsecatch")) NA_integer_ else 28L
  if (!is.na(expect) && !all(nchar(gex) == expect))
    warning(sprintf("Expected %d nt %s barcodes; got widths %s. Check that this matrix ",
                    expect, if (inherits(x, "tenxcatch")) "Cell Ranger" else "PIPseeker",
                    paste(unique(nchar(gex)), collapse = "/")),
            "comes from the same run, and for PIPseeker that it is not a translated ",
            "16 nt ID (the mapping is in metrics/barcodes/barcode_whitelist.txt).",
            call. = FALSE)

  ## per-cell metadata, one row per GEX cell -----------------------------------
  pc <- x$per_cell
  if (!is.null(filtered))
    pc <- pc %>% select(-any_of(c("True_UCIxcellID", "ValidCell"))) %>%
      left_join(filtered$cells %>% select(cellID, True_UCIxcellID, ValidCell), by = "cellID")
  meta <- data.frame(cellID = gex) %>%
    left_join(pc %>% select(cellID, shRNA = top_shRNA, gene, set = top_set,
                           UCI = top_UCI, barcode = top_barcode,
                           UMIxUCI = reads_top, reads_total, n_UCI, n_sets,
                           True_UCIxcellID, ValidCell), by = "cellID") %>%
    mutate(ValidCell = tidyr::replace_na(ValidCell, FALSE))
  if (!is.null(empty) && nrow(empty$empty_cells))
    meta <- meta %>% mutate(is_empty = cellID %in% empty$empty_cells$cellID)
  else meta$is_empty <- FALSE
  meta <- meta %>% mutate(catcheR_status = case_when(
    ValidCell                      ~ "valid",
    !is.na(shRNA) & n_UCI > 1      ~ "ambiguous",
    !is.na(shRNA)                  ~ "called_not_valid",
    is_empty                       ~ "empty",
    TRUE                           ~ "not_detected"))
  utils::write.csv(meta, op(sprintf("%s_cell_metadata_%s.csv", pfx, sample_id)), row.names = FALSE)

  if (verbose) {
    message(sprintf("GEX called cells: %s", format(length(gex), big.mark = ",")))
    message(sprintf("  with a catcheR call  : %s (%.1f%%)",
                    format(sum(!is.na(meta$shRNA)), big.mark = ","), 100 * mean(!is.na(meta$shRNA))))
    message(sprintf("  ValidCell            : %s (%.1f%%)",
                    format(sum(meta$ValidCell), big.mark = ","), 100 * mean(meta$ValidCell)))
    lost <- setdiff(x$per_cell$cellID, gex)
    message(sprintf("  catcheR cells NOT in the GEX list: %s (ambient / not called as cells)",
                    format(length(lost), big.mark = ",")))
    print(meta %>% count(catcheR_status, name = "cells") %>% arrange(desc(cells)) %>% as.data.frame())
  }

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    message("Seurat not installed - returning metadata only.")
    return(invisible(meta))
  }
  if (is.null(counts))
    counts <- Seurat::ReadMtx(mtx = mtx, features = features, cells = cells,
                              feature.column = feature.column)
  # ReadMtx takes the barcodes straight from the file, so the GEM-well suffix
  # survives here even though it was stripped from `gex` above. Left alone, the
  # columns and the metadata rows are keyed differently and CreateSeuratObject
  # silently attaches a frame of NAs - the per-cell summary printed above would
  # still look correct. Strip it on both sides.
  if (!all(keep_cells)) {
    # restrict the matrix to the same GEM well the metadata was built from
    sfx_c <- ifelse(grepl("-[0-9]+$", colnames(counts)), sub("^.*-", "", colnames(counts)), NA)
    counts <- counts[, !is.na(sfx_c) & sfx_c == gem.well, drop = FALSE]
  }
  colnames(counts) <- sub("-[0-9]+$", "", colnames(counts))
  if (anyDuplicated(colnames(counts)))
    stop("Stripping the GEM-well suffix left ", sum(duplicated(colnames(counts))),
         " duplicated barcode(s). Pass `gem.well=` to pick one well of an ",
         "aggregated matrix.", call. = FALSE)
  if (!is.null(protein.coding))
    counts <- counts[rownames(counts) %in% protein.coding, , drop = FALSE]

  md <- as.data.frame(meta); rownames(md) <- md$cellID
  obj <- Seurat::CreateSeuratObject(counts = counts, project = project, assay = "RNA",
                                    meta.data = md[colnames(counts), , drop = FALSE])
  obj$sample <- sample_id

  if (subset.valid) {
    # The matrix and catcheR do not cover the same cells: the matrix holds every
    # cell the chemistry called, catcheR only those whose cassette was read and
    # survived the filter. Say what is being dropped rather than silently
    # returning a smaller object.
    keep <- colnames(obj)[obj$ValidCell | obj$is_empty]
    if (verbose) {
      drop <- setdiff(colnames(obj), keep)
      message(sprintf("  subset.valid: keeping %s of %s cells (%s dropped - %s)",
                      format(length(keep), big.mark = ","),
                      format(ncol(obj), big.mark = ","),
                      format(length(drop), big.mark = ","),
                      paste(sprintf("%s %s", table(obj$catcheR_status[match(drop, colnames(obj))]),
                                    names(table(obj$catcheR_status[match(drop, colnames(obj))]))),
                            collapse = ", ")))
    }
    obj <- subset(obj, cells = keep)
  }
  if (rename) {
    nm <- ifelse(obj$ValidCell,
                 paste(colnames(obj), obj$UMIxUCI, obj$barcode, obj$gene, obj$UCI, sep = "_"),
                 ifelse(obj$is_empty, paste(colnames(obj), "?", "empty", "empty", "empty", sep = "_"),
                        colnames(obj)))
    obj <- Seurat::RenameCells(obj, new.names = nm)
  }
  if (isTRUE(write.matrix.csv)) {
    # catcheR's silencing_matrix_<sample>.csv. Genes x cells, renamed columns.
    # Dense and large: ~nrow x ncol numbers. Prefer the Seurat .rds unless a
    # downstream tool genuinely needs a flat file.
    m <- as.matrix(Seurat::GetAssayData(obj, layer = "counts"))
    utils::write.csv(m, op(sprintf("silencing_matrix_%s.csv", sample_id)))
    if (verbose) message(sprintf("  silencing_matrix_%s.csv written (%d x %d - dense)",
                                 sample_id, nrow(m), ncol(m)))
  }
  saveRDS(obj, op(sprintf("%s_seurat_%s.rds", pfx, sample_id)))
  if (verbose)
    message(sprintf("\nSeurat object: %d genes x %d cells -> %s",
                    nrow(obj), ncol(obj), op(sprintf("%s_seurat_%s.rds", pfx, sample_id))))
  invisible(obj)
}

# ── 5. catcheR_pip_report ───────────────────────────────────────────────────

#' Render the step-by-step pipcatch report
#'
#' Thin wrapper over the parameterised Rmd. Knitting re-runs the extraction, so
#' this is the reproducible record rather than a summary of an object already
#' in memory.
#' @export
catcheR_pip_report <- function(folder, fastq.read1, fastq.read2,
                               ref_csv,
                               set_names  = c("set1", "set2"),
                               set_sizes  = NULL,
                               outdir     = "pipcatch_out",
                               output_file = "catcheR_pipcatch_report.html",
                               rmd        = NULL,
                               percentage = 15, ratio = 2, umi_count = 0,
                               max_reads  = Inf, ...) {
  if (!requireNamespace("rmarkdown", quietly = TRUE))
    stop("rmarkdown is required.", call. = FALSE)
  if (is.null(rmd)) {
    rmd <- system.file("rmd", "catcheR_pipcatch_report.Rmd", package = "catcheR")
    if (!nzchar(rmd)) rmd <- "catcheR_pipcatch_report.Rmd"
  }
  if (!file.exists(rmd)) stop("Report template not found: ", rmd, call. = FALSE)
  pars <- list(fastq_r1 = file.path(folder, fastq.read1),
               fastq_r2 = file.path(folder, fastq.read2),
               ref_csv  = file.path(folder, ref_csv),
               set_names = set_names, outdir = outdir,
               percentage = percentage, ratio = ratio, umi_count = umi_count,
               max_reads = max_reads)
  if (!is.null(set_sizes)) pars$set_sizes <- as.integer(set_sizes)
  rmarkdown::render(rmd, output_file = output_file, params = pars,
                    envir = new.env(), ...)
}
