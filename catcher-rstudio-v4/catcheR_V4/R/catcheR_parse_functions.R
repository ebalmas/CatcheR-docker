## Part of catcheR. Downstream of [catcheR_parsecatch()].

# ── 2. catcheR_parse_filtercatch ────────────────────────────────────────────

#' Re-filter a parsecatch object at new thresholds
#'
#' Mirrors [catcheR_pip_filtercatch()] and adds the per-cell UMI floor, which
#' PIPseq has no use for (its molecular index is 3 nt and not deduplicable) but
#' which is the single most consequential parameter on Parse data.
#'
#' @section The three thresholds, exactly:
#' Two of them exclude with a STRICT `>` and one with an inclusive `>=`, which
#' is catcheR's own behaviour and is easy to get wrong:
#' \itemize{
#'   \item `UMI.count` - a (cell, UCI) pair is kept when `UMIxUCI > UMI.count`.
#'     It is therefore the highest count that is still **excluded**, not a
#'     minimum: `UMI.count = 1` keeps pairs of **2 or more** UMIs, and
#'     `UMI.count = 0` keeps everything down to a single UMI.
#'   \item `percentage` - kept when `UMIpercentagexUCI > percentage`, also
#'     strict. A UCI holding exactly 15% of its cell does **not** pass
#'     `percentage = 15`.
#'   \item `min.cell.UMI` - a cell enters the comparison when
#'     `UMIxcellID >= min.cell.UMI`, **inclusive**. `min.cell.UMI = 20` keeps
#'     cells of exactly 20 cassette UMIs.
#' }
#' `ratio` applies only once a cell is down to one surviving UCI: that UCI must
#' hold at least `ratio` times the UMIs of the runner-up, unless it holds 100%
#' of the cell, in which case there is no runner-up to compare against.
#'
#' Note that the floor IS carried forward by default: both the extraction entry
#' points and `catcheR_*_filtercatch()` default `min.cell.UMI` to 20, so a
#' re-filter applies the same floor unless you pass a different value (`0`
#' disables it). The two steps do NOT share state - the default simply matches.
#' @param x a `parsecatch` object or the path to its `.rds`.
#' @param UMI.count highest UMI count on a (cell, UCI) pair that is still
#'   EXCLUDED: the pair is kept when `UMIxUCI > UMI.count`, strict as in
#'   catcheR. `UMI.count = 1` therefore keeps pairs of 2 or more UMIs.
#' @param percentage share of the cell's UMIs a UCI must EXCEED, also strict:
#'   a UCI at exactly `percentage` does not pass.
#' @param ratio dominant-to-second ratio required once a cell is down to one
#'   surviving UCI. Defaults to 2, NOT to `catcheR_10Xcatch()`'s documented 5:
#'   at enrichment depth a 5:1 ratio needs at least 6 UMIs on the dominant UCI,
#'   which discards essentially every multi-UCI cell. Pass `ratio = 5` to
#'   reproduce the 10X behaviour exactly.
#' @param min.cell.UMI minimum cassette UMIs in the cell, INCLUSIVE
#'   (`UMIxcellID >= min.cell.UMI`). Defaults to 20, the same as the extraction
#'   entry points, so the floor keeps applying unless you pass another value;
#'   `0` disables it.
#' @param sweep extra `UMI.count` values to report alongside the chosen one.
#' @param to.scratch optional path to the experiment's `scratch/`. When given,
#'   the files a later step reads - `cells_with_at_least_1_UCI.csv`,
#'   `new_cell_names.csv`, `complete_table_fin_<sample>.csv` - are COPIED there
#'   with `_<chemistry>_<sample>` appended, so two pipelines can share one
#'   scratch. The originals stay in `output_dir`.
#' @param subdir if TRUE (default) results go to a subfolder of `output_dir`
#'   named for the thresholds, e.g. `filtercatch_p15_u1_r2_f20/`, so repeated
#'   re-filtering does not overwrite earlier attempts. The path is returned
#'   as `$output_dir`.
#' @param output_dir,jitter,plots,verbose as elsewhere.
#' @return a `parse_filtercatch` object, also saved as `.rds`.
#' @export
catcheR_parse_filtercatch <- function(x,
                                      UMI.count    = 1,
                                      percentage   = 15,
                                      ratio        = 2,
                                      min.cell.UMI = 20,
                                      output_dir   = NULL,
                                      subdir       = TRUE,
                                      to.scratch   = NULL,
                                      sweep        = c(0, 1),
                                      jitter       = "auto",
                                      plots        = TRUE,
                                      verbose      = TRUE) {
  if (is.character(x)) x <- readRDS(x)
  # Was the floor chosen, or is it just the default? 20 is a reasonable starting
  # point, not a determination, so when it is not supplied the function says so
  # and writes the figure that lets the value be chosen on this dataset rather
  # than taken on trust.
  floor_defaulted <- missing(min.cell.UMI)
  # every chemistry's object carries "pipcatch" as its base class
  stopifnot(inherits(x, "pipcatch"))
  if (is.null(output_dir)) output_dir <- x$output_dir
  pfx <- if (inherits(x, "tenxcatch")) "tenxcatch"
         else if (inherits(x, "parsecatch")) "parsecatch" else "pipcatch"
  # Re-filtering is a working step: each threshold set gets its own subfolder,
  # tagged with the parameters, so trying 15/1/2 and then 12/0/2 leaves two
  # readable result sets instead of one overwritten by the other. The `.rds`
  # and the plots inside are the ones for THOSE thresholds.
  if (isTRUE(subdir)) {
    tag <- sprintf("filtercatch_p%s_u%s_r%s", percentage, UMI.count, ratio)
    if (!is.null(min.cell.UMI) && isTRUE(min.cell.UMI > 0))
      tag <- paste0(tag, "_f", min.cell.UMI)
    output_dir <- file.path(output_dir, tag)
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  op <- function(f) file.path(output_dir, f)

  depth <- x$tab %>% group_by(cellID) %>%
    summarise(UMIxcellID = first(UMIxcellID), .groups = "drop")
  apply_floor <- function(v) {
    v$cells <- v$cells %>% left_join(depth, by = "cellID") %>%
      mutate(ValidCell = ValidCell & UMIxcellID >= min.cell.UMI)
    v$pairs <- v$pairs %>% select(-any_of("ValidCell")) %>%
      left_join(v$cells %>% select(cellID, ValidCell), by = "cellID")
    v
  }
  val <- apply_floor(.cat_validate(x$tab, UMI.count, percentage, ratio))
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
                    chem = pfx, sample = x$sample, verbose = verbose)
  # The input manifest is written by the EXTRACTION step, so it sits one level
  # up from this filter's own folder - copy it from there, not from here.
  if (!is.null(to.scratch))
    .cat_to_scratch(stats::setNames("input_files",
                                    sprintf("%s_%s_input_files.csv", x$sample, pfx)),
                    from = x$output_dir, to = to.scratch,
                    chem = pfx, sample = x$sample, verbose = verbose)

  thr_grid <- sort(unique(c(UMI.count, sweep, x$thresholds$bimodal)))
  thr_grid <- thr_grid[!is.na(thr_grid)]
  filter_returns <- bind_rows(lapply(thr_grid, function(t) {
    r <- apply_floor(.cat_validate(x$tab, t, percentage, ratio))
    data.frame(threshold = t, keeps_pairs_with_min_UMIs = t + 1L,
               TorF_pairs = sum(r$pairs$TorF),
               cells_with_true_UCI = sum(r$cells$True_UCIxcellID > 0),
               ValidCell = sum(r$cells$ValidCell))
  }))
  utils::write.csv(filter_returns, op(paste0(pfx, "_filter_returns.csv")), row.names = FALSE)
  if (plots) {
    pp <- modifyList(x$params, list(UMI.count = UMI.count, percentage = percentage, ratio = ratio))
    .cat_plots_stage2(ctf, val$cells, percentage, UMI.count, output_dir, jitter)
    .cat_plots_thresholds(x$tab, x$thresholds, ctf, pp, output_dir,
                          file = paste0(pfx, "_filtercatch_thresholds.png"))
    if (floor_defaulted) {
      # catcheR_call_accuracy() subsamples the deep cells and scores every
      # outcome correct / wrong / no call, so the panel shows what the floor
      # buys in accuracy against what it costs in yield. It legitimately cannot
      # run on a library with no cell deep enough to serve as ground truth - say
      # which panel is missing and why rather than omitting it silently.
      acc <- try(catcheR_call_accuracy(x, UMI.count = UMI.count,
                                       percentage = percentage, ratio = ratio,
                                       output_dir = output_dir, verbose = FALSE),
                 silent = TRUE)
      if (inherits(acc, "try-error")) {
        if (verbose)
          message("  min.cell.UMI not supplied, using ", min.cell.UMI,
                  "; no threshold-choice panel: ",
                  conditionMessage(attr(acc, "condition")))
      } else {
        .par_plot_yield_error(x, acc, output_dir, min.cell.UMI,
                              called = x$called_cells, prefix = pfx)
        if (verbose)
          message("  min.cell.UMI not supplied, using ", min.cell.UMI,
                  " - see ", paste0(pfx, "_yield_vs_error.png"), " and ",
                  paste0("catcheR_call_accuracy_", x$sample, ".csv"),
                  " in ", output_dir, " to confirm or replace it on this dataset")
      }
    }
  }
  out <- list(sample = x$sample, output_dir = output_dir,
              params = list(UMI.count = UMI.count, percentage = percentage,
                            ratio = ratio, min.cell.UMI = min.cell.UMI),
              complete_table_fin = ctf, cells = val$cells,
              valid_cells = valid, new_names = new_names, filter_returns = filter_returns)
  class(out) <- c("parse_filtercatch", "pip_filtercatch", "list")
  saveRDS(out, op(sprintf("%s_filtercatch_%s.rds", pfx, x$sample)))
  if (verbose) {
    message(sprintf("UMI.count=%s percentage=%s ratio=%s cell floor=%s",
                    UMI.count, percentage, ratio, min.cell.UMI))
    message(sprintf("  cells with >=1 true UCI : %s",
                    format(sum(val$cells$True_UCIxcellID > 0), big.mark = ",")))
    message(sprintf("  ValidCell               : %s", format(length(valid), big.mark = ",")))
    print(as.data.frame(filter_returns))
  }
  invisible(out)
}

# ── 3. catcheR_parse_nocatch ────────────────────────────────────────────────

#' Unperturbed (empty-vector) cells in a Parse library
#'
#' Identical in logic to [catcheR_pip_nocatch()]; kept as its own entry point
#' so the Parse workflow reads consistently and the output files are named for
#' the chemistry.
#'
#' @param x a `parsecatch` object or path to its `.rds`.
#' @param filtered optional `parse_filtercatch` object.
#' @param threshold minimum empty-reference UMIs in the cell.
#' @param output_dir,verbose as elsewhere.
#' @export
catcheR_parse_nocatch <- function(x, filtered = NULL, threshold = 0,
                                  output_dir = NULL, verbose = TRUE) {
  if (is.character(x)) x <- readRDS(x)
  # every chemistry's object carries "pipcatch" as its base class
  stopifnot(inherits(x, "pipcatch"))
  catcheR_pip_nocatch(x, filtered = filtered, threshold = threshold,
                      output_dir = output_dir, verbose = verbose)
}

# ── 4. catcheR_parse_seurat ─────────────────────────────────────────────────

#' Join parsecatch calls to the split-pipe expression matrix
#'
#' split-pipe writes `count_matrix.mtx` as CELLS x GENES, the transpose of the
#' 10X and PIPseeker convention, so `Seurat::ReadMtx()` cannot read it and the
#' matrix is transposed on load. Cells are keyed by `bc_wells`, which is
#' exactly what [catcheR_parsecatch()] emits, so the join needs no translation.
#'
#' @param x a `parsecatch` object or path to its `.rds`.
#' @param filtered,empty optional objects from the two functions above.
#' @param mtx,genes,cells the three split-pipe DGE files.
#' @param project,sample_id labels for the Seurat object.
#' @param protein.coding optional gene symbols to keep.
#' @param rename rename cells to catcheR's
#'   `cellID_UMIxUCI_barcode_gene_UCI` form.
#' @param subset.valid keep only ValidCell (and empty, if supplied) cells.
#' @param write.matrix.csv also write the dense genes x cells CSV.
#' @param output_dir,verbose as elsewhere.
#' @return a Seurat object if Seurat is installed, otherwise the metadata
#'   table. The metadata CSV is always written.
#' @export
catcheR_parse_seurat <- function(x, filtered = NULL, empty = NULL,
                                 mtx, genes, cells,
                                 project = "parsecatch", sample_id = NULL,
                                 protein.coding = NULL, rename = FALSE,
                                 subset.valid = FALSE, write.matrix.csv = FALSE,
                                 output_dir = NULL, verbose = TRUE) {
  if (is.character(x)) x <- readRDS(x)
  # every chemistry's object carries "pipcatch" as its base class
  stopifnot(inherits(x, "pipcatch"))
  if (is.character(filtered)) filtered <- readRDS(filtered)
  if (is.character(empty))    empty    <- readRDS(empty)
  if (is.null(output_dir)) output_dir <- x$output_dir
  pfx <- if (inherits(x, "tenxcatch")) "tenxcatch"
         else if (inherits(x, "parsecatch")) "parsecatch" else "pipcatch"
  if (is.null(sample_id))  sample_id  <- x$sample
  op <- function(f) file.path(output_dir, f)

  dge <- .par_read_dge(mtx, genes, cells)
  gex <- colnames(dge$counts)

  pc <- x$per_cell
  if (!is.null(filtered))
    pc <- pc %>% select(-any_of(c("True_UCIxcellID", "ValidCell"))) %>%
      left_join(filtered$cells %>% select(cellID, True_UCIxcellID, ValidCell), by = "cellID")
  meta <- data.frame(cellID = gex, stringsAsFactors = FALSE) %>%
    left_join(pc %>% select(cellID, sublibrary, shRNA = top_shRNA, gene, set = top_set,
                            UCI = top_UCI, barcode = top_barcode,
                            UMIxUCI = reads_top, UMI_total = reads_total,
                            n_UCI, n_sets, True_UCIxcellID, ValidCell), by = "cellID") %>%
    mutate(ValidCell = tidyr::replace_na(ValidCell, FALSE),
           UMI_total = tidyr::replace_na(UMI_total, 0))
  meta$is_empty <- if (!is.null(empty) && nrow(empty$empty_cells))
    meta$cellID %in% empty$empty_cells$cellID else FALSE
  meta <- meta %>% mutate(catcheR_status = case_when(
    ValidCell                 ~ "valid",
    !is.na(shRNA) & n_UCI > 1 ~ "ambiguous",
    !is.na(shRNA)             ~ "called_not_valid",
    is_empty                  ~ "empty",
    TRUE                      ~ "not_detected"))
  meta <- meta %>% left_join(dge$meta %>% select(cellID = bc_wells, sample, gene_count, tscp_count),
                             by = "cellID")
  utils::write.csv(meta, op(sprintf("%s_cell_metadata_%s.csv", pfx, sample_id)), row.names = FALSE)

  if (verbose) {
    message(sprintf("split-pipe called cells: %s", format(length(gex), big.mark = ",")))
    message(sprintf("  with a parsecatch call : %s (%.1f%%)",
                    format(sum(!is.na(meta$shRNA)), big.mark = ","), 100 * mean(!is.na(meta$shRNA))))
    message(sprintf("  ValidCell              : %s (%.1f%%)",
                    format(sum(meta$ValidCell), big.mark = ","), 100 * mean(meta$ValidCell)))
    lost <- setdiff(x$per_cell$cellID, gex)
    message(sprintf("  parsecatch barcodes NOT called as cells by split-pipe: %s (ambient)",
                    format(length(lost), big.mark = ",")))
    print(meta %>% count(catcheR_status, name = "cells") %>% arrange(desc(cells)) %>% as.data.frame())
  }

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    message("Seurat not installed - returning metadata only.")
    return(invisible(meta))
  }
  counts <- dge$counts
  if (!is.null(protein.coding))
    counts <- counts[rownames(counts) %in% protein.coding, , drop = FALSE]
  md <- as.data.frame(meta); rownames(md) <- md$cellID
  obj <- Seurat::CreateSeuratObject(counts = counts, project = project, assay = "RNA",
                                    meta.data = md[colnames(counts), , drop = FALSE])
  obj$sample_id <- sample_id
  if (subset.valid) obj <- subset(obj, cells = colnames(obj)[obj$ValidCell | obj$is_empty])
  if (rename) {
    nm <- ifelse(obj$ValidCell,
                 paste(colnames(obj), obj$UMIxUCI, obj$barcode, obj$gene, obj$UCI, sep = "_"),
                 ifelse(obj$is_empty,
                        paste(colnames(obj), "?", "empty", "empty", "empty", sep = "_"),
                        colnames(obj)))
    obj <- Seurat::RenameCells(obj, new.names = nm)
  }
  if (isTRUE(write.matrix.csv)) {
    m <- as.matrix(Seurat::GetAssayData(obj, layer = "counts"))
    utils::write.csv(m, op(sprintf("silencing_matrix_%s.csv", sample_id)))
  }
  saveRDS(obj, op(sprintf("%s_seurat_%s.rds", pfx, sample_id)))
  if (verbose)
    message(sprintf("\nSeurat object: %d genes x %d cells -> %s",
                    nrow(obj), ncol(obj), op(sprintf("%s_seurat_%s.rds", pfx, sample_id))))
  invisible(obj)
}

# ── 5. catcheR_parse_report ─────────────────────────────────────────────────

#' Render the Parse report template
#'
#' Thin wrapper over the parameterised Rmd shipped in `inst/rmd`. Knitting
#' re-runs the extraction, so the document is a reproducible record rather
#' than a summary of an object already in memory.
#'
#' @param folder the analysis folder; `scratch/` and `Output/` are built from it.
#' @param output_file,rmd output and template paths.
#' @param ... passed to the template's params.
#' @export
catcheR_parse_report <- function(folder = "analysis_EB",
                                 output_file = "catcheR_parsecatch_report.html",
                                 rmd = NULL, ...) {
  if (!requireNamespace("rmarkdown", quietly = TRUE))
    stop("rmarkdown is required.", call. = FALSE)
  if (is.null(rmd)) {
    rmd <- system.file("rmd", "catcheR_parsecatch_analysis_EB.Rmd", package = "catcheR")
    if (!nzchar(rmd)) rmd <- "catcheR_parsecatch_analysis_EB.Rmd"
  }
  if (!file.exists(rmd)) stop("Report template not found: ", rmd, call. = FALSE)
  rmarkdown::render(rmd, output_file = output_file,
                    params = utils::modifyList(list(analysis_dir = folder), list(...)),
                    envir = new.env())
}
