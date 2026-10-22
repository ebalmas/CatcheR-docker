## Part of catcheR. See R/catcheR_pipcatch.R for the entry point and
## inst/rmd/catcheR_pipcatch_report.Rmd for the report template.

# ── 1. catcheR_pipcatch ─────────────────────────────────────────────────────

#' Run the PIPseq CATCHER extraction
#'
#' @param folder          directory holding the FASTQs and the barcode CSVs
#' @param fastq.read1,fastq.read2 file names within `folder`
#' @param barcodes.genes  one or more barcode->shRNA CSVs (no header, 2 columns).
#'                        Orientation is resolved automatically per file.
#' @param set.names,set.sizes labels and library sizes for the pooled sets
#' @param reference,empty.reference,constant construct landmarks
#' @param UCI.length,barcode.length  field widths
#' @param percentage,ratio,UMI.count  catcheR filter parameters applied for the
#'        QC pass. `catcheR_pip_filtercatch()` re-applies refined values.
#'        NOTE the defaults here deviate from `catcheR_10Xcatch()`: `ratio` is 2
#'        rather than catcheR's documented 5, and `UMI.count` is 0. At PIPseq
#'        enrichment depth `ratio = 5` discards every multi-UCI cell (no cell
#'        has the >=6 reads needed to reach 5:1). At `UMI.count = 0` the ratio
#'        test is inert either way - every surviving cell sits at
#'        UMIpercentagexUCI == 100 and takes the escape clause - but the two
#'        values diverge as soon as the threshold is raised (on AB019b,
#'        threshold 1 gives 499 ValidCells at ratio 2 vs 469 at ratio 5).
#'        Pass `ratio = 5` to reproduce catcheR's documented behaviour.
#' @param zoom.x,zoom.y  axis limits for the zoomed `UMIxUCI_<x>_<y>` panel.
#' @param jitter          geometry for the 2D QC panels. `"auto"` (default) uses
#'        catcheR's `geom_count()` when the percentage/UMIxUCI plane is densely
#'        populated and jitters the points when it is not. At PIPseq enrichment
#'        depth UMIxUCI takes only a few integer values and the percentage only
#'        simple fractions, so thousands of pairs collapse onto a few dozen
#'        coordinates and `geom_count()` degenerates into a handful of very large
#'        bubbles; jittering shows the same data as a distribution. Force the
#'        choice with `TRUE` / `FALSE`.
#' @param called.cells   optional: the cells PIPseeker called, as a path to
#'        its `barcodes.tsv.gz` (the same 28 nt concatenated 4-tier barcode
#'        this function emits), a plain file of IDs, or a character vector.
#'        Supplying it unlocks the panels judged against a called-cell
#'        denominator; see [catcheR_parsecatch()].
#' @param mode            "noise" or "bimodal"; only used to report what the
#'        automatic estimator WOULD pick. The returned object always carries
#'        both, and at low depth both are unusable - see the printed note.
#' @return invisibly, a `pipcatch` list; also written to
#'         `<output_dir>/pipcatch_<sample_name>.rds`
#' @export
catcheR_pipcatch <- function(folder,
                             fastq.read1,
                             fastq.read2,
                             barcodes.genes  = "rc_barcodes_genes.csv",
                             set.names       = NULL,
                             set.sizes       = NULL,
                             output_dir      = "Output/",
                             subdir          = TRUE,
                             sample_name     = "1",
                             reference       = "GGCGCGTTCATCTGGGGGAGCCG",
                             empty.reference = "TACGCGTTCATCTGGGGGAGCCG",
                             constant        = "ACGCGTAGCTCGCTGATCAG",
                             UCI.length      = 6,
                             barcode.length  = 8,
                             percentage      = 15,
                             ratio           = 2,
                             UMI.count       = 0,
                             mode            = c("noise", "bimodal"),
                             gene.sep        = "_",
                             called.cells    = NULL,
                             zoom.x          = NULL,
                             zoom.y          = NULL,
                             jitter          = "auto",
                             chunk.reads     = 1e6,
                             max.reads       = Inf,
                             plots           = TRUE,
                             verbose         = TRUE) {

  mode <- match.arg(mode); ptm <- proc.time()
  if (!dir.exists(folder)) stop("Folder not found: ", folder, call. = FALSE)
  folder <- normalizePath(folder, mustWork = TRUE)
  f1 <- file.path(folder, fastq.read1); f2 <- file.path(folder, fastq.read2)
  refp <- file.path(folder, barcodes.genes)
  for (f in c(f1, f2, refp)) if (!file.exists(f)) stop("Missing file: ", f, call. = FALSE)
  if (is.null(set.names)) set.names <- if (length(refp) == 1) "set1" else paste0("set", seq_along(refp))
  # One subfolder per chemistry. All three pipelines write the 13 plots that
  # carry catcheR's own file names (barcode_distribution, UMIxUCI, ...), which
  # have no chemistry in the filename by design, so running 10x and Parse into
  # the same Output/ used to have the second silently replace the first.
  # Everything downstream reads `x$output_dir`, so nesting once here propagates.
  if (isTRUE(subdir)) output_dir <- file.path(output_dir, "pipcatch")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  op <- function(f) file.path(output_dir, f)

  RCREF    <- .cat_rc(reference)
  RC_EMPTY <- .cat_rc(empty.reference)
  ANCH     <- substr(RCREF, 1, nchar(RCREF) - barcode.length)
  UL <- UCI.length; BL <- barcode.length

  if (verbose) message("Streaming FASTQ pair ...")
  st <- .cat_stream(f1, f2, .pip_adapter(), ANCH, RC_EMPTY, constant, UL, BL,
                   chunk.reads, max.reads, verbose)
  if (!nrow(st$cassette)) stop("No cassette-bearing reads found. Check `reference` and `constant`.", call. = FALSE)
  st$cassette$UCI     <- .cat_rc(st$cassette$UCI_rc)
  st$cassette$barcode <- .cat_rc(st$cassette$bc_rc)

  if (verbose) message("Resolving reference orientation ...")
  ref <- .cat_load_reference(refp, set.names, unique(st$cassette$barcode), gene.sep)
  utils::write.table(ref[, c("barcode", "shRNA")], op("rc_barcodes_genes_merged.csv"), sep = ",", row.names = FALSE,
                     col.names = FALSE, quote = FALSE)

  colfile <- file.path(folder, "colors.csv")
  cols <- if (file.exists(colfile))
    utils::read.csv(colfile, header = FALSE, col.names = c("shRNA", "colors"),
                    colClasses = "character") else NULL

  T_reads <- st$cassette %>% left_join(ref, by = "barcode")
  assign_tbl <- T_reads %>% filter(!is.na(shRNA)) %>%
    select(cellID, index, barcode, UCI, shRNA, set, gene, construct)   # one row per READ
  utils::write.csv(assign_tbl, op("pipcatch_cell_assignments.csv"), row.names = FALSE)

  ## read-class summary -------------------------------------------------------
  assigned_cells <- unique(assign_tbl$cellID)
  split_cass <- c("assigned" = sum(!is.na(T_reads$shRNA)),
                  "cassette present, barcode not in reference" = sum(is.na(T_reads$shRNA)))
  read_class <- bind_rows(
    data.frame(class = names(st$class_counts), reads = as.integer(unname(st$class_counts))) %>%
      filter(class != "cassette") %>% mutate(reads = as.integer(reads)),
    data.frame(class = names(split_cass), reads = as.integer(unname(split_cass))) %>% mutate(reads = as.integer(reads))
  ) %>% arrange(desc(reads)) %>% mutate(pct = round(100 * reads / sum(reads), 2))
  utils::write.csv(read_class, op("pipcatch_read_class_summary.csv"), row.names = FALSE)

  ## per-cell fates (deterministic tie-breaking) ------------------------------
  CLASS_PRIORITY <- c("assigned", "cassette present, barcode not in reference",
                      "empty vector (no insertion)", "primer + poly-A chimera",
                      "primer, no cassette, no poly-A", "non-construct insert (poly-A)",
                      "non-construct insert")
  cass_split <- T_reads %>%
    transmute(cell = cellID,
              cls = ifelse(is.na(shRNA), "cassette present, barcode not in reference", "assigned")) %>%
    count(cell, cls, name = "n")
  cell_class <- bind_rows(filter(st$cell_class, cls != "cassette"), cass_split)
  nocall <- cell_class %>% filter(!cell %in% assigned_cells) %>%
    mutate(prio = match(cls, CLASS_PRIORITY)) %>%
    group_by(cell) %>% arrange(desc(n), prio, .by_group = TRUE) %>% slice(1) %>% ungroup() %>%
    count(cls, name = "cells") %>% arrange(desc(cells)) %>%
    mutate(pct = round(100 * cells / sum(cells), 1))

  ## catcheR tables -----------------------------------------------------------
  tb  <- .cat_tabulate(assign_tbl); cu <- tb$cu; tab <- tb$tab
  thr <- .cat_thresholds(tab)
  val <- .cat_validate(tab, UMI.count, percentage, ratio)
  complete_table_fin <- val$pairs %>%
    left_join(cu %>% select(cellID, UCI, barcode, shRNA, set, gene), by = c("cellID", "UCI"))
  utils::write.csv(complete_table_fin, op(sprintf("complete_table_fin_%s.csv", sample_name)), row.names = FALSE)

  per_cell <- cu %>% group_by(cellID) %>% arrange(desc(reads), .by_group = TRUE) %>%
    summarise(top_barcode = first(barcode), top_UCI = first(UCI), top_shRNA = first(shRNA),
              top_set = first(set), gene = first(gene),
              reads_top = first(reads), reads_second = if (n() > 1) sort(reads, decreasing = TRUE)[2] else 0L,
              reads_total = sum(reads), n_UCI = n(), n_sets = n_distinct(set), .groups = "drop") %>%
    mutate(pct_top = 100 * reads_top / reads_total,
           ratio_top = ifelse(reads_second > 0, reads_top / reads_second, Inf)) %>%
    left_join(val$cells %>% select(cellID, True_UCIxcellID, ValidCell), by = "cellID")
  utils::write.csv(per_cell, op("pipcatch_per_cell_table.csv"), row.names = FALSE)

  uci_tab <- per_cell %>% count(n_UCI, name = "cells") %>%
    mutate(pct = round(100 * cells / sum(cells), 1))
  clone_summary <- assign_tbl %>% group_by(set, gene, shRNA, construct) %>%
    summarise(clones = n_distinct(UCI), cells = n_distinct(cellID), .groups = "drop") %>%
    arrange(desc(clones), desc(cells))
  utils::write.csv(clone_summary, op("pipcatch_shRNA_clone_summary.csv"), row.names = FALSE)

  res <- list(sample = sample_name, output_dir = output_dir,
              params = list(reference = reference, empty.reference = empty.reference,
                            constant = constant, UCI.length = UL, barcode.length = BL,
                            percentage = percentage, ratio = ratio, UMI.count = UMI.count,
                            mode = mode, gene.sep = gene.sep,
                            set.names = set.names, set.sizes = set.sizes,
                            colors = cols),
              total_pairs = st$total_pairs,
              valid_cell_reads = st$valid_cell_reads,
              ## PIPseq is a single library, but the shared analysis panels
              ## expect the same per-library read funnel the multi-sublibrary
              ## chemistries produce, so give them a one-row version.
              per_sublibrary = data.frame(
                sublibrary     = "all",
                pairs          = st$total_pairs,
                cell_bc_reads  = st$valid_cell_reads,
                cassette_reads = nrow(st$cassette)),
              reference_table = ref,
              cell_class = cell_class, read_class = read_class, nocall = nocall,
              assignments = assign_tbl, empty_reads = st$empty,
              cu = cu, tab = tab, thresholds = thr,
              complete_table_fin = complete_table_fin, per_cell = per_cell,
              uci_table = uci_tab, clone_summary = clone_summary,
              n_cell_barcodes = n_distinct(cell_class$cell),
              n_cells_called = length(assigned_cells),
              n_valid = sum(per_cell$ValidCell, na.rm = TRUE))
  class(res) <- c("pipcatch", "list")

  res$input_files <- .cat_input_manifest(
    list("fastq read 1" = f1, "fastq read 2" = f2, "barcode-to-gene reference" = refp,
         "called-cell list" = if (is.character(called.cells)) called.cells else NULL),
    sample_name, "pipcatch", res$output_dir)

  if (plots) .cat_plots_qc(res, zoom.x = zoom.x, zoom.y = zoom.y, jitter = jitter)
  # Keep the resolved called-cell list on the object: the re-filter draws the
  # threshold-choice panel and without this it would fall back to every observed
  # barcode, ambient floor included, which overstates coverage.
  res$called_cells <- if (!is.null(called.cells)) .par_called(called.cells) else NULL
  if (plots) .par_plots_analysis(res, called = called.cells, verbose = verbose)
  saveRDS(res, op(sprintf("pipcatch_%s.rds", sample_name)))

  if (verbose) {
    message(sprintf("\nValidated cell barcodes in R1 : %s", format(res$n_cell_barcodes, big.mark = ",")))
    message(sprintf("  with an shRNA call          : %s (%.1f%%)",
                    format(res$n_cells_called, big.mark = ","),
                    100 * res$n_cells_called / res$n_cell_barcodes))
    message(sprintf("  ValidCell at UMI.count=%s     : %s", UMI.count, format(res$n_valid, big.mark = ",")))
    message(sprintf("Constructs detected           : %d of %d", nrow(clone_summary), nrow(ref)))
    message(sprintf("Clones (barcode-UCI)          : %d", sum(clone_summary$clones)))
    message(sprintf("\nAutomatic thresholds: noise = %s, bimodal = %s  (%.1f%% of UCIs are singletons)",
                    ifelse(is.na(thr$noise), "UNDEFINED (>74.07% singletons)", thr$noise),
                    ifelse(is.na(thr$bimodal), "no valley", thr$bimodal), thr$pct_singletons))
    if (is.na(thr$noise) || isTRUE(thr$bimodal > 1))
      message("  -> depth is too low for the automatic modes; keep UMI.count = 0 and sequence deeper.")
    message(sprintf("\nDone in %.1f min.", (proc.time() - ptm)[["elapsed"]] / 60))
  }
  invisible(res)
}

#' @export
print.pipcatch <- function(x, ...) {
  cat(sprintf("<pipcatch> sample %s | %s read pairs | %s validated cell barcodes\n",
              x$sample, format(x$total_pairs, big.mark = ","),
              format(x$n_cell_barcodes, big.mark = ",")))
  cat(sprintf("  %s cells with a call, %s ValidCell, %d constructs, %d clones\n",
              format(x$n_cells_called, big.mark = ","), format(x$n_valid, big.mark = ","),
              nrow(x$clone_summary), sum(x$clone_summary$clones)))
  invisible(x)
}

