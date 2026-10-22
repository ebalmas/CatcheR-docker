## Part of catcheR. 10x Genomics entry point, native R.

# ── catcheR_tenxcatch ───────────────────────────────────────────────────────

#' Extract CATCHER clone identities from a 10x barcode-enrichment library
#'
#' The native-R equivalent of [catcheR_10Xcatch()], which drives the bundled
#' shell scripts through Docker. Same construct logic and the same catcheR
#' thresholds; the difference is that this one streams the FASTQ pair in R, so
#' it needs no container, and it shares the extraction layer with
#' [catcheR_pipcatch()] and [catcheR_parsecatch()] - which means the three
#' chemistries produce directly comparable tables and figures.
#'
#' Read 1 is sliced at fixed offsets (`[cell barcode][UMI]`), exactly as the
#' shell pipeline's `cut -c 1-16` / `cut -c 17-28`. Read 2 carries the cassette
#' with the reference reverse-complemented at the read start.
#'
#' @param folder directory holding the FASTQ files and the reference CSVs.
#' @param fastq.read1,fastq.read2 file names. Several lanes of ONE sample may
#'   be given; they are streamed in turn and pooled, since the cell barcode is
#'   shared across lanes. Do not pass several samples - their barcodes collide.
#' @param barcodes.genes one or more barcode-to-gene CSVs. Orientation is
#'   resolved against the data, so forward or reverse-complement both work.
#' @param set.names,set.sizes labels and nominal sizes for those CSVs.
#' @param cell.whitelist optional 10x barcode whitelist, as a character vector
#'   or a path to `filtered_feature_bc_matrix.h5` / `barcodes.tsv.gz` / a plain
#'   text file. Strongly recommended: without it every structurally valid
#'   16-mer is accepted, so the observed-barcode count includes the ambient
#'   floor plus one spurious barcode per sequencing error, and no correction is
#'   possible. Cell Ranger's own `-1` suffix is stripped.
#' @param called.cells optional: the cells Cell Ranger called, same formats.
#'   Defaults to `cell.whitelist` when that is a Cell Ranger file, since the
#'   filtered matrix IS the called-cell list. Unlocks the panels judged against
#'   a called-cell denominator.
#' @param cb.length,umi.length read-1 field widths. 16/12 is v3 (and v2 used
#'   16/10); check your chemistry rather than trusting the default.
#' @param reference,empty.reference,constant construct sequences.
#' @param UCI.length,barcode.length cassette field widths.
#' @param max.edit.distance correction radius against `cell.whitelist`.
#' @param percentage,ratio,UMI.count catcheR's calling thresholds. Note
#'   `ratio` defaults to 2 here, not to `catcheR_10Xcatch()`'s documented 5;
#'   see [catcheR_pip_filtercatch()] for why.
#' @param min.cell.UMI minimum cassette UMIs for a cell to be callable,
#'   INCLUSIVE (`UMIxcellID >= min.cell.UMI`). See the threshold section below;
#'   `catcheR_*_filtercatch()` defaults to the same 20, so it keeps applying.
#' @param mode,gene.sep,zoom.x,zoom.y,jitter,chunk.reads,max.reads,plots,verbose
#'   as in [catcheR_parsecatch()].
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
#' @param output_dir,sample_name where to write and how to label.
#' @param subdir if TRUE (default) everything is written to
#'   `<output_dir>/tenxcatch/`, so the three chemistries cannot
#'   overwrite each other's plots. Set FALSE for the flat layout.
#' @return a `tenxcatch` object, also saved as `.rds`.
#' @export
catcheR_tenxcatch <- function(folder,
                              fastq.read1,
                              fastq.read2,
                              barcodes.genes  = "rc_barcodes_genes.csv",
                              set.names       = NULL,
                              set.sizes       = NULL,
                              cell.whitelist  = NULL,
                              called.cells    = NULL,
                              output_dir      = "Output/",
                              subdir          = TRUE,
                              sample_name     = "1",
                              reference       = "GGCGCGTTCATCTGGGGGAGCCG",
                              empty.reference = "TACGCGTTCATCTGGGGGAGCCG",
                              constant        = "ACGCGTAGCTCGCTGATCAG",
                              UCI.length      = 6,
                              barcode.length  = 8,
                              cb.length       = 16,
                              umi.length      = 12,
                              max.edit.distance = 1,
                              percentage      = 15,
                              ratio           = 2,
                              UMI.count       = 1,
                              min.cell.UMI    = 20,
                              mode            = c("noise", "bimodal"),
                              gene.sep        = "_",
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
  if (length(fastq.read1) != length(fastq.read2))
    stop("fastq.read1 and fastq.read2 must have the same length (one per lane).",
         call. = FALSE)
  f1 <- file.path(folder, fastq.read1); f2 <- file.path(folder, fastq.read2)
  refp <- file.path(folder, barcodes.genes)
  for (f in c(f1, f2, refp))
    if (!file.exists(f)) stop("Missing file: ", f, call. = FALSE)
  if (is.null(set.names))
    set.names <- if (length(refp) == 1) "set1" else paste0("set", seq_along(refp))
  # One subfolder per chemistry. All three pipelines write the 13 plots that
  # carry catcheR's own file names (barcode_distribution, UMIxUCI, ...), which
  # have no chemistry in the filename by design, so running 10x and Parse into
  # the same Output/ used to have the second silently replace the first.
  # Everything downstream reads `x$output_dir`, so nesting once here propagates.
  if (isTRUE(subdir)) output_dir <- file.path(output_dir, "tenxcatch")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  op <- function(f) file.path(output_dir, f)

  ## cell-barcode whitelist --------------------------------------------------
  wl <- if (is.character(cell.whitelist) && length(cell.whitelist) == 1L &&
            file.exists(cell.whitelist)) .tenx_called(cell.whitelist)
        else if (!is.null(cell.whitelist)) unique(sub("-[0-9]+$", "", cell.whitelist))
        else NULL
  cmap <- NULL
  if (!is.null(wl)) {
    base <- stats::setNames(wl, wl)
    cmap <- .par_corr_map(base, max.edit.distance = as.integer(max.edit.distance))
    if (verbose)
      message(sprintf("Cell barcodes: %s whitelist entries -> %s accepted with correction",
                      format(length(wl), big.mark = ","), format(length(cmap), big.mark = ",")))
  } else {
    warning("No cell.whitelist supplied. Every structurally valid ", cb.length,
            "-mer will be accepted, so the observed-barcode count will include ",
            "the ambient floor and one spurious barcode per sequencing error, ",
            "and no correction is applied.", call. = FALSE)
  }
  if (is.null(called.cells) && is.character(cell.whitelist) &&
      length(cell.whitelist) == 1L && file.exists(cell.whitelist))
    called.cells <- cell.whitelist

  RCREF    <- .cat_rc(reference)
  RC_EMPTY <- .cat_rc(empty.reference)
  ANCH     <- substr(RCREF, 1, nchar(RCREF) - barcode.length)
  UL <- UCI.length; BL <- barcode.length

  ## stream each lane --------------------------------------------------------
  ad <- .tenx_adapter(as.integer(cb.length), as.integer(umi.length), cmap)
  parts <- lapply(seq_along(f1), function(i) {
    if (verbose) message(sprintf("Streaming %s ...", basename(f1[i])))
    .cat_stream(f1[i], f2[i], ad, ANCH, RC_EMPTY, constant, UL, BL,
                chunk.reads, max.reads, verbose)
  })
  st <- list(
    total_pairs  = sum(vapply(parts, function(p) p$total_pairs, numeric(1))),
    class_counts = tapply(unlist(lapply(parts, function(p) unname(p$class_counts))),
                          unlist(lapply(parts, function(p) names(p$class_counts))), sum),
    cell_class   = bind_rows(lapply(parts, `[[`, "cell_class")) %>%
                   count(cell, cls, wt = n, name = "n"),
    cassette     = bind_rows(lapply(parts, `[[`, "cassette")),
    empty        = bind_rows(lapply(parts, `[[`, "empty")),
    valid_cell_reads = sum(vapply(parts, function(p) p$valid_cell_reads, numeric(1))))
  ## lanes of one sample share their cells, so the read funnel pools them into
  ## a single row rather than pretending they are independent libraries
  st$per_sublib <- data.frame(
    sublibrary     = "all",
    pairs          = st$total_pairs,
    cell_bc_reads  = st$valid_cell_reads,
    cassette_reads = nrow(st$cassette))
  if (!nrow(st$cassette))
    stop("No cassette-bearing reads found. Check `reference`, `constant` and the read ",
         "order (10x carries the cassette on read 2).", call. = FALSE)
  st$cassette$UCI     <- .cat_rc(st$cassette$UCI_rc)
  st$cassette$barcode <- .cat_rc(st$cassette$bc_rc)

  if (verbose) message("Resolving reference orientation ...")
  ref <- .cat_load_reference(refp, set.names, unique(st$cassette$barcode), gene.sep)
  utils::write.table(ref[, c("barcode", "shRNA")], op("rc_barcodes_genes_merged.csv"),
                     sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  colfile <- file.path(folder, "colors.csv")
  cols <- if (file.exists(colfile))
    utils::read.csv(colfile, header = FALSE, col.names = c("shRNA", "colors"),
                    colClasses = "character") else NULL

  T_reads <- st$cassette %>% left_join(ref, by = "barcode")

  ## UMI deduplication: `index` is a real UMI, so one row per MOLECULE --------
  assign_tbl <- T_reads %>% filter(!is.na(shRNA)) %>%
    select(cellID, index, barcode, UCI, shRNA, set, gene, construct) %>%
    distinct()
  n_reads_assigned <- sum(!is.na(T_reads$shRNA))
  mol_reads <- T_reads %>% filter(!is.na(shRNA)) %>%
    count(cellID, index, barcode, UCI, name = "reads")
  assign_tbl <- assign_tbl %>%
    left_join(mol_reads, by = c("cellID", "index", "barcode", "UCI"))
  utils::write.csv(assign_tbl, op("tenxcatch_cell_assignments.csv"), row.names = FALSE)

  ## read-class summary ------------------------------------------------------
  assigned_cells <- unique(assign_tbl$cellID)
  split_cass <- c("assigned" = n_reads_assigned,
                  "cassette present, barcode not in reference" = sum(is.na(T_reads$shRNA)))
  read_class <- bind_rows(
    data.frame(class = names(st$class_counts), reads = as.integer(unname(st$class_counts))) %>%
      filter(class != "cassette"),
    data.frame(class = names(split_cass), reads = as.integer(unname(split_cass)))
  ) %>% arrange(desc(reads)) %>% mutate(pct = round(100 * reads / sum(reads), 2))
  utils::write.csv(read_class, op("tenxcatch_read_class_summary.csv"), row.names = FALSE)

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

  ## catcheR tables ----------------------------------------------------------
  tb  <- .cat_tabulate(assign_tbl); cu <- tb$cu; tab <- tb$tab
  thr <- .cat_thresholds(tab)
  val <- .cat_validate(tab, UMI.count, percentage, ratio)
  depth <- tab %>% group_by(cellID) %>% summarise(UMIxcellID = first(UMIxcellID), .groups = "drop")
  val$cells <- val$cells %>% left_join(depth, by = "cellID") %>%
    mutate(ValidCell = ValidCell & UMIxcellID >= min.cell.UMI)
  val$pairs <- val$pairs %>% select(-any_of("ValidCell")) %>%
    left_join(val$cells %>% select(cellID, ValidCell), by = "cellID")
  complete_table_fin <- val$pairs %>%
    left_join(cu %>% select(cellID, UCI, barcode, shRNA, set, gene), by = c("cellID", "UCI"))
  utils::write.csv(complete_table_fin, op(sprintf("complete_table_fin_%s.csv", sample_name)),
                   row.names = FALSE)

  per_cell <- cu %>% group_by(cellID) %>% arrange(desc(reads), .by_group = TRUE) %>%
    summarise(top_barcode = first(barcode), top_UCI = first(UCI), top_shRNA = first(shRNA),
              top_set = first(set), gene = first(gene),
              reads_top = first(reads),
              reads_second = if (n() > 1) sort(reads, decreasing = TRUE)[2] else 0L,
              reads_total = sum(reads), n_UCI = n(), n_sets = n_distinct(set), .groups = "drop") %>%
    mutate(sublibrary = "all",
           pct_top = 100 * reads_top / reads_total,
           ratio_top = ifelse(reads_second > 0, reads_top / reads_second, Inf)) %>%
    left_join(val$cells %>% select(cellID, True_UCIxcellID, ValidCell), by = "cellID")
  utils::write.csv(per_cell, op("tenxcatch_per_cell_table.csv"), row.names = FALSE)

  uci_tab <- per_cell %>% count(n_UCI, name = "cells") %>%
    mutate(pct = round(100 * cells / sum(cells), 1))
  clone_summary <- assign_tbl %>% group_by(set, gene, shRNA, construct) %>%
    summarise(clones = n_distinct(UCI), cells = n_distinct(cellID), .groups = "drop") %>%
    arrange(desc(clones), desc(cells))
  utils::write.csv(clone_summary, op("tenxcatch_shRNA_clone_summary.csv"), row.names = FALSE)

  res <- list(sample = sample_name, output_dir = output_dir, chemistry = "tenx",
              params = list(reference = reference, empty.reference = empty.reference,
                            constant = constant, UCI.length = UL, barcode.length = BL,
                            cb.length = cb.length, umi.length = umi.length,
                            max.edit.distance = max.edit.distance,
                            percentage = percentage, ratio = ratio, UMI.count = UMI.count,
                            min.cell.UMI = min.cell.UMI, mode = mode, gene.sep = gene.sep,
                            set.names = set.names, set.sizes = set.sizes,
                            whitelisted = !is.null(wl), colors = cols),
              total_pairs = st$total_pairs, per_sublibrary = st$per_sublib,
              reference_table = ref,
              cell_class = cell_class, read_class = read_class, nocall = nocall,
              assignments = assign_tbl, empty_reads = st$empty,
              reads_assigned = n_reads_assigned,
              valid_cell_reads = st$valid_cell_reads,
              cu = cu, tab = tab, thresholds = thr,
              complete_table_fin = complete_table_fin, per_cell = per_cell,
              uci_table = uci_tab, clone_summary = clone_summary,
              n_cell_barcodes = n_distinct(cell_class$cell),
              n_cells_called = length(assigned_cells),
              n_valid = sum(val$cells$ValidCell),
              elapsed = (proc.time() - ptm)[["elapsed"]])
  class(res) <- c("tenxcatch", "pipcatch", "list")
  saveRDS(res, op(sprintf("tenxcatch_%s.rds", sample_name)))
  res$input_files <- .cat_input_manifest(
    list("fastq read 1" = f1, "fastq read 2" = f2, "barcode-to-gene reference" = refp,
         "called cells / whitelist" = if (is.character(cell.whitelist)) cell.whitelist else NULL,
         "called-cell list" = if (is.character(called.cells)) called.cells else NULL),
    sample_name, "tenxcatch", res$output_dir)

  if (plots) .cat_plots_qc(res, zoom.x = zoom.x, zoom.y = zoom.y, jitter = jitter)
  # Keep the resolved called-cell list on the object: the re-filter draws the
  # threshold-choice panel and without this it would fall back to every observed
  # barcode, ambient floor included, which overstates coverage.
  res$called_cells <- if (!is.null(called.cells)) .par_called(called.cells) else NULL
  if (plots) .par_plots_analysis(res, called = called.cells, prefix = "tenxcatch",
                                 verbose = verbose)
  if (verbose) print(res)
  invisible(res)
}

#' @export
print.tenxcatch <- function(x, ...) {
  cat(sprintf("<tenxcatch> sample %s | %s read pairs | %s cell barcodes%s\n",
              x$sample, format(x$total_pairs, big.mark = ","),
              format(x$n_cell_barcodes, big.mark = ","),
              if (isTRUE(x$params$whitelisted)) " (whitelisted)" else " (NO whitelist)"))
  cat(sprintf("  %s UMI-deduplicated events, %s barcodes with a call, %s ValidCell\n",
              format(nrow(x$assignments), big.mark = ","),
              format(x$n_cells_called, big.mark = ","),
              format(x$n_valid, big.mark = ",")))
  cat(sprintf("  %d constructs, %s clones | UMI.count %s, percentage %s, ratio %s, cell floor %s\n",
              nrow(x$clone_summary), format(sum(x$clone_summary$clones), big.mark = ","),
              x$params$UMI.count, x$params$percentage, x$params$ratio, x$params$min.cell.UMI))
  invisible(x)
}
# ── 10x downstream: thin, chemistry-named wrappers ──────────────────────────
#
# The re-filter, empty-cell and Seurat steps are chemistry-independent: they
# act on the (cell, UCI) table that every entry point produces, not on reads.
# They are implemented once in R/catcheR_parse_functions.R. These wrappers
# exist so a 10x script reads as 10x throughout rather than calling a
# function named after a different chemistry.

#' Re-filter a tenxcatch object at new thresholds
#'
#' Chemistry-named wrapper around [catcheR_parse_filtercatch()], which holds the
#' shared implementation. All arguments are passed through unchanged.
#' @inheritParams catcheR_parse_filtercatch
#' @return a `parse_filtercatch` object; files are written with a `tenxcatch_`
#'   prefix because the prefix follows the class of `x`.
#' @export
catcheR_tenx_filtercatch <- function(x, ...) catcheR_parse_filtercatch(x, ...)

#' Identify empty-vector (unperturbed) cells in a tenxcatch object
#'
#' Chemistry-named wrapper around [catcheR_parse_nocatch()].
#' @inheritParams catcheR_parse_nocatch
#' @export
catcheR_tenx_nocatch <- function(x, ...) catcheR_parse_nocatch(x, ...)

#' Join tenxcatch calls to a Cell Ranger expression matrix
#'
#' Chemistry-named wrapper around [catcheR_pip_seurat()], which reads a
#' MatrixMarket triple with `Seurat::ReadMtx()` - the format Cell Ranger and
#' PIPseeker both emit. Pass Cell Ranger's `matrix.mtx.gz` /
#' `features.tsv.gz` / `barcodes.tsv.gz` as `mtx` / `features` / `cells`.
#'
#' NOT routed through [catcheR_parse_seurat()]: that one reads split-pipe's
#' `all_genes.csv.gz` / `cell_metadata.csv.gz`, i.e. CSVs with `gene_name` and
#' `bc_wells` columns, which a Cell Ranger `features.tsv.gz` does not have.
#' @inheritParams catcheR_pip_seurat
#' @export
catcheR_tenx_seurat <- function(x, ...) catcheR_pip_seurat(x, ...)
