## Part of catcheR. Parse (Evercode) entry point. See
## inst/rmd/catcheR_parsecatch_analysis_EB.Rmd for the worked template.

# ── 1. catcheR_parsecatch ───────────────────────────────────────────────────

#' Extract CATCHER clone identities from a Parse (Evercode) barcode library
#'
#' The Parse equivalent of [catcheR_pipcatch()]. Streams each sublibrary's read
#' pair once, pulls the shRNA barcode and UCI off read 1, parses and error-
#' corrects the three-round cell barcode on read 2, deduplicates by UMI, and
#' produces the same tables and plots as the 10X and PIPseq paths.
#'
#' @section Why the barcode file is required:
#' `barcode.data` must be split-pipe's own `barcode_data.csv`. A whitelist
#' inferred from the data looks convincing and is wrong: round 1 carries two
#' sequences per well, one oligo-dT and one random hexamer, and the hexamer
#' barcodes are several-fold less abundant in a 3' enrichment library. Taking
#' the most abundant sequences per round recovers only the dT set and silently
#' discards or mis-corrects every hexamer-primed read.
#'
#' @section Sublibraries:
#' Parse reuses the same well combinations across sublibraries, so the same
#' 24 nt barcode in two sublibraries is two different cells. Pass one entry per
#' sublibrary in `fastq.read1`, `fastq.read2` and `sublibraries`; cell IDs are
#' emitted in split-pipe's `bc_wells` form (`01_01_59__s1`) so they join
#' directly to `cell_metadata.csv`.
#'
#' @param folder directory holding the FASTQ files and the reference CSVs.
#' @param fastq.read1,fastq.read2 file names, one per sublibrary. Read 1
#'   carries the cassette, read 2 the cell barcode - the opposite of PIPseq.
#' @param sublibraries suffixes matching split-pipe's, e.g. `c("s1","s2")`.
#' @param barcode.data split-pipe's `barcode_data.csv`.
#' @param barcodes.genes one or more barcode-to-gene CSVs. Orientation is
#'   resolved per file against the reads, so pooled libraries written in
#'   opposite orientations are handled without being told.
#' @param set.names,set.sizes labels and nominal sizes for those CSVs.
#' @param reference,empty.reference,constant construct sequences. Defaults are
#'   the CATCHER vector's.
#' @param UCI.length,barcode.length cassette field widths.
#' @param umi.length,linker.a,linker.b read-2 geometry. Defaults are Parse v3/v4.
#' @param max.edit.distance barcode correction radius, 1 (default) or 2.
#'   split-pipe corrects to 2; the ambiguous fraction grows with it.
#' @param percentage,ratio,UMI.count catcheR's calling thresholds. The Parse
#'   defaults differ from the 10X ones and are evidence-based: see
#'   [catcheR_call_accuracy()].
#' @param min.cell.UMI minimum cassette UMIs for a cell to be callable,
#'   INCLUSIVE (`UMIxcellID >= min.cell.UMI`). See the threshold section below;
#'   `catcheR_*_filtercatch()` defaults to the same 20, so it keeps applying.
#' @param called.cells optional: the cell list the chemistry's own pipeline
#'   called, as a path to `cell_metadata.csv[.gz]` (its `bc_wells` column holds
#'   the same cell ID this function emits), a plain file of IDs, or a character
#'   vector. Supplying it unlocks the panels judged against a called-cell
#'   denominator, which is the honest one: the barcodes this function observes
#'   include an ambient floor the chemistry never called as cells.
#'   Below roughly 20, one call in ten is of the wrong clone.
#' @param mode,gene.sep,zoom.x,zoom.y,jitter,chunk.reads,max.reads,plots,verbose
#'   as in [catcheR_pipcatch()].
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
#'   `<output_dir>/parsecatch/`, so the three chemistries cannot
#'   overwrite each other's plots. Set FALSE for the flat layout.
#'
#' @return an object of class `parsecatch`, also saved as an `.rds`.
#' @export
catcheR_parsecatch <- function(folder,
                               fastq.read1,
                               fastq.read2,
                               sublibraries    = NULL,
                               barcode.data    = "barcode_data.csv",
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
                               umi.length      = 10,
                               linker.a        = "ATGAGGGGTCAG",
                               linker.b        = "TCCAACCACCTC",
                               max.edit.distance = 1,
                               percentage      = 15,
                               ratio           = 2,
                               UMI.count       = 1,
                               min.cell.UMI    = 20,
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
  if (length(fastq.read1) != length(fastq.read2))
    stop("fastq.read1 and fastq.read2 must have the same length (one per sublibrary).",
         call. = FALSE)
  if (is.null(sublibraries)) sublibraries <- paste0("s", seq_along(fastq.read1))
  if (length(sublibraries) != length(fastq.read1))
    stop("sublibraries must have one entry per FASTQ pair.", call. = FALSE)
  f1 <- file.path(folder, fastq.read1); f2 <- file.path(folder, fastq.read2)
  refp <- file.path(folder, barcodes.genes); bdp <- file.path(folder, barcode.data)
  for (f in c(f1, f2, refp, bdp))
    if (!file.exists(f)) stop("Missing file: ", f, call. = FALSE)
  if (is.null(set.names))
    set.names <- if (length(refp) == 1) "set1" else paste0("set", seq_along(refp))
  # One subfolder per chemistry. All three pipelines write the 13 plots that
  # carry catcheR's own file names (barcode_distribution, UMIxUCI, ...), which
  # have no chemistry in the filename by design, so running 10x and Parse into
  # the same Output/ used to have the second silently replace the first.
  # Everything downstream reads `x$output_dir`, so nesting once here propagates.
  if (isTRUE(subdir)) output_dir <- file.path(output_dir, "parsecatch")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  op <- function(f) file.path(output_dir, f)

  ## whitelists --------------------------------------------------------------
  wl <- .par_whitelist(bdp)
  cmaps <- lapply(wl$maps, .par_corr_map, max.edit.distance = as.integer(max.edit.distance))
  if (verbose) {
    message(sprintf("Parse barcodes: %s", paste(sub("^#\\s*", "", wl$header[1]), collapse = "")))
    for (r in names(wl$maps))
      message(sprintf("  %s: %d whitelist sequences over %d wells -> %d accepted with correction",
                      r, length(wl$maps[[r]]), length(unique(wl$maps[[r]])), length(cmaps[[r]])))
  }

  RCREF    <- .cat_rc(reference)
  RC_EMPTY <- .cat_rc(empty.reference)
  ANCH     <- substr(RCREF, 1, nchar(RCREF) - barcode.length)
  UL <- UCI.length; BL <- barcode.length

  ## stream each sublibrary --------------------------------------------------
  parts <- lapply(seq_along(f1), function(i) {
    if (verbose) message(sprintf("Streaming sublibrary %s ...", sublibraries[i]))
    ad <- .par_adapter(cmaps, sublibraries[i], linker.a, linker.b, umi.length, BL)
    s <- .cat_stream(f1[i], f2[i], ad, ANCH, RC_EMPTY, constant, UL, BL,
                     chunk.reads, max.reads, verbose)
    s$sublibrary <- sublibraries[i]
    if (nrow(s$cassette)) s$cassette$sublibrary <- sublibraries[i]
    s
  })
  st <- list(
    total_pairs  = sum(vapply(parts, function(p) p$total_pairs, numeric(1))),
    class_counts = tapply(unlist(lapply(parts, function(p) unname(p$class_counts))),
                          unlist(lapply(parts, function(p) names(p$class_counts))), sum),
    cell_class   = bind_rows(lapply(parts, `[[`, "cell_class")) %>%
                   count(cell, cls, wt = n, name = "n"),
    cassette     = bind_rows(lapply(parts, `[[`, "cassette")),
    empty        = bind_rows(lapply(parts, `[[`, "empty")),
    valid_cell_reads = sum(vapply(parts, function(p) p$valid_cell_reads, numeric(1))),
    per_sublib   = data.frame(sublibrary = sublibraries,
                              pairs = vapply(parts, function(p) p$total_pairs, numeric(1)),
                              cell_bc_reads = vapply(parts, function(p) p$valid_cell_reads, numeric(1)),
                              cassette_reads = vapply(parts, function(p) nrow(p$cassette), integer(1))))
  if (!nrow(st$cassette))
    stop("No cassette-bearing reads found. Check `reference`, `constant` and the read order ",
         "(Parse carries the cassette on read 1).", call. = FALSE)
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

  ## UMI deduplication -------------------------------------------------------
  ## This is the substantive difference from PIPseq: `index` here is a real
  ## 10 nt UMI, so one row per MOLECULE rather than per read. Everything the
  ## shared layer calls UMIxUCI is then genuinely a UMI count.
  assign_tbl <- T_reads %>% filter(!is.na(shRNA)) %>%
    select(cellID, index, barcode, UCI, shRNA, set, gene, construct) %>%
    distinct()
  n_reads_assigned <- sum(!is.na(T_reads$shRNA))
  ## Reads supporting each molecule, carried as a column on the molecule table.
  ## The dedup above collapses them, so the tally has to be taken here or the
  ## rarefaction curve cannot be built: it needs to know, molecule by molecule,
  ## which ones would drop out at reduced depth and which cell they belong to.
  mol_reads <- T_reads %>% filter(!is.na(shRNA)) %>%
    count(cellID, index, barcode, UCI, name = "reads")
  assign_tbl <- assign_tbl %>%
    left_join(mol_reads, by = c("cellID", "index", "barcode", "UCI"))
  utils::write.csv(assign_tbl, op("parsecatch_cell_assignments.csv"), row.names = FALSE)

  ## read-class summary ------------------------------------------------------
  assigned_cells <- unique(assign_tbl$cellID)
  split_cass <- c("assigned" = n_reads_assigned,
                  "cassette present, barcode not in reference" = sum(is.na(T_reads$shRNA)))
  read_class <- bind_rows(
    data.frame(class = names(st$class_counts), reads = as.integer(unname(st$class_counts))) %>%
      filter(class != "cassette"),
    data.frame(class = names(split_cass), reads = as.integer(unname(split_cass)))
  ) %>% arrange(desc(reads)) %>% mutate(pct = round(100 * reads / sum(reads), 2))
  utils::write.csv(read_class, op("parsecatch_read_class_summary.csv"), row.names = FALSE)

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
  ## the per-cell UMI floor is applied here, not inside the shared validator,
  ## because it is a Parse-depth decision rather than part of catcheR's rule.
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
    mutate(sublibrary = sub("^.*__", "", cellID),
           pct_top = 100 * reads_top / reads_total,
           ratio_top = ifelse(reads_second > 0, reads_top / reads_second, Inf)) %>%
    left_join(val$cells %>% select(cellID, True_UCIxcellID, ValidCell), by = "cellID")
  utils::write.csv(per_cell, op("parsecatch_per_cell_table.csv"), row.names = FALSE)
  uci_tab <- per_cell %>% count(n_UCI, name = "cells") %>%
    mutate(pct = round(100 * cells / sum(cells), 1))
  clone_summary <- assign_tbl %>% group_by(set, gene, shRNA, construct) %>%
    summarise(clones = n_distinct(UCI), cells = n_distinct(cellID), .groups = "drop") %>%
    arrange(desc(clones), desc(cells))
  utils::write.csv(clone_summary, op("parsecatch_shRNA_clone_summary.csv"), row.names = FALSE)

  res <- list(sample = sample_name, output_dir = output_dir, chemistry = "parse",
              params = list(reference = reference, empty.reference = empty.reference,
                            constant = constant, UCI.length = UL, barcode.length = BL,
                            umi.length = umi.length, linker.a = linker.a, linker.b = linker.b,
                            max.edit.distance = max.edit.distance,
                            percentage = percentage, ratio = ratio, UMI.count = UMI.count,
                            min.cell.UMI = min.cell.UMI, mode = mode, gene.sep = gene.sep,
                            set.names = set.names, set.sizes = set.sizes,
                            sublibraries = sublibraries, colors = cols),
              total_pairs = st$total_pairs, per_sublibrary = st$per_sublib,
              reference_table = ref, whitelist = wl$table,
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
  class(res) <- c("parsecatch", "pipcatch", "list")
  saveRDS(res, op(sprintf("parsecatch_%s.rds", sample_name)))
  res$input_files <- .cat_input_manifest(
    list("fastq read 1" = f1, "fastq read 2" = f2, "barcode-to-gene reference" = refp,
         "barcode round definitions" = bdp,
         "called-cell list" = if (is.character(called.cells)) called.cells else NULL),
    sample_name, "parsecatch", res$output_dir)

  if (plots) .cat_plots_qc(res, zoom.x = zoom.x, zoom.y = zoom.y, jitter = jitter)
  # Keep the resolved called-cell list on the object: the re-filter draws the
  # threshold-choice panel and without this it would fall back to every observed
  # barcode, ambient floor included, which overstates coverage.
  res$called_cells <- if (!is.null(called.cells)) .par_called(called.cells) else NULL
  if (plots) .par_plots_analysis(res, called = called.cells, verbose = verbose)
  if (verbose) print(res)
  invisible(res)
}

#' @export
print.parsecatch <- function(x, ...) {
  cat(sprintf("<parsecatch> sample %s | %s read pairs over %d sublibrar%s | %s cell barcodes\n",
              x$sample, format(x$total_pairs, big.mark = ","),
              nrow(x$per_sublibrary), if (nrow(x$per_sublibrary) == 1) "y" else "ies",
              format(x$n_cell_barcodes, big.mark = ",")))
  cat(sprintf("  %s UMI-deduplicated events, %s barcodes with a call, %s ValidCell\n",
              format(nrow(x$assignments), big.mark = ","),
              format(x$n_cells_called, big.mark = ","),
              format(x$n_valid, big.mark = ",")))
  cat(sprintf("  %d constructs, %s clones | UMI.count %s, percentage %s, ratio %s, cell floor %s\n",
              nrow(x$clone_summary), format(sum(x$clone_summary$clones), big.mark = ","),
              x$params$UMI.count, x$params$percentage, x$params$ratio, x$params$min.cell.UMI))
  invisible(x)
}
