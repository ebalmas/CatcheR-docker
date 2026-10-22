## catcheR_pip_functions.R ---------------------------------------------------
## PIPseq (Fluent) port of the catcheR 10X entry points.
##
##   catcheR_pipcatch()        extraction + QC + threshold suggestion
##   catcheR_pip_filtercatch() re-filter with refined thresholds
##   catcheR_pip_nocatch()     empty / no-insertion cells
##   catcheR_pip_seurat()      join calls to the PIPseeker matrix, build Seurat
##   catcheR_pip_report()      render the step-by-step Rmd report
##
## Differences from the 10X pipeline, all forced by the chemistry:
##   * No docker and no bash scripts: everything runs in R.
##   * R1 is not sliced at fixed offsets. PIPseq puts 0-3 phasing bases before a
##     4-tier barcode (8 ATG 6 GAG 6 TCGAG 8), so the fields move read to read.
##   * The 3 nt molecular index is NOT a UMI (64 values, all observed). Every
##     metric counts READS. `UMIxUCI` is kept as a column name for continuity
##     with catcheR but holds a read count.
##   * The barcode/gene CSVs are orientation-checked against the reads instead
##     of being trusted as given.
## ---------------------------------------------------------------------------

## Part of catcheR. See R/catcheR_pipcatch.R for the entry point and
## inst/rmd/catcheR_pipcatch_report.Rmd for the report template.

#' @importFrom dplyr %>% filter mutate group_by summarise arrange bind_rows
#'   case_when n_distinct pull rename ungroup count left_join select distinct
#'   transmute slice first nth n desc any_of all_of
#' @importFrom tidyr replace_na
#' @importFrom ggplot2 ggplot aes geom_col geom_text geom_tile geom_point
#'   geom_histogram geom_vline facet_grid labs coord_flip position_stack
#'   scale_fill_manual scale_colour_manual scale_fill_gradient scale_size_area
#'   scale_x_continuous scale_y_log10 expansion theme theme_set theme_minimal
#'   element_text element_blank ggsave
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom utils read.csv write.csv write.table
#' @importFrom stats setNames reorder median
#' @name catcheR-pipcatch-internals
#' @keywords internal
NULL

# ── internals ───────────────────────────────────────────────────────────────

.cat_rc <- function(x) vapply(strsplit(chartr("ACGTN", "TGCAN", x), NULL),
                              function(z) paste(rev(z), collapse = ""), character(1))

.cat_theme <- function() theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 9, hjust = 0))

#' Parse a PIPseq read 1 into its 4-tier cell barcode.
#'
#' `valid` requires the three linkers AND the downstream 3 nt index, the
#' constant ACG and the oligo-dT. Matching the linkers alone admits ~15x more
#' reads, nearly all of them spurious hits in construct or adapter sequence.
.pip_parse_r1 <- function(s1) {
  n <- length(s1); phase <- rep(NA_integer_, n)
  for (p in 0:3) {
    hit <- is.na(phase) &
      substr(s1, p +  9, p + 11) == "ATG" &
      substr(s1, p + 18, p + 20) == "GAG" &
      substr(s1, p + 27, p + 31) == "TCGAG"
    phase[which(hit)] <- p
  }
  ok <- !is.na(phase); ph <- ifelse(ok, phase, 0L)
  dT <- substr(s1, ph + 46, ph + 53)
  list(linker_only = ok,
       valid = ok & substr(s1, ph + 43, ph + 45) == "ACG" &
               nchar(gsub("[^T]", "", dT)) >= 7,
       cell  = paste0(substr(s1, ph +  1, ph +  8), substr(s1, ph + 12, ph + 17),
                      substr(s1, ph + 21, ph + 26), substr(s1, ph + 32, ph + 39)),
       index = substr(s1, ph + 40, ph + 42))
}

#' PIPseq adapter for [.cat_stream()].
#'
#' Supplies the two chemistry-specific pieces the shared streamer needs: which
#' read carries the cell barcode, and how to parse it. The cassette side is
#' chemistry-independent - in every chemistry so far the UCI and the shRNA
#' barcode sit immediately 5' of the vector constant, in reverse complement -
#' so the streamer slices it the same way for all of them.
.pip_adapter <- function() {
  list(name      = "pipseq",
       cell_read = 1L,       # cell barcode on read 1, cassette on read 2
       has_umi   = FALSE,    # 3 nt index, not a deduplicable UMI
       lead      = 0L,       # cassette starts at the very start of read 2
       cell_fun  = function(s) {
         pr <- .pip_parse_r1(s)
         list(valid = pr$valid, cell = pr$cell, umi = pr$index)
       })
}

#' Load the barcode/gene CSVs and resolve each file's orientation against the
#' reads. Pooled experiments are routinely written in opposite orientations.
.cat_load_reference <- function(paths, set_names, observed, gene.sep = "_") {
  refs <- lapply(seq_along(paths), function(i) {
    r <- utils::read.csv(paths[i], header = FALSE, col.names = c("barcode", "shRNA"),
                         strip.white = TRUE, colClasses = "character")
    r$barcode <- toupper(trimws(r$barcode)); r$shRNA <- trimws(r$shRNA)
    hf <- sum(r$barcode %in% observed); hr <- sum(.cat_rc(r$barcode) %in% observed)
    if (hr > hf) { r$barcode <- .cat_rc(r$barcode); or <- "reverse-complemented" } else or <- "as given"
    message(sprintf("  %s: %d barcodes, %s (fwd %d / rc %d)", set_names[i], nrow(r), or, hf, hr))
    r$set <- set_names[i]; r$csv_orientation <- or; r
  })
  ref <- bind_rows(refs)
  if (any(duplicated(ref$barcode)))
    stop("Barcode collision after merging the reference files: ",
         paste(ref$barcode[duplicated(ref$barcode)], collapse = ", "), call. = FALSE)
  # NOTE: catcheR's remove_after_dot() assumes shRNA names look like GENE.n.
  # These CSVs use GENE_n, so the separator is a parameter here.
  ref %>% mutate(
    gene = sub(paste0("\\", gene.sep, "[^", gene.sep, "]*$"), "", shRNA),
    idx  = suppressWarnings(as.integer(sub(paste0(".*\\", gene.sep), "", shRNA))),
    construct = ifelse(shRNA %in% shRNA[duplicated(shRNA)],
                       paste0(shRNA, " (", set, ")"), shRNA))
}

#' Collapse reads to one row per (cellID, UCI), the catcheR unit of evidence.
.cat_tabulate <- function(assign_tbl) {
  cu <- assign_tbl %>% count(cellID, barcode, UCI, shRNA, set, gene, name = "reads")
  tab <- cu %>% group_by(cellID, UCI) %>%
    summarise(UMIxUCI = sum(reads), .groups = "drop") %>%
    left_join(cu %>% count(cellID, wt = reads, name = "UMIxcellID"), by = "cellID") %>%
    mutate(UMIpercentagexUCI = 100 * UMIxUCI / UMIxcellID)
  list(cu = cu, tab = tab)
}

#' catcheR's two automatic threshold estimators, with the guard the original
#' lacks. `noise` divides by (1 - perc_1/100 * 1.35); once more than
#' 1/1.35 = 74.07% of UCIs are singletons that term goes negative and the
#' original subsets with a negative index, which errors in R.
.cat_thresholds <- function(tab) {
  U <- nrow(tab); one <- sum(tab$UMIxUCI == 1); perc_1 <- 100 * one / U
  p <- perc_1 / 100 * 1.35
  thr_noise <- if (p >= 1) NA_integer_ else
    min(sort(tab$UMIxUCI, decreasing = TRUE)[seq_len(ceiling(U * (1 - p)))])
  h <- as.data.frame(table(tab$UMIxUCI), stringsAsFactors = FALSE)
  names(h) <- c("level", "count"); h$level <- as.integer(h$level); v <- h$count
  val <- if (length(v) >= 3)
    which(c(FALSE, v[-c(1, length(v))] < v[-c(length(v) - 1, length(v))] &
                   v[-c(1, length(v))] <= v[-c(1, 2)], FALSE)) else integer(0)
  # catcheR uses findValleys(...)[1] directly, i.e. the histogram POSITION.
  # That equals the count only when levels are contiguous. Read the level.
  list(hist = h, pairs = U, singletons = one, pct_singletons = perc_1,
       noise = thr_noise,
       bimodal = if (length(val)) h$level[val[1]] else NA_integer_)
}

#' catcheR's per-cell validity rule, verbatim:
#'   TorF            = UMIpercentagexUCI > percentage & UMIxUCI > UMI.count
#'   True_UCIxcellID = number of TorF UCIs in the cell
#'   ValidCell       = True_UCIxcellID == 1 AND (top/second >= ratio OR pct == 100)
#' The ratio is consulted only after a cell is already down to one surviving
#' UCI; it is not an alternative to being single-UCI.
.cat_validate <- function(tab, UMI.count, percentage, ratio) {
  t <- tab %>% mutate(TorF = UMIpercentagexUCI > percentage & UMIxUCI > UMI.count) %>%
    group_by(cellID) %>% mutate(True_UCIxcellID = sum(TorF)) %>% ungroup()
  v <- t %>% group_by(cellID) %>%
    summarise(True_UCIxcellID = first(True_UCIxcellID),
              mx = max(UMIxUCI),
              second = if (n() > 1) sort(UMIxUCI, decreasing = TRUE)[2] else NA_real_,
              pct_at_max = UMIpercentagexUCI[which.max(UMIxUCI)], .groups = "drop") %>%
    mutate(ValidCell = True_UCIxcellID == 1 &
             ((!is.na(second) & mx / second >= ratio) | pct_at_max == 100))
  list(pairs = t %>% left_join(v %>% select(cellID, ValidCell), by = "cellID"),
       cells = v)
}

# ── plotting ────────────────────────────────────────────────────────────────
# Two families of output:
#   * the 13 plots catcheR_10Xcatch produces, written under catcheR's own file
#     names (.jpg + .pdf) so they can be referenced exactly as in the 10X
#     figures;
#   * the additional PIPseq-specific panels (cell detection, no-call fates,
#     clone matrix) that have no 10X equivalent, written as .png.
#
# `jitter` controls the 2D panels. catcheR uses geom_count(), which sizes each
# marker by how many (cell, UCI) pairs share the coordinate. That is the right
# geometry at 10X depth, where UMIxUCI spans 1-200+ and the percentages are
# near-continuous, so coordinates rarely collide. At PIPseq enrichment depth
# UMIxUCI takes a handful of integer values and the percentage only simple
# fractions, so thousands of pairs collapse onto a few dozen coordinates and
# geom_count degenerates into three or four huge bubbles. `jitter = "auto"`
# detects that case and spreads the pairs instead, which is what reveals the
# real density. The data plotted is identical either way.

.cat_save_both <- function(p, dir, name, w = 7, h = 9) {
  ggplot2::ggsave(file.path(dir, paste0(name, ".jpg")), plot = p,
                  width = w, height = h, dpi = 300)
  grDevices::pdf(file.path(dir, paste0(name, ".pdf")), width = w, height = h)
  print(p); grDevices::dev.off()
  invisible(NULL)
}

#' @return TRUE when the (percentage, UMIxUCI) plane is so sparsely populated
#'   that geom_count collapses the data into a few markers.
.cat_jitter_on <- function(tab, jitter = c("auto", "TRUE", "FALSE")) {
  if (isTRUE(jitter)) return(TRUE)
  if (identical(jitter, FALSE)) return(FALSE)
  jitter <- match.arg(as.character(jitter)[1], c("auto", "TRUE", "FALSE"))
  if (jitter == "TRUE")  return(TRUE)
  if (jitter == "FALSE") return(FALSE)
  # "auto" now always answers FALSE: the aggregated-bubble geometry handles the
  # sparse-coordinate case correctly on its own, so there is no depth at which
  # jittering is the better default. Kept as a function rather than deleted so
  # jitter = TRUE still works for anyone who wants the point cloud.
  FALSE
}

#' catcheR's 2D percentage-vs-UMIxUCI panel.
.cat_2d <- function(tab, percentage, UMI.count, jitter = FALSE,
                    colour_var = NULL, pal = NULL, alpha = 0.3,
                    thresholds = TRUE) {
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = UMIpercentagexUCI, y = UMIxUCI))
  if (!is.null(colour_var))
    p <- p + ggplot2::aes(colour = factor(.data[[colour_var]]))
  # Aggregated bubbles are the default everywhere. catcheR uses geom_count();
  # its default size scale makes a handful of enormous markers once a few
  # thousand pairs share one coordinate, so the radius is put on an area scale
  # capped at 9 pt, which is what the multi-panel PNG figures already did. The
  # two code paths used to disagree - the .jpg panels jittered while the PNG
  # panels aggregated - and this is the single geometry both now use.
  #
  # `jitter = TRUE` remains available but is no longer selected automatically:
  # on a log y axis the vertical jitter smears the enormous y = 1 and y = 2
  # stacks into solid opaque blocks, which hides exactly the density it is
  # meant to reveal.
  # a list, not `a + b`: two ggproto objects cannot be added to each other,
  # only to a plot
  p <- p + if (jitter)
    list(ggplot2::geom_jitter(alpha = alpha, size = 0.35, width = 2.2, height = 0.03))
  else
    list(ggplot2::geom_count(alpha = 0.55),
         ggplot2::scale_size_area(max_size = 9, guide = "none"))
  p <- p + ggplot2::scale_y_continuous(trans = "log10") + .cat_theme()
  if (thresholds)
    p <- p + ggplot2::geom_vline(xintercept = percentage, colour = "red") +
             ggplot2::geom_hline(yintercept = UMI.count, colour = "red")
  if (!is.null(pal))
    p <- p + ggplot2::scale_colour_manual(values = pal, name = colour_var)
  else if (!is.null(colour_var))
    p <- p + ggplot2::labs(colour = colour_var)
  if (jitter)
    p <- p + ggplot2::labs(caption = paste("jitter = TRUE: one marker per (cell, UCI) pair,",
                                           "positions spread; the default aggregates instead"))
  p
}

# ── the 7 plots catcheR writes at the explorative-analysis stage ─────────────

.cat_plots_stage1 <- function(x, dir, zoom.x = NULL, zoom.y = NULL, jitter = "auto") {
  a <- x$assignments; tab <- x$tab; thr <- x$thresholds; p <- x$params
  ggplot2::theme_set(.cat_theme())
  jit <- .cat_jitter_on(tab, jitter)

  # 1-2. barcode_distribution / gene_distribution. catcheR keeps the categories
  #      in alphabetical order and fills by category; an optional colors.csv in
  #      the input folder overrides the palette.
  bd <- dplyr::count(a, construct, name = "reads")
  gd <- dplyr::count(a, gene, name = "reads")
  rot <- ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = .5, hjust = 1))
  pb <- ggplot2::ggplot(bd, ggplot2::aes(construct, reads, fill = construct)) +
    ggplot2::geom_bar(stat = "identity", show.legend = FALSE) +
    ggplot2::labs(x = NULL, y = "reads") + rot +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 4))
  if (!is.null(p$colors) && nrow(p$colors))
    pb <- pb + ggplot2::scale_fill_manual(
      values = stats::setNames(p$colors$colors, p$colors$shRNA), na.value = "grey70")
  .cat_save_both(pb, dir, "barcode_distribution", 7, 9)
  pg <- ggplot2::ggplot(gd, ggplot2::aes(gene, reads, fill = gene)) +
    ggplot2::geom_bar(stat = "identity", show.legend = FALSE) +
    ggplot2::labs(x = NULL, y = "reads") + rot
  .cat_save_both(pg, dir, "gene_distribution", 7, 9)

  # 3-4. UMIxUCI histogram, full and zoomed
  pu <- ggplot2::ggplot(tab, ggplot2::aes(x = UMIxUCI)) +
    ggplot2::geom_histogram(binwidth = 1, alpha = .7) +
    ggplot2::labs(x = "UMIsxUCI (reads per cell-UCI)") +
    { if (!is.na(thr$bimodal)) ggplot2::geom_vline(xintercept = thr$bimodal, col = "red") }
  .cat_save_both(pu, dir, "UMIxUCI", 7, 9)
  cnt <- sort(thr$hist$count, decreasing = TRUE)
  # The point of the zoomed panel is to clip the singleton peak (y) and the
  # long sparse tail (x) so the body of the distribution is readable. Using
  # max(UMIxUCI) for x is fine when the range is 1-11, as on a read-counted
  # PIPseq library, but on a UMI-deduplicated Parse library the maximum is in
  # the thousands and the "zoom" reproduces the unzoomed panel. Cut at the
  # 99th percentile instead, with a floor so a shallow library still shows a
  # few bars.
  if (is.null(zoom.x))
    zoom.x <- max(10, ceiling(stats::quantile(tab$UMIxUCI, 0.99, names = FALSE)))
  if (is.null(zoom.y)) zoom.y <- ceiling(if (length(cnt) > 1) cnt[2] * 1.2 else cnt[1])
  # NOTE catcheR names the zoomed PDF with `y` twice (UMIxUCI_<y>_<y>.pdf) while
  # the JPG uses <x>_<y>; both are written as <x>_<y> here.
  .cat_save_both(pu + ggplot2::coord_cartesian(xlim = c(0, zoom.x), ylim = c(0, zoom.y)),
                 dir, sprintf("UMIxUCI_%s_%s", zoom.x, zoom.y), 7, 9)

  # 5. UCIxcell
  uc <- dplyr::count(dplyr::distinct(a, cellID, UCI), cellID, name = "Freq")
  pc <- ggplot2::ggplot(uc, ggplot2::aes(x = Freq)) +
    ggplot2::geom_histogram(bins = 70) + ggplot2::labs(x = "UCIsxcell")
  .cat_save_both(pc, dir, "UCIxcell", 7, 9)

  # 6. percentage_of_UMIxUCI_dist. catcheR de-duplicates on (UCI, percentage)
  #    and so drops cellID, merging distinct cells that share a UCI at the same
  #    percentage; every (cell, UCI) pair is kept here instead.
  pp <- ggplot2::ggplot(tab, ggplot2::aes(x = UMIpercentagexUCI)) +
    ggplot2::geom_histogram(bins = 100) +
    ggplot2::geom_vline(xintercept = p$percentage, col = "red") +
    ggplot2::labs(x = "UMIpercentagexUCI")
  .cat_save_both(pp, dir, "percentage_of_UMIxUCI_dist", 7, 9)

  # 7. 2D, no thresholds drawn
  .cat_save_both(.cat_2d(tab, p$percentage, p$UMI.count, jit, alpha = .3,
                         thresholds = FALSE),
                 dir, "2D_percentage_of_UMIxUCI_UMI_count", 10, 6)
  invisible(jit)
}

# ── the 6 plots catcheR writes at the cell-filtering stage ───────────────────

.cat_plots_stage2 <- function(pairs, cells, percentage, UMI.count, dir, jitter = "auto") {
  ggplot2::theme_set(.cat_theme())
  jit <- .cat_jitter_on(pairs, jitter)

  # 8. 2D with both thresholds drawn
  .cat_save_both(.cat_2d(pairs, percentage, UMI.count, jit, alpha = .2),
                 dir, "2D_percentage_of_UMIxUCI_UMI_count_thresholds", 10, 6)

  # 9. true_UCIxcell  (+ the table catcheR writes alongside it)
  tu <- dplyr::count(cells, True_UCIxcellID, name = "cells")
  utils::write.csv(tu, file.path(dir, "table_true_UCIxcell.csv"), row.names = FALSE)
  .cat_save_both(ggplot2::ggplot(tu, ggplot2::aes(factor(True_UCIxcellID), cells)) +
                   ggplot2::geom_col() +
                   ggplot2::geom_text(ggplot2::aes(label = cells), vjust = -.3, size = 2.6) +
                   ggplot2::labs(x = "True_UCIxcellID", y = "cells"),
                 dir, "true_UCIxcell", 5, 3)

  # 10. 2D coloured by True_UCIxcellID (catcheR palette: 0 grey, 1 orange, >=2 blue)
  pal_t <- c(`0` = "#C0C0C0", `1` = "#E69F00", `2` = "#0072B2", `3` = "#0072B2",
             `4` = "#0072B2", `5` = "#0072B2", `6` = "#0072B2")
  .cat_save_both(.cat_2d(pairs, percentage, UMI.count, jit,
                         colour_var = "True_UCIxcellID", pal = pal_t, alpha = .3),
                 dir, "2D_percentage_of_UMIxUCI_UMI_count_trueorfalse", 10, 6)

  # 11-12. cellsxclone, raw and log10
  cx <- pairs %>% dplyr::filter(TorF) %>% dplyr::distinct(cellID, UCI) %>%
    dplyr::count(UCI, name = "cellsxUCI")
  .cat_save_both(ggplot2::ggplot(cx, ggplot2::aes(x = cellsxUCI)) + ggplot2::geom_bar() +
                   ggplot2::labs(x = "cells per clone (UCI)", y = "clones"),
                 dir, "cellsxclone", 5, 3)
  .cat_save_both(ggplot2::ggplot(cx, ggplot2::aes(x = log10(cellsxUCI))) + ggplot2::geom_bar() +
                   ggplot2::labs(x = "log10(cells per clone)", y = "clones"),
                 dir, "cellsxclone_log10", 5, 3)

  # 13. 2D coloured by ValidCell (catcheR palette: FALSE blue, TRUE orange)
  .cat_save_both(.cat_2d(pairs, percentage, UMI.count, jit,
                         colour_var = "ValidCell",
                         pal = c(`FALSE` = "#0072B2", `TRUE` = "#E69F00"), alpha = .3),
                 dir, "2D_percentage_of_UMIxUCI_UMI_count_ValidCells", 10, 6)
  invisible(list(jitter = jit, cellsxclone = cx, true_UCIxcell = tu))
}

.cat_plots_qc <- function(x, zoom.x = NULL, zoom.y = NULL, jitter = "auto",
                          prefix = NULL) {
  # PNG panels are named for the chemistry that produced them, so a Parse run
  # does not litter its Output/ with files called "pipcatch_*".
  if (is.null(prefix))
    prefix <- if (inherits(x, "tenxcatch")) "tenxcatch"
              else if (inherits(x, "parsecatch")) "parsecatch" else "pipcatch"
  # Regenerating panels from a saved object is a normal thing to do, and the
  # entry point that created the directory is not in the picture then.
  dir.create(x$output_dir, showWarnings = FALSE, recursive = TRUE)
  op <- function(f) file.path(x$output_dir, sub("^pipcatch_", paste0(prefix, "_"), f))

  # The 13 plots that carry catcheR's own file names have no sample in the
  # name - that is deliberate, so they can be referenced exactly as in the
  # published figures - which means a second sample written into the SAME
  # output_dir silently overwrites the first. The .rds and the tables do carry
  # the sample, so the mismatch is easy to miss. Refuse to do it quietly.
  prior <- list.files(x$output_dir, pattern = "catch_[^.]+\\.rds$")
  prior <- setdiff(sub("^.*catch_(.*)\\.rds$", "\\1", prior), x$sample)
  if (length(prior))
    warning("output_dir already holds results for sample(s) ",
            paste(prior, collapse = ", "), ". The plots written under catcheR's ",
            "own names (barcode_distribution, UMIxUCI, true_UCIxcell, ...) carry no ",
            "sample in the filename and are about to be overwritten. Use one ",
            "output_dir per sample.", call. = FALSE, immediate. = TRUE)
  ggplot2::theme_set(.cat_theme())
  # the 13 plots catcheR_10Xcatch produces, under catcheR's own file names
  .cat_plots_stage1(x, x$output_dir, zoom.x, zoom.y, jitter)
  .cat_plots_stage2(x$complete_table_fin, x$per_cell, x$params$percentage,
                    x$params$UMI.count, x$output_dir, jitter)
  den <- x$n_cell_barcodes; asn <- x$n_cells_called
  pc <- x$per_cell; p <- x$params

  ## detection + library of origin
  d1 <- data.frame(grp = factor(c("shRNA barcode assigned", "no shRNA call"),
                                    levels = c("shRNA barcode assigned", "no shRNA call")),
                       cells = c(asn, den - asn))
  pD <- ggplot(d1, aes(1, cells, fill = grp)) + geom_col(width = .45) +
    geom_text(aes(label = sprintf("%s\n%.1f%%", format(cells, big.mark = ","), 100 * cells / den)),
              position = position_stack(vjust = .5), size = 2.6, colour = "white") +
    scale_fill_manual(values = c("#2E6E8E", "#9AA5AD"), name = NULL) + coord_flip() +
    labs(x = NULL, y = "validated PIPseq cell barcodes in read 1",
         title = sprintf("%s of %s cell barcodes carry an shRNA call",
                         format(asn, big.mark = ","), format(den, big.mark = ","))) +
    theme(axis.text.y = element_blank(), legend.position = "top")
  setc <- x$assignments %>% distinct(cellID, set) %>% count(set, name = "cells")
  if (!is.null(p$set.sizes)) setc <- setc %>%
    left_join(data.frame(set = p$set.names, n_shRNA = p$set.sizes), by = "set")
  pE <- ggplot(setc, aes(set, cells, fill = set)) + geom_col(width = .55, show.legend = FALSE) +
    geom_text(aes(label = if ("n_shRNA" %in% names(setc))
                    sprintf("%s (%.0f%%)\n%.1f per shRNA", format(cells, big.mark = ","),
                            100 * cells / sum(setc$cells), cells / n_shRNA)
                  else format(cells, big.mark = ",")), vjust = -.25, size = 2.5) +
    scale_fill_manual(values = rep(c("#2E6E8E", "#D1A13A"), length.out = nrow(setc))) +
    scale_y_continuous(expand = expansion(mult = c(0, .28))) +
    labs(x = NULL, y = "cells with an shRNA call", title = "Pooled library of origin")
  ggplot2::ggsave(op("pipcatch_cell_detection.png"),
         patchwork::wrap_plots(pD, pE, widths = c(1.5, 1)),
         width = 9.6, height = 2.9, dpi = 300)

  ## UCIs per cell + no-call fates
  # One bar per distinct value is readable when n_UCI tops out around 10, as it
  # does on read-counted enrichment data. On a UMI-deduplicated library the
  # ambient floor pushes it into the hundreds, so every level becomes an
  # unreadable tick with an overlapping label. Collapse the tail into one
  # bucket; the bars that carry the signal are all at the low end anyway.
  CAP <- 10L
  ut <- x$uci_table
  ut$grp <- ifelse(ut$n_UCI >= CAP, sprintf("\u2265%d", CAP), as.character(ut$n_UCI))
  ut <- stats::aggregate(list(cells = ut$cells), by = list(grp = ut$grp), FUN = sum)
  ut$grp <- factor(ut$grp, levels = c(as.character(seq_len(CAP - 1L)), sprintf("\u2265%d", CAP)))
  ut <- ut[!is.na(ut$grp), , drop = FALSE]
  pA <- ggplot(ut, aes(grp, cells)) + geom_col(fill = "#2E6E8E", width = .6) +
    geom_text(aes(label = format(cells, big.mark = ",")), vjust = -.4, size = 2.5) +
    scale_y_log10(expand = expansion(mult = c(0, .12))) +
    labs(x = "distinct UCIs in the cell", y = "cells",
         title = sprintf("%s cells with a call", format(asn, big.mark = ",")))
  pB <- ggplot(x$nocall, aes(cells, reorder(cls, cells))) + geom_col(fill = "#B4432F", width = .65) +
    geom_text(aes(label = sprintf("%s (%.1f%%)", format(cells, big.mark = ","), pct)),
              hjust = -.08, size = 2.4) +
    scale_x_continuous(expand = expansion(mult = c(0, .28))) +
    labs(x = "cells (classified by their dominant read)", y = NULL,
         title = sprintf("%s cells without a call", format(den - asn, big.mark = ",")))
  ggplot2::ggsave(op("pipcatch_nocall_fates.png"),
         patchwork::wrap_plots(pA, pB, widths = c(1, 2)),
         width = 10, height = 3.2, dpi = 300)

  ## clone matrix
  grid <- x$reference_table %>% select(set, gene, shRNA, idx) %>%
    left_join(x$clone_summary %>% select(set, gene, shRNA, clones, cells), by = c("set", "gene", "shRNA")) %>%
    mutate(clones = tidyr::replace_na(clones, 0L), cells = tidyr::replace_na(cells, 0L),
           gene = factor(gene, levels = rev(unique(x$reference_table$gene[
             order(match(x$reference_table$set, p$set.names))]))),
           label = ifelse(clones == 0, "n.d.", as.character(clones)),
           dark = cells > max(cells) * .55)
  pM <- ggplot(grid, aes(factor(idx), gene, fill = cells)) +
    geom_tile(colour = "white", linewidth = .5) +
    geom_text(aes(label = label, colour = dark), size = 2.1, show.legend = FALSE) +
    scale_fill_gradient(low = "#F7FBFF", high = "#08306B", name = "cells") +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#333333")) +
    facet_grid(set ~ ., scales = "free_y", space = "free_y") +
    labs(x = "shRNA within gene", y = NULL,
         title = sprintf("%d of %d constructs detected, %d clones",
                         sum(grid$clones > 0), nrow(grid), sum(grid$clones)),
         subtitle = "number = distinct clones (UCIs); n.d. = in library, not detected") +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6))
  ggplot2::ggsave(op("pipcatch_shRNA_clone_matrix.png"), pM, width = 5.4, height = 8.2, dpi = 300)

  ## barcode / gene distributions
  bd <- x$assignments %>% count(construct, name = "reads")
  gd <- x$assignments %>% count(gene, name = "reads")
  pbd <- ggplot(bd, aes(reorder(construct, -reads), reads)) + geom_col(fill = "#2E6E8E") +
    labs(x = NULL, y = "reads", title = sprintf("barcode distribution - %d constructs of %d, %s reads",
         nrow(bd), nrow(x$reference_table), format(nrow(x$assignments), big.mark = ","))) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 3.6))
  pgd <- ggplot(gd, aes(reorder(gene, -reads), reads)) + geom_col(fill = "#D1A13A") +
    labs(x = NULL, y = "reads", title = sprintf("gene distribution - %d targets", nrow(gd))) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 5))
  ggplot2::ggsave(op("pipcatch_QC_distributions.png"),
         patchwork::wrap_plots(pbd, pgd, ncol = 1), width = 10, height = 6.4, dpi = 300)

  .cat_plots_thresholds(x$tab, x$thresholds, x$complete_table_fin, x$params, x$output_dir,
                        file = paste0(prefix, "_QC_thresholds.png"))
  invisible(NULL)
}

.cat_plots_thresholds <- function(tab, thr, ctf, p, output_dir, file = "pipcatch_QC_thresholds.png") {
  ggplot2::theme_set(.cat_theme()); h <- thr$hist
  p1 <- ggplot(h, aes(level, count)) + geom_col(fill = "grey45", width = .7) +
    # One break per level is right at enrichment depth, where UMIxUCI spans a
    # handful of integers. On a UMI-deduplicated Parse library it spans into
    # the thousands and labelling every level makes the axis unreadable, so
    # fall back to pretty breaks once there are more than 25 of them.
    scale_y_log10() +
    scale_x_continuous(breaks = if (nrow(h) <= 25) h$level else waiver()) +
    { if (!is.na(thr$bimodal)) geom_vline(xintercept = thr$bimodal + .5, colour = "#B4432F") } +
    labs(x = "UMIxUCI (reads per cell-UCI)", y = "(cell, UCI) pairs", title = "UMIxUCI distribution")
  p2 <- ggplot(tab, aes(UMIpercentagexUCI)) + geom_histogram(binwidth = 2, fill = "grey45") +
    scale_y_log10() + geom_vline(xintercept = p$percentage, colour = "#B4432F") +
    labs(x = "UMIpercentagexUCI", y = "(cell, UCI) pairs",
         title = sprintf("percentage per UCI (minimum %.0f%%)", min(tab$UMIpercentagexUCI)))
  d2 <- ctf %>% count(UMIpercentagexUCI, UMIxUCI, True_UCIxcellID, ValidCell, name = "n")
  p3 <- ggplot(d2, aes(UMIpercentagexUCI, UMIxUCI, size = n, colour = factor(True_UCIxcellID))) +
    geom_point(alpha = .55) + scale_y_log10() + scale_size_area(max_size = 9, guide = "none") +
    geom_vline(xintercept = p$percentage, colour = "#B4432F", linewidth = .3) +
    labs(colour = "True_UCIxcellID", title = "coloured by true UCIs per cell")
  p4 <- ggplot(d2, aes(UMIpercentagexUCI, UMIxUCI, size = n, colour = ValidCell)) +
    geom_point(alpha = .55) + scale_y_log10() + scale_size_area(max_size = 9, guide = "none") +
    scale_colour_manual(values = c(`TRUE` = "#D1A13A", `FALSE` = "#8FA8BF")) +
    geom_vline(xintercept = p$percentage, colour = "#B4432F", linewidth = .3) +
    labs(title = "coloured by ValidCell")
  fig <- patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2) +
    patchwork::plot_annotation(
    title = sprintf("UMI.count %s, percentage %s, ratio %s - %s of %s cells valid",
                    p$UMI.count, p$percentage, p$ratio,
                    format(sum(ctf$ValidCell[!duplicated(ctf$cellID)]), big.mark = ","),
                    format(n_distinct(ctf$cellID), big.mark = ",")),
    theme = theme(plot.title = element_text(size = 9)))
  ggplot2::ggsave(file.path(output_dir, file), fig, width = 9.6, height = 6.6, dpi = 300)
  invisible(NULL)
}

