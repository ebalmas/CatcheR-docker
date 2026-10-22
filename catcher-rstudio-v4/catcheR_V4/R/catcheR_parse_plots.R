## Part of catcheR. Analysis panels for the Parse path.
##
## These are NOT the catcheR_10Xcatch plots - those 13 are reproduced under
## their original names by the shared .cat_plots_* layer. The panels here
## answer the questions that only arise once a run is UMI-deduplicated and
## has a called-cell list to be judged against: where the reads go, where the
## barcode-rank knee falls, whether more sequencing would help, how clone
## composition varies with depth, and how the two filter parameters trade
## yield against accuracy.
##
## Every one of them is written as a .png beside the .jpg/.pdf pairs.

.par_sublib_of <- function(cellID)
  ifelse(grepl("__", cellID, fixed = TRUE), sub("^.*__", "", cellID), "all")

#' Midpoint of each segment of a stacked bar, for in-bar labels.
#'
#' `d` must already be ordered by group and then by the stacking factor.
#' Assignment is BY INDEX, not by the order `split()` happens to return its
#' groups in - getting that wrong silently puts each label on the wrong
#' segment while leaving the bars themselves correct.
.par_stack_mid <- function(d, value, grp) {
  out <- rep(NA_real_, nrow(d))
  for (i in split(seq_len(nrow(d)), grp, drop = TRUE)) {
    v <- d[[value]][i]
    out[i] <- rev(cumsum(rev(v)) - rev(v) / 2)
  }
  out
}

#' Write one analysis panel. PNG only - these are composed figures with text
#' annotations, not the jpg/pdf pairs the catcheR-named plots ship as.
.cat_save_png <- function(p, dir, file, w, h) {
  ggplot2::ggsave(file.path(dir, file), plot = p, width = w, height = h, dpi = 300)
  invisible(NULL)
}

#' Read fates, per sublibrary and pooled.
#'
#' The funnel is strictly nested: each stage is a subset of the one above, so
#' the bars are absolute read counts and never renormalised within a stage.
.par_plot_funnel <- function(x, dir, prefix = "catcheR") {
  ps <- x$per_sublibrary
  asn <- x$assignments
  asn$sublibrary <- .par_sublib_of(asn$cellID)
  # chemistries with a real UMI carry a per-molecule read count; PIPseq's 3 nt
  # index is not deduplicable, so its assignment table is one row per read and
  # the funnel counts rows instead
  asn$.reads <- if (is.null(asn$reads)) 1L else asn$reads
  ref_reads <- stats::aggregate(list(reads = asn$.reads),
                                by = list(sublibrary = asn$sublibrary), FUN = sum)
  d <- merge(ps, ref_reads, by = "sublibrary", all.x = TRUE)
  d$reads[is.na(d$reads)] <- 0

  lv <- c("read pairs", "cell barcode parsed + corrected",
          "cassette found", "shRNA barcode in reference")
  long <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) data.frame(
    sublibrary = d$sublibrary[i],
    stage = factor(lv, levels = rev(lv)),
    reads = c(d$pairs[i], d$cell_bc_reads[i], d$cassette_reads[i], d$reads[i]),
    stringsAsFactors = FALSE)))
  pooled <- stats::aggregate(list(reads = long$reads),
                             by = list(stage = long$stage), FUN = sum)
  pooled$sublibrary <- "pooled"
  long <- rbind(long, pooled[, names(long)])
  long$sublibrary <- factor(long$sublibrary, levels = c(d$sublibrary, "pooled"))
  # share of that sublibrary's own read pairs, so the label is interpretable
  tot <- stats::setNames(long$reads[long$stage == lv[1]],
                         long$sublibrary[long$stage == lv[1]])
  long$pct <- 100 * long$reads / as.numeric(tot[as.character(long$sublibrary)])
  long$lab <- sprintf("%s (%.1f%%)", format(round(long$reads), big.mark = ","), long$pct)

  p <- ggplot2::ggplot(long, ggplot2::aes(x = stage, y = reads)) +
    ggplot2::geom_col(fill = "#2E6E8E", width = .68) +
    ggplot2::geom_text(ggplot2::aes(label = lab), hjust = -0.05, size = 2.5) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ sublibrary, ncol = 1, scales = "free_x") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .38))) +
    ggplot2::labs(x = NULL, y = "read pairs",
                  title = "Where the reads go",
                  subtitle = "each stage is a subset of the one above") +
    .cat_theme()
  .cat_save_png(p, dir, paste0(prefix, "_detection_funnel.png"), 7.2,
                2.1 * nlevels(long$sublibrary) + 0.8)
  invisible(long)
}

#' Barcode-rank ("knee") curve.
#'
#' Barcodes sorted by cassette UMIs. The knee is where real cells stop and
#' ambient barcodes begin; comparing it to the cell count the chemistry's own
#' pipeline called is the cheapest check that the extraction is sane.
.par_plot_rank <- function(x, dir, min.cell.UMI = 20, n_called = NULL,
                           prefix = "catcheR") {
  pc <- x$per_cell
  u <- sort(pc$reads_total, decreasing = TRUE)
  d <- data.frame(rank = seq_along(u), umis = u)
  # thin the flat tail; plotting 50k identical points is slow and adds nothing
  keep <- unique(c(seq_len(min(500, nrow(d))),
                   round(exp(seq(log(min(500, nrow(d))), log(nrow(d)), length.out = 2000)))))
  d <- d[keep[keep <= nrow(d)], ]
  above <- sum(u >= min.cell.UMI)

  p <- ggplot2::ggplot(d, ggplot2::aes(rank, umis)) +
    ggplot2::geom_line(colour = "#2E6E8E", linewidth = .6) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::geom_hline(yintercept = min.cell.UMI, colour = "#B4432F", linetype = "dashed",
                        linewidth = .4) +
    ggplot2::annotate("text", x = 1.2, y = min.cell.UMI * 1.35, hjust = 0, size = 2.6,
                      colour = "#B4432F",
                      label = sprintf("min.cell.UMI = %s  \u2192  %s barcodes",
                                      min.cell.UMI, format(above, big.mark = ","))) +
    ggplot2::labs(x = "barcode rank", y = "cassette UMIs",
                  title = "Barcode-rank curve") +
    .cat_theme()
  if (!is.null(n_called) && is.finite(n_called))
    p <- p + ggplot2::geom_vline(xintercept = n_called, colour = "#D1A13A", linewidth = .4) +
      # the line can sit at the far right of the rank axis (10x, where the
      # whitelist IS the called-cell list), so flip the label inside the panel
      ggplot2::annotate("text", x = n_called, y = max(d$umis),
                        hjust = if (n_called > 0.6 * max(d$rank)) 1.05 else -0.05,
                        size = 2.6, colour = "#D1A13A",
                        label = sprintf("%s cells called", format(n_called, big.mark = ",")))
  .cat_save_png(p, dir, paste0(prefix, "_barcode_rank.png"), 5.6, 3.6)
  invisible(d)
}

#' Exact rarefaction: would more sequencing help?
#'
#' Subsamples the READS supporting each molecule, which is why
#' `catcheR_parsecatch()` keeps a per-molecule read count. The detection
#' probability is `1 - (1 - t)^reads`, which is exact and returns the observed
#' counts at t = 1 - unlike `1 - exp(-t * reads)`, which asymptotes below the
#' data because it cannot rediscover molecules that were never seen.
#' Extrapolation beyond the data is deliberately not offered.
.par_plot_saturation <- function(x, dir, min.cell.UMI = 20, percentage = 15,
                                 fracs = seq(0.1, 1, by = 0.1), seed = 1,
                                 prefix = "catcheR") {
  a <- x$assignments
  if (is.null(a$reads)) {
    warning("No per-molecule read counts on this object; skipping the saturation curve. ",
            "Re-run catcheR_parsecatch() to produce them.", call. = FALSE)
    return(invisible(NULL))
  }
  set.seed(seed)
  cid <- as.integer(factor(a$cellID)); uid <- as.integer(factor(paste(a$barcode, a$UCI)))
  k <- a$reads; ncell <- max(cid)
  rows <- lapply(fracs, function(t) {
    det <- if (t >= 1) rep(TRUE, length(k)) else stats::runif(length(k)) < (1 - (1 - t)^k)
    um <- tabulate(cid[det], nbins = ncell)
    deep <- which(um >= min.cell.UMI)
    g <- stats::aggregate(list(n = rep(1, sum(det))),
                          by = list(cid = cid[det], uid = uid[det]), FUN = sum)
    g$share <- 100 * g$n / um[g$cid]
    dom <- unique(g$cid[g$share >= percentage])
    data.frame(frac = t, reads_M = round(t * sum(k) / 1e6, 2),
               molecules = sum(det),
               cells_deep = length(deep),
               with_dominant = length(intersect(dom, deep)),
               constructs = length(unique(a$shRNA[det])))
  })
  R <- do.call(rbind, rows)
  utils::write.csv(R, file.path(dir, paste0(prefix, "_saturation.csv")), row.names = FALSE)

  long <- rbind(
    data.frame(reads_M = R$reads_M, n = R$cells_deep,
               series = sprintf("barcodes \u2265%s UMIs", min.cell.UMI)),
    data.frame(reads_M = R$reads_M, n = R$with_dominant,
               series = "with a dominant clone"))
  p1 <- ggplot2::ggplot(long, ggplot2::aes(reads_M, n, colour = series)) +
    ggplot2::geom_line(linewidth = .6) + ggplot2::geom_point(size = 1.1) +
    ggplot2::scale_colour_manual(values = c("#8FA8BF", "#2E6E8E"), name = NULL) +
    ggplot2::labs(x = "read pairs used (millions)", y = "barcodes",
                  title = "Rarefaction") +
    .cat_theme() + ggplot2::theme(legend.position = "bottom")

  mid <- (R$reads_M[-1] + R$reads_M[-nrow(R)]) / 2
  d2 <- data.frame(reads_M = mid,
                   slope = diff(R$with_dominant) / diff(R$reads_M))
  p2 <- ggplot2::ggplot(d2, ggplot2::aes(reads_M, slope)) +
    ggplot2::geom_line(colour = "#2E6E8E", linewidth = .6) +
    ggplot2::geom_point(colour = "#2E6E8E", size = 1.1) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(x = "read pairs used (millions)", y = "new calls per million reads",
                  title = sprintf("Marginal return (%.0f \u2192 %.0f)",
                                  d2$slope[1], d2$slope[nrow(d2)])) +
    .cat_theme()
  .cat_save_png(patchwork::wrap_plots(p1, p2, nrow = 1), dir,
                paste0(prefix, "_saturation.png"), 8.4, 3.4)
  invisible(R)
}

#' Clone composition per cell, at a ladder of per-cell UMI cutoffs.
#'
#' Absolute counts, never renormalised within a cutoff, with the total number
#' of barcodes at that cutoff printed above each stack - otherwise a rising
#' "single clone" percentage hides the cells being thrown away to achieve it.
.par_plot_uci_per_cell <- function(x, dir, cutoffs = c(1, 5, 10, 20, 50),
                                   called = NULL, prefix = "catcheR",
                                   file = NULL) {
  if (is.null(file)) file <- paste0(prefix, "_UCI_per_cell.png")
  pc <- x$per_cell
  called <- .par_called(called)
  unit <- "barcodes"
  scope <- "every barcode observed, ambient floor included"
  if (!is.null(called)) {
    pc <- pc[pc$cellID %in% called, , drop = FALSE]
    unit <- "cells"
    scope <- "restricted to the cells called by the chemistry's pipeline"
    if (!nrow(pc)) {
      warning("No called cell carries a cassette UMI; skipping ", file, call. = FALSE)
      return(invisible(NULL))
    }
  }
  pc$sublibrary <- .par_sublib_of(pc$cellID)
  cls <- function(n) ifelse(n == 0, "no dominant UCI",
                     ifelse(n == 1, "1 clone", ifelse(n == 2, "2 clones", "\u22653 clones")))
  parts <- lapply(cutoffs, function(k) {
    s <- pc[pc$reads_total >= k, , drop = FALSE]
    if (!nrow(s)) return(NULL)
    one <- data.frame(cutoff = k, sublibrary = s$sublibrary, cl = cls(s$True_UCIxcellID))
    if (length(unique(pc$sublibrary)) > 1L)
      one <- rbind(one, data.frame(cutoff = k, sublibrary = "pooled",
                                   cl = cls(s$True_UCIxcellID)))
    one
  })
  d <- do.call(rbind, parts)
  d <- stats::aggregate(list(cells = rep(1, nrow(d))),
                        by = list(cutoff = d$cutoff, sublibrary = d$sublibrary, cl = d$cl),
                        FUN = sum)
  d$cl <- factor(d$cl, levels = c("1 clone", "2 clones", "\u22653 clones", "no dominant UCI"))
  # sublibraries in run order, pooled last - alphabetical would put "pooled" first
  d$sublibrary <- factor(d$sublibrary,
                         levels = c(setdiff(unique(pc$sublibrary), "pooled"), "pooled"))
  tots <- stats::aggregate(list(tot = d$cells),
                           by = list(cutoff = d$cutoff, sublibrary = d$sublibrary), FUN = sum)
  tots$sublibrary <- factor(tots$sublibrary, levels = levels(d$sublibrary))
  # percentage of the bar, written inside the segment, with the absolute total
  # still printed above it: the percentage denominator shifts from bar to bar
  # here (that is the point of the panel), so the count has to stay visible.
  d <- merge(d, tots, by = c("cutoff", "sublibrary"))
  d$pct <- 100 * d$cells / d$tot
  d <- d[order(d$sublibrary, d$cutoff, d$cl), ]
  d$ymid <- .par_stack_mid(d, "cells", interaction(d$sublibrary, d$cutoff, drop = TRUE))
  # a label needs both a meaningful share AND enough pixels to sit in; the
  # facets share one y scale, so height is judged against the tallest bar
  d$lab_ok <- d$pct >= 8 & d$cells >= 0.03 * max(tots$tot)

  p <- ggplot2::ggplot(d, ggplot2::aes(factor(cutoff), cells, fill = cl)) +
    ggplot2::geom_col(width = .7) +
    ggplot2::geom_text(data = d[d$lab_ok, ],
                       ggplot2::aes(y = ymid, label = sprintf("%.0f%%", pct)),
                       size = 2.2, colour = "white") +
    ggplot2::geom_text(data = tots, ggplot2::aes(factor(cutoff), tot, label =
                         format(tot, big.mark = ",")), inherit.aes = FALSE,
                       vjust = -0.4, size = 2.3) +
    ggplot2::facet_wrap(~ sublibrary, nrow = 1) +
    ggplot2::scale_fill_manual(values = c("1 clone" = "#2E6E8E", "2 clones" = "#8FA8BF",
                                          "\u22653 clones" = "#C7D3DC",
                                          "no dominant UCI" = "#D1A13A"), name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .12))) +
    ggplot2::labs(x = "minimum cassette UMIs in the cell", y = unit,
                  title = "Clone composition against the per-cell UMI floor",
                  subtitle = sprintf("absolute counts, %% of the bar inside each segment; %s", scope)) +
    .cat_theme() + ggplot2::theme(legend.position = "bottom")
  .cat_save_png(p, dir, file, 9.0, 3.8)
  invisible(d)
}

#' Resolve the called-cell list supplied by the user.
#'
#' Accepts the chemistry's own cell table (`cell_metadata.csv[.gz]`, whose
#' `bc_wells` column holds the same cell ID this package emits), a plain file
#' of IDs, or a character vector. Returns NULL when nothing was supplied, which
#' every caller treats as "no called-cell list available".
.par_called <- function(called) {
  if (is.null(called)) return(NULL)
  if (length(called) == 1L && file.exists(called)) {
    if (grepl("\\.csv(\\.gz)?$", called)) {
      d <- utils::read.csv(called, stringsAsFactors = FALSE)
      col <- intersect(c("bc_wells", "cellID", "cell_barcode", "barcode"), names(d))
      if (!length(col))
        stop("Could not find a cell-ID column in ", basename(called),
             ". Expected one of bc_wells / cellID / cell_barcode / barcode; got: ",
             paste(names(d), collapse = ", "), call. = FALSE)
      called <- d[[col[1]]]
    } else {
      called <- readLines(called)
    }
  }
  unique(as.character(called))
}

#' Clone composition of the cells the chemistry's own pipeline called.
#'
#' The other composition panels are over every barcode this package observed,
#' which includes the ambient floor. This one uses the called-cell list as a
#' FIXED denominator, so a percentage axis is honest here: it does not shift
#' from bar to bar. The absolute cell count is printed above each bar anyway.
#'
#' A called cell absent from `per_cell` carried no cassette UMI at all, and is
#' shown as its own category rather than dropped - it is a real outcome of the
#' enrichment, not a missing value.
.par_plot_called_composition <- function(x, called, dir, file = NULL,
                                        prefix = "catcheR") {
  called <- .par_called(called)
  if (is.null(called)) return(invisible(NULL))
  pc <- x$per_cell
  LV <- c("1 clone", "2 clones", "\u22653 clones", "no dominant UCI", "no cassette UMI")
  n <- match(called, pc$cellID)
  cls <- ifelse(is.na(n), "no cassette UMI", {
    k <- pc$True_UCIxcellID[n]
    ifelse(is.na(k) | k == 0, "no dominant UCI",
    ifelse(k == 1, "1 clone", ifelse(k == 2, "2 clones", "\u22653 clones")))
  })
  sub <- .par_sublib_of(called)
  # with one library the pooled panel would be an identical twin of it
  d <- data.frame(sublibrary = sub, cls = cls, stringsAsFactors = FALSE)
  if (length(unique(sub)) > 1L)
    d <- rbind(d, data.frame(sublibrary = "pooled", cls = cls, stringsAsFactors = FALSE))
  d <- stats::aggregate(list(cells = rep(1L, nrow(d))),
                        by = list(sublibrary = d$sublibrary, cls = d$cls), FUN = sum)
  tot <- stats::aggregate(list(tot = d$cells), by = list(sublibrary = d$sublibrary), FUN = sum)
  d <- merge(d, tot, by = "sublibrary")
  d$pct <- 100 * d$cells / d$tot
  d$cls <- factor(d$cls, levels = LV)
  lev <- c(setdiff(sort(unique(sub)), "pooled"), "pooled")
  d$sublibrary   <- factor(d$sublibrary, levels = lev)
  tot$sublibrary <- factor(tot$sublibrary, levels = lev)
  d <- d[order(d$sublibrary, d$cls), ]
  # label position: centre of each stacked segment, counting from the top of
  # the stack because geom_col stacks the first factor level last
  d$ymid <- .par_stack_mid(d, "pct", d$sublibrary)
  # Headline over ALL called cells. It used to read the "pooled" row, which a
  # single-library chemistry no longer has, and silently printed NA.
  one <- d$cls == "1 clone"
  hd <- if (any(d$sublibrary == "pooled")) d$sublibrary == "pooled" else rep(TRUE, nrow(d))
  ttl <- if (any(one & hd))
    sprintf("%.0f%% of the called cells get exactly one clone",
            sum(d$pct[one & hd] * d$tot[one & hd]) / sum(d$tot[one & hd]))
  else "Composition of the called cells"

  p <- ggplot2::ggplot(d, ggplot2::aes(sublibrary, pct, fill = cls)) +
    ggplot2::geom_col(width = .62) +
    ggplot2::geom_text(data = d[d$pct >= 4, ],
                       ggplot2::aes(y = ymid, label = sprintf("%.0f%%", pct)),
                       size = 2.6, colour = "white") +
    ggplot2::geom_text(data = tot, ggplot2::aes(sublibrary, 103,
                       label = sprintf("%s cells", format(tot, big.mark = ","))),
                       inherit.aes = FALSE, size = 2.5, colour = "grey35") +
    ggplot2::scale_fill_manual(values = c("1 clone" = "#2E6E8E", "2 clones" = "#D1A13A",
                                          "\u22653 clones" = "#B4432F",
                                          "no dominant UCI" = "#C7D3DC",
                                          "no cassette UMI" = "#EFEFEF"), name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 20),
                                expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(x = NULL, y = "% of the cells called by the chemistry's pipeline",
                  title = ttl,
                  subtitle = "fixed denominator per panel; absolute cell count above each bar") +
    .cat_theme() + ggplot2::theme(legend.position = "bottom")
  if (is.null(file)) file <- paste0(prefix, "_called_composition.png")
  .cat_save_png(p, dir, file, 1.5 + 1.3 * nlevels(d$sublibrary), 4.0)
  utils::write.csv(d[, c("sublibrary", "cls", "cells", "tot", "pct")],
                   file.path(dir, sub("\\.png$", ".csv", file)), row.names = FALSE)
  invisible(d)
}

#' Yield against expected error, as a function of the per-cell UMI floor.
#'
#' The one panel that makes the threshold decision directly: the left axis is
#' how many cells still receive a clone call, the right axis is the share of
#' those calls expected to be WRONG, taken from [catcheR_call_accuracy()] at
#' the matching depth. A floor is worth raising only while the error is still
#' falling; past that point it costs cells and buys nothing.
.par_plot_yield_error <- function(x, acc, dir, min.cell.UMI = 20, called = NULL,
                                  prefix = "catcheR", file = NULL) {
  if (is.null(file)) file <- paste0(prefix, "_yield_vs_error.png")
  if (is.null(acc) || !nrow(acc)) return(invisible(NULL))
  pc <- x$per_cell
  called <- .par_called(called)
  scope <- if (is.null(called)) "all observed barcodes" else {
    pc <- pc[pc$cellID %in% called, , drop = FALSE]
    "cells called by the chemistry's pipeline"
  }
  floors <- sort(unique(c(0, acc$UMIs_in_cell)))
  yield <- vapply(floors, function(k)
    sum(pc$reads_total >= k & pc$True_UCIxcellID == 1), integer(1))
  err <- acc$error_among_calls[match(floors, acc$UMIs_in_cell)]
  d <- data.frame(floor = floors, yield = yield, err = err)

  # second axis: map the error range onto the yield range
  yr <- range(d$yield); er <- range(d$err, na.rm = TRUE)
  if (diff(er) <= 0) er <- er + c(0, 1)
  to_y <- function(e) yr[1] + (e - er[1]) / diff(er) * diff(yr)
  to_e <- function(y) er[1] + (y - yr[1]) / diff(yr) * diff(er)

  lab <- d[d$floor %in% c(0, min.cell.UMI, max(d$floor)), ]
  p <- ggplot2::ggplot(d, ggplot2::aes(floor, yield)) +
    ggplot2::geom_vline(xintercept = min.cell.UMI, colour = "grey55",
                        linetype = "dotted", linewidth = .4) +
    ggplot2::geom_line(colour = "#2E6E8E", linewidth = .7) +
    ggplot2::geom_point(colour = "#2E6E8E", size = 1.6) +
    ggplot2::geom_line(ggplot2::aes(y = to_y(err)), colour = "#B4432F",
                       linewidth = .6, linetype = "dashed", na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(y = to_y(err)), colour = "#B4432F",
                       size = 1.4, shape = 15, na.rm = TRUE) +
    ggplot2::geom_text(data = lab, ggplot2::aes(label = format(yield, big.mark = ",")),
                       vjust = -1, size = 2.6, colour = "#2E6E8E") +
    ggplot2::annotate("text", x = min.cell.UMI, y = yr[2], hjust = -0.1, vjust = 1.4,
                      size = 2.5, colour = "grey40",
                      label = sprintf("%s UMIs", min.cell.UMI)) +
    ggplot2::scale_y_continuous(
      name = "cells with a clone call",
      expand = ggplot2::expansion(mult = c(0.08, 0.12)),
      sec.axis = ggplot2::sec_axis(~ to_e(.), name = "expected wrong calls (%)")) +
    ggplot2::labs(x = "minimum cassette UMIs required in the cell",
                  title = "Where the error stops falling and the yield keeps dropping",
                  subtitle = scope) +
    .cat_theme() +
    ggplot2::theme(axis.title.y.left  = ggplot2::element_text(colour = "#2E6E8E"),
                   axis.text.y.left   = ggplot2::element_text(colour = "#2E6E8E"),
                   axis.title.y.right = ggplot2::element_text(colour = "#B4432F"),
                   axis.text.y.right  = ggplot2::element_text(colour = "#B4432F"))
  .cat_save_png(p, dir, file, 6.2, 4.0)
  utils::write.csv(d, file.path(dir, sub("\\.png$", ".csv", file)), row.names = FALSE)
  invisible(d)
}

#' Yield and accuracy against the two filter parameters.
.par_plot_diagnostics <- function(acc = NULL, sweep = NULL, dir,
                                  percentage = 15, min.cell.UMI = 20,
                                  prefix = "catcheR") {
  if (!is.null(acc)) {
    long <- rbind(
      data.frame(UMIs = acc$UMIs_in_cell, pct = acc$correct, what = "correct call"),
      data.frame(UMIs = acc$UMIs_in_cell, pct = acc$wrong,   what = "WRONG call"),
      data.frame(UMIs = acc$UMIs_in_cell, pct = acc$no_call, what = "no call"))
    p <- ggplot2::ggplot(long, ggplot2::aes(UMIs, pct, colour = what)) +
      ggplot2::geom_line(linewidth = .6) + ggplot2::geom_point(size = 1.1) +
      ggplot2::geom_vline(xintercept = min.cell.UMI, colour = "#B4432F",
                          linetype = "dashed", linewidth = .4) +
      ggplot2::scale_x_log10() +
      ggplot2::scale_colour_manual(values = c("correct call" = "#2E6E8E",
                                              "WRONG call" = "#B4432F",
                                              "no call" = "#C7D3DC"), name = NULL) +
      ggplot2::labs(x = "cassette UMIs in the cell (subsampled)", y = "% of reference cells",
                    title = "Call accuracy against depth",
                    subtitle = "a wrong call is worse than no call: it adds a cross-clone assignment") +
      .cat_theme() + ggplot2::theme(legend.position = "bottom")
    .cat_save_png(p, dir, paste0(prefix, "_call_accuracy.png"), 6.0, 3.8)
  }
  if (!is.null(sweep)) {
    s <- sweep[sweep$min.cell.UMI == min.cell.UMI, , drop = FALSE]
    # a floor tuned for a UMI-deduplicated library can exclude every cell of a
    # read-counted one, which draws an all-zero chart instead of saying so
    floor_used <- min.cell.UMI
    if (!nrow(s) || all(s$cells_kept == 0)) {
      alt <- sweep[sweep$cells_kept > 0, , drop = FALSE]
      if (!nrow(alt)) return(invisible(NULL))
      # Pick the floor NEAREST the one asked for, ties going to the lower - not
      # the highest available, which is the most restrictive and was leaving a
      # panel drawn on a handful of cells. On a read-counted library the
      # requested floor is often absent from the grid entirely (2 against a grid
      # of 0/5/10/20/50), and 0 is both the nearest and the informative choice.
      cand <- sort(unique(alt$min.cell.UMI))
      floor_used <- cand[which.min(abs(cand - min.cell.UMI))]
      s <- sweep[sweep$min.cell.UMI == floor_used, , drop = FALSE]
    }
    long <- rbind(
      data.frame(percentage = s$percentage, UMI.count = s$UMI.count,
                 cells = s$one_clone, what = "1 clone"),
      data.frame(percentage = s$percentage, UMI.count = s$UMI.count,
                 cells = s$ambiguous, what = "ambiguous"),
      data.frame(percentage = s$percentage, UMI.count = s$UMI.count,
                 cells = s$none_pass, what = "no UCI passes"))
    p <- ggplot2::ggplot(long, ggplot2::aes(percentage, cells, colour = what)) +
      ggplot2::geom_line(linewidth = .6) + ggplot2::geom_point(size = 1) +
      ggplot2::facet_wrap(~ paste("UMI.count =", UMI.count), nrow = 1) +
      ggplot2::scale_colour_manual(values = c("1 clone" = "#2E6E8E",
                                              "ambiguous" = "#D1A13A",
                                              "no UCI passes" = "#C7D3DC"), name = NULL) +
      ggplot2::labs(x = "percentage", y = "cells",
                    # This asserted "lowering the percentage mostly converts
                    # 'no call' into 'ambiguous'", which is what the Parse
                    # library shows and is plainly false on a read-counted one,
                    # where the lines are flat. State the measured effect at the
                    # most permissive UMI.count instead of a fixed conclusion.
                    title = local({
                      u  <- s[s$UMI.count == min(s$UMI.count), , drop = FALSE]
                      u  <- u[order(u$percentage), ]
                      d1 <- u$one_clone[1] - u$one_clone[nrow(u)]
                      da <- u$ambiguous[1] - u$ambiguous[nrow(u)]
                      if (max(abs(c(d1, da))) < 0.02 * max(u$cells_kept, 1))
                        sprintf("percentage %s to %s barely moves this library (%+d one-clone cells)",
                                max(u$percentage), min(u$percentage), d1)
                      else
                        sprintf("percentage %s to %s: %+d cells with one clone, %+d ambiguous",
                                max(u$percentage), min(u$percentage), d1, da)
                    }),
                    subtitle = sprintf("cells with \u2265%s cassette UMIs%s", floor_used,
                                       if (floor_used != min.cell.UMI)
                                         sprintf(" - the requested floor of %s keeps no cell on this library",
                                                 min.cell.UMI)
                                       else "")) +
      .cat_theme() + ggplot2::theme(legend.position = "bottom")
    .cat_save_png(p, dir, paste0(prefix, "_threshold_sweep.png"), 9.0, 3.4)
  }
  invisible(NULL)
}

#' Run every analysis panel that the supplied object can support.
.par_plots_analysis <- function(x, dir = x$output_dir, min.cell.UMI = NULL,
                                percentage = NULL, n_called = NULL, called = NULL,
                                prefix = NULL, verbose = TRUE) {
  if (is.null(prefix))
    prefix <- if (inherits(x, "tenxcatch")) "tenxcatch"
              else if (inherits(x, "parsecatch")) "parsecatch" else "pipcatch"
  called <- .par_called(called)
  if (!is.null(called) && is.null(n_called)) n_called <- length(called)
  if (is.null(min.cell.UMI))
    # PIPseq has no deduplicable UMI, so its per-cell counts are reads and top
    # out around 10; a floor of 20 would exclude every cell. Objects that do
    # carry the parameter (Parse, 10x) use their own.
    min.cell.UMI <- if (!is.null(x$params$min.cell.UMI)) x$params$min.cell.UMI else 2
  if (is.null(percentage)) percentage <- x$params$percentage
  .par_plot_funnel(x, dir, prefix)
  .par_plot_rank(x, dir, min.cell.UMI, n_called, prefix)
  .par_plot_uci_per_cell(x, dir, prefix = prefix)
  # the same panel restricted to the called cells, plus their composition on a
  # fixed denominator - only possible when a called-cell list was supplied
  .par_plot_uci_per_cell(x, dir, called = called, prefix = prefix,
                         file = paste0(prefix, "_UCI_per_cell_called.png"))
  .par_plot_called_composition(x, called, dir, prefix = prefix)
  .par_plot_saturation(x, dir, min.cell.UMI, percentage, prefix = prefix)
  # A diagnostic can legitimately be impossible on a given library - the
  # accuracy test needs cells deep enough to serve as ground truth, the
  # saturation curve needs a deduplicable index. Say which panel is missing and
  # why, rather than silently writing one file fewer than the docs promise.
  acc <- try(catcheR_call_accuracy(x, output_dir = dir, verbose = FALSE), silent = TRUE)
  if (inherits(acc, "try-error") && verbose)
    message("  skipping call accuracy and yield-vs-error: ",
            conditionMessage(attr(acc, "condition")))
  swp <- try(catcheR_threshold_sweep(x, output_dir = dir, verbose = FALSE), silent = TRUE)
  if (inherits(swp, "try-error") && verbose)
    message("  skipping the threshold sweep: ",
            conditionMessage(attr(swp, "condition")))
  .par_plot_yield_error(x, if (inherits(acc, "try-error")) NULL else acc, dir,
                        min.cell.UMI, called, prefix)
  .par_plot_diagnostics(if (inherits(acc, "try-error")) NULL else acc,
                        if (inherits(swp, "try-error")) NULL else swp,
                        dir, percentage, min.cell.UMI, prefix)
  if (verbose) message("  analysis panels written to ", dir)
  invisible(NULL)
}
