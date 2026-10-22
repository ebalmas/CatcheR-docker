## Part of catcheR. Cross-chemistry comparison.
##
## Everything here is derived from the objects the entry points already return,
## so the comparison costs seconds: run 10x, PIPseq and Parse once each, then
## point this at the three `.rds` files.

#' Compare barcode yield across chemistries
#'
#' Builds the one table and the one figure that put a `tenxcatch`, `pipcatch`
#' and `parsecatch` run side by side. Nothing is re-streamed: every quantity
#' comes from the objects the entry points return.
#'
#' Two choices in here are deliberate and worth knowing before reading the
#' output:
#'
#' * **The denominator is the platform's own called-cell list**, not the
#'   barcodes this package observed. The latter includes the ambient floor and
#'   overstates coverage - on a Parse enrichment library by more than tenfold.
#'   Supply `called` or the per-cell columns are left `NA` rather than being
#'   silently computed against the wrong denominator.
#' * **Per-cell depth is not on the same footing across chemistries.** A
#'   chemistry whose molecular index is long enough to deduplicate reports
#'   molecules; one whose index is too short (PIPseq's is 3 nt) reports reads.
#'   The returned `umi_dedup` column says which, and the figure says so in the
#'   panel subtitle. Do not quote the two as one number.
#'
#' @param x named list of objects, or of paths to the `.rds` each entry point
#'   writes into its `output_dir`. Names are used as figure labels.
#' @param called named character vector or list, parallel to `x`, giving the
#'   called-cell list for each run: Cell Ranger `barcodes.tsv.gz` or `.h5`,
#'   PIPseeker `barcodes.tsv.gz`, split-pipe `cell_metadata.csv.gz`, or a plain
#'   text file of IDs. `NA` or absent for a run leaves its per-cell columns `NA`.
#' @param output_dir where to write; created if absent.
#' @param file basename for the `.csv` and `.png` (no extension).
#' @param plots write the figure as well as the table.
#' @param verbose print the table.
#' @return a data frame, invisibly; also written as `<file>.csv`.
#' @export
catcheR_compare_chemistries <- function(x,
                                        called     = NULL,
                                        output_dir = "Output",
                                        file       = "chemistry_comparison",
                                        plots      = TRUE,
                                        verbose    = TRUE) {
  if (!is.list(x)) x <- as.list(x)
  if (is.null(names(x)) || any(!nzchar(names(x))))
    stop("`x` must be a NAMED list - the names label the chemistries.", call. = FALSE)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  obj <- lapply(x, function(o) if (is.character(o)) readRDS(o) else o)
  for (i in seq_along(obj))
    if (!inherits(obj[[i]], "pipcatch"))
      stop("Element ", i, " (", names(obj)[i], ") is not a catcheR run object.", call. = FALSE)

  cl <- rep(NA_character_, length(obj)); names(cl) <- names(obj)
  if (!is.null(called)) {
    called <- unlist(called)
    cl[names(called)] <- called
  }

  rows <- lapply(seq_along(obj), function(i) {
    o <- obj[[i]]
    pipe <- if (inherits(o, "tenxcatch")) "Cell Ranger"
            else if (inherits(o, "parsecatch")) "split-pipe" else "PIPseeker"
    # a per-molecule read count exists only where the index is long enough to
    # deduplicate on; its absence is how we know the depth column means reads
    dedup <- !is.null(o$assignments$reads)
    asn <- o$read_class$reads[o$read_class$class == "assigned"]
    if (!length(asn)) asn <- NA_integer_
    ccl <- if (is.na(cl[i])) NULL else
      if (inherits(o, "tenxcatch")) .tenx_called(cl[i]) else .par_called(cl[i])
    pcc <- if (is.null(ccl)) o$per_cell[0, , drop = FALSE] else
      o$per_cell[o$per_cell$cellID %in% ccl, , drop = FALSE]
    data.frame(
      chemistry            = names(obj)[i],
      pipeline             = pipe,
      umi_dedup            = dedup,
      read_pairs           = o$total_pairs,
      cellbc_reads         = o$valid_cell_reads,
      cassette_reads       = sum(o$per_sublibrary$cassette_reads),
      assigned_reads       = asn,
      events               = nrow(o$assignments),
      called_cells         = if (is.null(ccl)) NA_integer_ else length(ccl),
      called_with_cassette = if (is.null(ccl)) NA_integer_ else nrow(pcc),
      called_one_clone     = if (is.null(ccl)) NA_integer_ else
                               sum(pcc$True_UCIxcellID == 1, na.rm = TRUE),
      median_depth_called  = if (is.null(ccl) || !nrow(pcc)) NA_real_ else
                               stats::median(pcc$reads_total),
      constructs           = nrow(o$clone_summary),
      clones               = sum(o$clone_summary$clones),
      stringsAsFactors     = FALSE)
  })
  C <- do.call(rbind, rows)
  C$pct_cellbc            <- round(100 * C$cellbc_reads / C$read_pairs, 2)
  C$pct_assigned          <- round(100 * C$assigned_reads / C$read_pairs, 2)
  C$pct_called_one_clone  <- round(100 * C$called_one_clone / C$called_cells, 1)
  utils::write.csv(C, file.path(output_dir, paste0(file, ".csv")), row.names = FALSE)

  if (plots) .cat_plot_comparison(C, obj, cl, output_dir, paste0(file, ".png"))
  if (verbose) {
    message(sprintf("Comparison of %d runs -> %s",
                    nrow(C), file.path(output_dir, paste0(file, ".csv"))))
    print(C[, c("chemistry", "read_pairs", "pct_cellbc", "pct_assigned",
                "called_cells", "called_one_clone", "pct_called_one_clone",
                "median_depth_called", "constructs", "clones")])
  }
  invisible(C)
}

#' Three-panel comparison figure.
#'
#' Panel A is on a log axis on purpose: the libraries can differ by three
#' orders of magnitude in the share of reads that reach a barcode, and that
#' spread is the finding rather than a presentation problem.
.cat_plot_comparison <- function(C, obj, cl, dir, file) {
  lv  <- C$chemistry
  pal <- stats::setNames(rep(c("#2E6E8E", "#B4432F", "#D1A13A", "#8FA8BF"),
                             length.out = length(lv)), lv)
  C$chemistry <- factor(C$chemistry, levels = lv)

  st <- c("read pairs" = "read_pairs", "cell barcode parsed" = "cellbc_reads",
          "cassette found" = "cassette_reads", "barcode in reference" = "assigned_reads")
  A <- do.call(rbind, lapply(names(st), function(s) data.frame(
    chemistry = C$chemistry, stage = s, pct = 100 * C[[st[s]]] / C$read_pairs)))
  A$stage <- factor(A$stage, levels = names(st))
  end <- A[A$stage == "barcode in reference", ]
  pA <- ggplot2::ggplot(A, ggplot2::aes(stage, pct, colour = chemistry, group = chemistry)) +
    ggplot2::geom_line(linewidth = .7) + ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = pal, name = NULL) +
    ggplot2::scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                           labels = c("0.01%", "0.1%", "1%", "10%", "100%")) +
    # runs that land within a few points of each other collide at the end of
    # the line, so stagger the end labels above and below
    ggplot2::geom_text(data = end,
                       ggplot2::aes(label = sprintf("%.2f%%", pct),
                                    vjust = ifelse(pct > stats::median(end$pct), -0.9, 1.9)),
                       hjust = 1.1, size = 2.6, show.legend = FALSE) +
    ggplot2::labs(x = NULL, y = "% of that library's read pairs",
                  title = "A. Where the reads go",
                  subtitle = "each stage is a subset of the one above") +
    .cat_theme() + ggplot2::theme(legend.position = "bottom",
                                  axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

  if (all(is.na(C$called_cells))) {
    .cat_save_png(pA, dir, file, 5.0, 4.2)
    return(invisible(NULL))
  }
  B <- do.call(rbind, lapply(seq_len(nrow(C)), function(i) data.frame(
    chemistry = C$chemistry[i],
    cls = factor(c("1 clone", "other / no call"),
                 levels = c("1 clone", "other / no call")),
    cells = c(C$called_one_clone[i], C$called_cells[i] - C$called_one_clone[i]))))
  tot <- data.frame(chemistry = C$chemistry, tot = C$called_cells, pct = C$pct_called_one_clone)
  pB <- ggplot2::ggplot(B, ggplot2::aes(chemistry, cells, fill = cls)) +
    ggplot2::geom_col(width = .6) +
    ggplot2::geom_text(data = tot, ggplot2::aes(chemistry, tot,
                       label = sprintf("%s cells", format(tot, big.mark = ","))),
                       inherit.aes = FALSE, vjust = -0.5, size = 2.5, colour = "grey35") +
    ggplot2::geom_text(data = tot, ggplot2::aes(chemistry, tot * 0.5,
                       label = sprintf("%.0f%%", pct)),
                       inherit.aes = FALSE, size = 3.2, colour = "white") +
    ggplot2::scale_fill_manual(values = c("1 clone" = "#2E6E8E",
                                          "other / no call" = "#C7D3DC"), name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .14))) +
    ggplot2::labs(x = NULL, y = "cells called by the platform pipeline",
                  title = "B. Clone calls per called cell",
                  subtitle = "% = cells resolving to exactly one clone") +
    .cat_theme() + ggplot2::theme(legend.position = "bottom",
                                  axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

  D <- do.call(rbind, lapply(seq_along(obj), function(i) {
    if (is.na(cl[i])) return(NULL)
    o <- obj[[i]]
    ccl <- if (inherits(o, "tenxcatch")) .tenx_called(cl[i]) else .par_called(cl[i])
    pc <- o$per_cell[o$per_cell$cellID %in% ccl, , drop = FALSE]
    if (!nrow(pc)) return(NULL)
    data.frame(chemistry = C$chemistry[i], depth = pc$reads_total)
  }))
  D$chemistry <- factor(as.character(D$chemistry), levels = lv)
  med <- stats::aggregate(list(med = D$depth),
                          by = list(chemistry = D$chemistry), FUN = stats::median)
  readsonly <- C$chemistry[!C$umi_dedup]
  sub <- if (length(readsonly))
    sprintf("%s values are READS - index too short to deduplicate",
            paste(readsonly, collapse = ", ")) else "UMI-deduplicated molecules"
  pC <- ggplot2::ggplot(D, ggplot2::aes(chemistry, depth, fill = chemistry)) +
    ggplot2::geom_violin(scale = "width", linewidth = .2, colour = "white") +
    ggplot2::geom_boxplot(width = .13, outlier.shape = NA, linewidth = .3, fill = "white") +
    ggplot2::geom_text(data = med, ggplot2::aes(chemistry, med,
                       label = format(med, big.mark = ",")),
                       inherit.aes = FALSE, hjust = -0.9, size = 2.6, colour = "grey25") +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = NULL, y = "cassette molecules per called cell",
                  title = "C. Depth per called cell", subtitle = sub) +
    .cat_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

  fig <- patchwork::wrap_plots(pA, pB, pC, nrow = 1) +
    patchwork::plot_annotation(
      title = "Barcode yield across chemistries",
      subtitle = sprintf("catcheR %s | percentage %s, ratio %s",
                         as.character(utils::packageVersion("catcheR")),
                         obj[[1]]$params$percentage, obj[[1]]$params$ratio),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 11, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 8, colour = "grey35")))
  .cat_save_png(fig, dir, file, 4.2 * nrow(C), 4.8)
  invisible(NULL)
}
