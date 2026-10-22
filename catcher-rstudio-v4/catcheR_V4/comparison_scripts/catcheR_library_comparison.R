## Comparison of the pooled BC libraries and of the sequencing sublibraries,
## for one catcheR run. Standalone: `source()` this file, no package changes.
##
##   res <- readRDS("analysis_EB_PIP/Output/pipcatch_AB019b.rds")
##   L <- compare_bc_libraries(res, called = ".../barcodes.tsv.gz", output_dir = "out")
##   S <- compare_sublibraries(res, called = ".../cell_metadata.csv.gz", output_dir = "out")
##
## Works on a tenxcatch, pipcatch or parsecatch object. Needs ggplot2 and
## patchwork; hdf5r only if you pass a Cell Ranger .h5 as the cell list.
## Nothing here calls catcheR internals, so it keeps working if the package
## changes underneath it.

suppressMessages({ library(ggplot2); library(patchwork) })

# --- small helpers, local copies so the file stands on its own ---------------

cat_theme <- function() theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 9, hjust = 0))

cat_save_png <- function(p, dir, file, w, h) {
  ggsave(file.path(dir, file), plot = p, width = w, height = h, dpi = 300)
  invisible(NULL)
}

# the sublibrary tag Parse appends to the cell ID; "" when there is none
cat_sublib_of <- function(cellID) sub("^.*__", "", cellID)

#' Read a called-cell list: Cell Ranger barcodes.tsv.gz or .h5, PIPseeker
#' barcodes.tsv.gz, split-pipe cell_metadata.csv.gz, or a plain text file of
#' IDs. The `-1` GEM-well suffix is stripped, as catcheR does.
cat_read_called <- function(called) {
  if (is.null(called)) return(NULL)
  if (length(called) == 1L && file.exists(called)) {
    if (grepl("\\.h5$", called)) {
      if (!requireNamespace("hdf5r", quietly = TRUE))
        stop("Reading a .h5 needs the hdf5r package; pass barcodes.tsv.gz instead.",
             call. = FALSE)
      h <- hdf5r::H5File$new(called, mode = "r"); on.exit(h$close_all())
      nm <- if (h$exists("matrix/barcodes")) "matrix/barcodes" else "GRCh38/barcodes"
      called <- h[[nm]]$read()
    } else if (grepl("\\.csv(\\.gz)?$", called)) {
      d <- utils::read.csv(called, stringsAsFactors = FALSE)
      col <- intersect(c("bc_wells", "cellID", "cell_barcode", "barcode"), names(d))
      if (!length(col))
        stop("No cell-ID column in ", basename(called), "; got: ",
             paste(names(d), collapse = ", "), call. = FALSE)
      called <- d[[col[1]]]
    } else {
      called <- readLines(if (grepl("\\.gz$", called)) gzfile(called) else called)
    }
  }
  unique(sub("-[0-9]+$", "", as.character(called)))
}

## Part of catcheR. Comparison of the pooled shRNA libraries inside one run.
##
## When two or more barcode-to-gene CSVs are pooled before loading, every
## construct is tagged with the set it came from. This asks how each library
## behaved in that pool: how many of its constructs came back, how evenly, and
## how many cells it accounts for.

# Compare the pooled shRNA libraries within a run
#
# Works on any run object (`tenxcatch`, `pipcatch`, `parsecatch`). The
# denominator per library is the number of constructs that library *contains*,
# taken from the run's own reference table, so a library that contributed 60
# constructs is never scored against 100.
#
# Two things this is for:
#
# * **Dropout** — constructs present in the reference and never detected. A
#   library with systematically more dropout than its partner was either
#   under-represented in the pool or lost in cloning, and that is visible here
#   before any biology is interpreted.
# * **Evenness** — how skewed the representation is within each library.
#   Reported as the coefficient of variation of cells per construct, the share
#   of cells held by the top decile of constructs, and the Gini coefficient.
#   A library can detect every construct and still be unusable if 10% of them
#   hold most of the cells.
#
# Cells carrying calls from more than one library are counted once, under the
# library of their dominant clone, and reported separately as `cross_set_cells`
# (an attribute of the returned table) because they are candidate doublets
# rather than evidence about either library.
#
# @param x a run object or the path to its `.rds`.
# @param called optional called-cell list; when given, the per-cell columns are
#   restricted to cells the platform's own pipeline called. Strongly
#   recommended - on an enrichment library the observed barcodes include the
#   ambient floor.
# @param output_dir where to write; created if absent.
# @param file basename for the `.csv` and `.png`.
# @param plots write the figure as well as the table.
# @param verbose print the table.
# @return a data frame, invisibly, with `cross_set_cells` and `per_construct`
#   attached as attributes; also written as `<file>.csv`.
compare_bc_libraries <- function(x,
                                      called     = NULL,
                                      output_dir = NULL,
                                      file       = NULL,
                                      plots      = TRUE,
                                      verbose    = TRUE) {
  if (is.character(x)) x <- readRDS(x)
  stopifnot(is.list(x), !is.null(x$per_cell), !is.null(x$clone_summary))
  if (is.null(output_dir)) output_dir <- x$output_dir
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  pfx <- if (inherits(x, "tenxcatch")) "tenxcatch"
         else if (inherits(x, "parsecatch")) "parsecatch" else "pipcatch"
  if (is.null(file)) file <- paste0(pfx, "_library_comparison")

  ref <- x$reference_table
  if (is.null(ref) || !"set" %in% names(ref))
    stop("This run has no per-library tagging: its reference table carries no ",
         "`set` column. Pass several barcode CSVs with `set.names` to get it.",
         call. = FALSE)
  sets <- unique(ref$set)
  if (length(sets) < 2L)
    warning("Only one library (", sets, ") in this run; the comparison has ",
            "nothing to compare it against.", call. = FALSE)

  pc <- x$per_cell
  ccl <- cat_read_called(called)
  if (!is.null(ccl)) pc <- pc[pc$cellID %in% ccl, , drop = FALSE]

  # cells per construct, zero-filled over EVERY construct in the reference, so
  # dropouts are rows of zero rather than missing rows
  # Key on (shRNA, set), never on shRNA alone: the control constructs - SCR,
  # LUCIFERASE, B2M, ARHGDIA - carry the SAME shRNA name in both libraries, so
  # merging on the name credits one detection to both and inflates every
  # per-library count. `construct` is the disambiguated name; `top_set` says
  # which library the cell's dominant call actually came from.
  cpc <- stats::aggregate(list(cells = pc$cellID),
                          by = list(shRNA = pc$top_shRNA, set = pc$top_set),
                          FUN = function(z) length(unique(z)))
  per_construct <- merge(ref[, c("construct", "shRNA", "gene", "set")], cpc,
                         by = c("shRNA", "set"), all.x = TRUE)
  per_construct$cells[is.na(per_construct$cells)] <- 0L

  gini <- function(v) {
    v <- sort(v); n <- length(v)
    if (!n || sum(v) == 0) return(NA_real_)
    sum((2 * seq_len(n) - n - 1) * v) / (n * sum(v))
  }
  rows <- lapply(sets, function(s) {
    pcs <- per_construct[per_construct$set == s, , drop = FALSE]
    cs  <- x$clone_summary[x$clone_summary$set == s, , drop = FALSE]
    asn <- x$assignments[x$assignments$set == s, , drop = FALSE]
    v   <- pcs$cells
    top <- sort(v, decreasing = TRUE)[seq_len(max(1L, ceiling(length(v) / 10)))]
    data.frame(
      library              = s,
      constructs_expected  = nrow(pcs),
      constructs_detected  = sum(v > 0),
      dropout              = sum(v == 0),
      pct_detected         = round(100 * sum(v > 0) / nrow(pcs), 1),
      genes_expected       = length(unique(pcs$gene)),
      genes_detected       = length(unique(pcs$gene[v > 0])),
      cells                = sum(v),
      cells_per_construct  = round(sum(v) / nrow(pcs), 1),
      median_cells         = stats::median(v),
      clones               = sum(cs$clones),
      events               = nrow(asn),
      cv_cells             = round(stats::sd(v) / mean(v), 2),
      top_decile_share     = round(100 * sum(top) / max(1, sum(v)), 1),
      gini_cells           = round(gini(v), 2),
      stringsAsFactors     = FALSE)
  })
  L <- do.call(rbind, rows)
  cross <- sum(pc$n_sets > 1, na.rm = TRUE)
  attr(L, "cross_set_cells") <- cross
  attr(L, "per_construct")   <- per_construct
  utils::write.csv(L, file.path(output_dir, paste0(file, ".csv")), row.names = FALSE)
  utils::write.csv(per_construct,
                   file.path(output_dir, paste0(file, "_per_construct.csv")), row.names = FALSE)

  if (plots) plot_bc_libraries(L, per_construct, cross, output_dir,
                                 paste0(file, ".png"), !is.null(ccl))
  if (verbose) {
    message(sprintf("%d libraries | %s cells carrying calls from more than one (candidate doublets)",
                    nrow(L), format(cross, big.mark = ",")))
    print(L[, c("library", "constructs_expected", "constructs_detected", "dropout",
                "cells", "cells_per_construct", "clones", "cv_cells",
                "top_decile_share", "gini_cells")])
  }
  invisible(L)
}

# Three-panel library comparison figure.
plot_bc_libraries <- function(L, per_construct, cross, dir, file, on_called) {
  pal <- stats::setNames(rep(c("#2E6E8E", "#D1A13A", "#B4432F", "#8FA8BF"),
                             length.out = nrow(L)), L$library)
  den <- if (on_called) "cells called by the chemistry's pipeline" else
    "every barcode observed, ambient floor included"

  # A. detection against each library's OWN size - absolute counts, with the
  # library size drawn as the reference rather than renormalising to 100%
  A <- rbind(
    data.frame(library = L$library, n = L$constructs_detected, what = "detected"),
    data.frame(library = L$library, n = L$dropout,             what = "not detected"))
  A$what <- factor(A$what, levels = c("detected", "not detected"))
  A$library <- factor(A$library, levels = L$library)
  lab <- data.frame(library = factor(L$library, levels = L$library),
                    n = L$constructs_expected,
                    txt = sprintf("%d / %d\n(%.0f%%)", L$constructs_detected,
                                  L$constructs_expected, L$pct_detected))
  pA <- ggplot2::ggplot(A, ggplot2::aes(library, n, fill = what)) +
    ggplot2::geom_col(width = .6) +
    ggplot2::geom_text(data = lab, ggplot2::aes(library, n, label = txt),
                       inherit.aes = FALSE, vjust = -0.25, size = 2.7, colour = "grey30") +
    ggplot2::scale_fill_manual(values = c(detected = "#2E6E8E",
                                          `not detected` = "#E3E3E3"), name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .22))) +
    ggplot2::labs(x = NULL, y = "constructs in the library",
                  title = "A. Constructs recovered",
                  subtitle = "bar height = that library's own size") +
    cat_theme() + ggplot2::theme(legend.position = "bottom")

  # B. evenness: the ranked cells-per-construct curve of each library
  pcs <- per_construct[order(per_construct$set, -per_construct$cells), ]
  pcs$rank <- stats::ave(pcs$cells, pcs$set, FUN = function(z) seq_along(z))
  pcs$set <- factor(pcs$set, levels = L$library)
  pB <- ggplot2::ggplot(pcs, ggplot2::aes(rank, cells + 1, colour = set)) +
    ggplot2::geom_step(linewidth = .7) +
    ggplot2::scale_colour_manual(values = pal, name = NULL) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "construct, ranked within its library",
                  y = "cells + 1 (log)", title = "B. Evenness",
                  subtitle = sprintf("CV %s | top decile holds %s",
                                     paste(L$cv_cells, collapse = " vs "),
                                     paste0(L$top_decile_share, "%", collapse = " vs "))) +
    cat_theme() + ggplot2::theme(legend.position = "bottom")

  # C. yield per construct - the fair way to compare pools of different size
  pC <- ggplot2::ggplot(L, ggplot2::aes(stats::reorder(library, -cells_per_construct),
                                        cells_per_construct, fill = library)) +
    ggplot2::geom_col(width = .55) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", cells_per_construct)),
                       vjust = -0.4, size = 2.8, colour = "grey30") +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .18))) +
    ggplot2::labs(x = NULL, y = "cells per construct",
                  title = "C. Yield per construct",
                  subtitle = sprintf("%s cells carry calls from both libraries",
                                     format(cross, big.mark = ","))) +
    cat_theme()

  fig <- patchwork::wrap_plots(pA, pB, pC, nrow = 1) +
    patchwork::plot_annotation(
      title = "The pooled shRNA libraries, compared within one run",
      subtitle = paste("denominator:", den),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 11, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 8, colour = "grey35")))
  cat_save_png(fig, dir, file, 11.0, 4.2)
  invisible(NULL)
}

## Part of catcheR. Comparison of the sublibraries inside one run.
##
## A sublibrary is a separate sequencing aliquot of the same cell pool: Parse
## splits a run into several, and the cell ID carries the tag. This asks whether
## the aliquots behaved alike, which is a different question from
## `catcheR_compare_libraries()` - that one compares the pooled shRNA libraries,
## this one compares the sequencing aliquots.

# Compare the sublibraries within a run
#
# Works on any run object that was given more than one FASTQ pair with
# `sublibraries`. For each aliquot: how much of its sequencing reached a
# barcode, how many cells it accounts for, how deep those cells are, and how
# many constructs and clones it saw.
#
# The one quantity that is not a per-aliquot repeat of the main QC is
# **clone sharing**: how many barcode-UCI clones appear in more than one
# aliquot. Aliquots of the same pool legitimately share clones - a clone is a
# population of cells, and its members can be distributed across aliquots - so
# a high overlap is expected and is reported as a number, not flagged. What
# would be worth attention is the opposite: near-zero sharing between aliquots
# of a pool that should be the same cells.
#
# @param x a run object or the path to its `.rds`.
# @param called optional called-cell list; when given, per-cell columns are
#   restricted to cells the platform's own pipeline called.
# @param output_dir where to write; created if absent.
# @param file basename for the `.csv` and `.png`.
# @param plots write the figure as well as the table.
# @param verbose print the table.
# @return a data frame, invisibly, with `shared_clones` attached as an
#   attribute; also written as `<file>.csv`.
compare_sublibraries <- function(x,
                                         called     = NULL,
                                         output_dir = NULL,
                                         file       = NULL,
                                         plots      = TRUE,
                                         verbose    = TRUE) {
  if (is.character(x)) x <- readRDS(x)
  stopifnot(is.list(x), !is.null(x$per_cell), !is.null(x$clone_summary))
  if (is.null(output_dir)) output_dir <- x$output_dir
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  pfx <- if (inherits(x, "tenxcatch")) "tenxcatch"
         else if (inherits(x, "parsecatch")) "parsecatch" else "pipcatch"
  if (is.null(file)) file <- paste0(pfx, "_sublibrary_comparison")

  ps <- x$per_sublibrary
  if (is.null(ps) || nrow(ps) < 2L)
    stop("This run has a single sequencing library, so there is nothing to ",
         "compare. Pass several FASTQ pairs with `sublibraries` to get ",
         "per-aliquot tables.", call. = FALSE)

  asn <- x$assignments
  asn$sublibrary <- cat_sublib_of(asn$cellID)
  pc <- x$per_cell
  pc$sublibrary <- cat_sublib_of(pc$cellID)
  ccl <- cat_read_called(called)
  if (!is.null(ccl)) pc <- pc[pc$cellID %in% ccl, , drop = FALSE]

  # a clone is (barcode, UCI); count the aliquots each one appears in
  cl <- unique(data.frame(clone = paste(asn$barcode, asn$UCI, sep = "_"),
                          sublibrary = asn$sublibrary, stringsAsFactors = FALSE))
  n_sub <- table(cl$clone)
  shared <- sum(n_sub > 1L)

  subs <- ps$sublibrary
  rows <- lapply(subs, function(s) {
    r  <- ps[ps$sublibrary == s, , drop = FALSE]
    a  <- asn[asn$sublibrary == s, , drop = FALSE]
    p  <- pc[pc$sublibrary == s, , drop = FALSE]
    own <- cl$clone[cl$sublibrary == s]
    data.frame(
      sublibrary          = s,
      read_pairs          = r$pairs,
      cellbc_reads        = r$cell_bc_reads,
      cassette_reads      = r$cassette_reads,
      pct_cassette        = round(100 * r$cassette_reads / r$pairs, 2),
      events              = nrow(a),
      cells               = nrow(p),
      cells_one_clone     = sum(p$True_UCIxcellID == 1, na.rm = TRUE),
      pct_one_clone       = round(100 * sum(p$True_UCIxcellID == 1, na.rm = TRUE) /
                                    max(1, nrow(p)), 1),
      median_depth        = if (nrow(p)) stats::median(p$reads_total) else NA_real_,
      constructs          = length(unique(a$shRNA)),
      clones              = length(own),
      clones_shared       = sum(own %in% names(n_sub)[n_sub > 1L]),
      stringsAsFactors    = FALSE)
  })
  S <- do.call(rbind, rows)
  S$pct_clones_shared <- round(100 * S$clones_shared / pmax(1, S$clones), 1)
  attr(S, "shared_clones") <- shared
  utils::write.csv(S, file.path(output_dir, paste0(file, ".csv")), row.names = FALSE)

  if (plots) plot_sublibraries(S, pc, shared, output_dir,
                                    paste0(file, ".png"), !is.null(ccl))
  if (verbose) {
    message(sprintf("%d sublibraries | %s of %s clones seen in more than one",
                    nrow(S), format(shared, big.mark = ","),
                    format(length(n_sub), big.mark = ",")))
    print(S[, c("sublibrary", "read_pairs", "pct_cassette", "cells",
                "pct_one_clone", "median_depth", "constructs", "clones",
                "pct_clones_shared")])
  }
  invisible(S)
}

# Three-panel sublibrary comparison figure.
plot_sublibraries <- function(S, pc, shared, dir, file, on_called) {
  lv  <- S$sublibrary
  pal <- stats::setNames(rep(c("#2E6E8E", "#D1A13A", "#B4432F", "#8FA8BF"),
                             length.out = length(lv)), lv)
  S$sublibrary <- factor(S$sublibrary, levels = lv)
  den <- if (on_called) "cells called by the chemistry's pipeline" else
    "every barcode observed, ambient floor included"

  # A. absolute cell counts, split by whether the cell resolves to one clone;
  # the aliquots differ in size, so the count is the honest bar and the share
  # is written on it
  A <- rbind(
    data.frame(sublibrary = S$sublibrary, n = S$cells_one_clone, what = "1 clone"),
    data.frame(sublibrary = S$sublibrary, n = S$cells - S$cells_one_clone,
               what = "other / no call"))
  A$what <- factor(A$what, levels = c("1 clone", "other / no call"))
  lab <- data.frame(sublibrary = S$sublibrary, n = S$cells,
                    txt = sprintf("%s cells\n%.0f%% one clone",
                                  format(S$cells, big.mark = ","), S$pct_one_clone))
  pA <- ggplot2::ggplot(A, ggplot2::aes(sublibrary, n, fill = what)) +
    ggplot2::geom_col(width = .6) +
    ggplot2::geom_text(data = lab, ggplot2::aes(sublibrary, n, label = txt),
                       inherit.aes = FALSE, vjust = -0.25, size = 2.6, colour = "grey30") +
    ggplot2::scale_fill_manual(values = c(`1 clone` = "#2E6E8E",
                                          `other / no call` = "#C7D3DC"), name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .24))) +
    ggplot2::labs(x = NULL, y = "cells", title = "A. Cells and calls",
                  subtitle = "bar height = that aliquot's own cell count") +
    cat_theme() + ggplot2::theme(legend.position = "bottom")

  # B. depth distribution per aliquot - the thing most likely to differ if one
  # aliquot was sequenced harder than the other
  pcf <- pc[pc$sublibrary %in% lv & pc$reads_total > 0, , drop = FALSE]
  pcf$sublibrary <- factor(pcf$sublibrary, levels = lv)
  med <- stats::aggregate(list(med = pcf$reads_total),
                          by = list(sublibrary = pcf$sublibrary), FUN = stats::median)
  pB <- ggplot2::ggplot(pcf, ggplot2::aes(sublibrary, reads_total, fill = sublibrary)) +
    ggplot2::geom_violin(scale = "width", linewidth = .2, colour = "white") +
    ggplot2::geom_boxplot(width = .13, outlier.shape = NA, linewidth = .3, fill = "white") +
    ggplot2::geom_text(data = med, ggplot2::aes(sublibrary, med,
                       label = format(med, big.mark = ",")),
                       inherit.aes = FALSE, hjust = -0.9, size = 2.6, colour = "grey25") +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = NULL, y = "cassette molecules per cell",
                  title = "B. Depth per cell", subtitle = "log scale") +
    cat_theme()

  # C. clone sharing. Aliquots of one pool are EXPECTED to share clones, so this
  # is reported, not flagged.
  pC <- ggplot2::ggplot(S, ggplot2::aes(sublibrary, clones, fill = sublibrary)) +
    ggplot2::geom_col(width = .55) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%s clones\n%.0f%% shared",
                                     format(clones, big.mark = ","), pct_clones_shared)),
                       vjust = -0.25, size = 2.6, colour = "grey30") +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .26))) +
    ggplot2::labs(x = NULL, y = "distinct barcode-UCI clones",
                  title = "C. Clones, and how many are shared",
                  # the long form of this ran off the panel; the reasoning is
                  # in the header of this file and in the template text
                  subtitle = sprintf("%s clones in >1 aliquot (expected)",
                                     format(shared, big.mark = ","))) +
    cat_theme()

  fig <- patchwork::wrap_plots(pA, pB, pC, nrow = 1) +
    patchwork::plot_annotation(
      title = "The sequencing sublibraries, compared within one run",
      subtitle = paste("denominator:", den),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 11, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 8, colour = "grey35")))
  cat_save_png(fig, dir, file, 11.0, 4.2)
  invisible(NULL)
}
