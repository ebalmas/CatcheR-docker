## Part of catcheR. Threshold diagnostics, shared by every chemistry.
##
## These exist because yield is the wrong criterion for choosing a threshold.
## Maximising the number of cells that receive a call argues for no filter at
## all; what a threshold is actually for is keeping the calls CORRECT. The two
## functions here measure that directly.

# ── catcheR_call_accuracy ───────────────────────────────────────────────────

#' How deep does a cell have to be before its clone call can be trusted?
#'
#' Takes the deep, unambiguous cells as ground truth, subsamples each one's UMI
#' distribution down to a target depth, re-applies catcheR's calling rule, and
#' scores the outcome as the correct clone, the WRONG clone, or no call. The
#' wrong-call rate among the calls actually made is the number that should set
#' the per-cell floor; the no-call rate only costs yield.
#'
#' @section Reading the result:
#' A non-monotonic `correct` column at the shallow end is not noise. When a
#' cell holds fewer than `100 / percentage` UMIs, every UCI supported by a
#' single UMI already clears the percentage cut, so several UCIs pass and the
#' cell is rejected as ambiguous. Raising `UMI.count` to 1 removes most of it.
#'
#' @param x a `pipcatch` or `parsecatch` object, or the path to its `.rds`.
#' @param depths target depths to test.
#' @param min.truth.UMI minimum UMIs for a cell to serve as ground truth.
#' @param UMI.count,percentage,ratio the rule to test. Defaults to the object's
#'   own parameters.
#' @param n.rep replicate subsamplings per depth; results are averaged.
#' @param seed RNG seed, for a reproducible table.
#' @param output_dir where to write the CSV. `NULL` uses the object's.
#' @param verbose print the table.
#' @return a data frame, one row per depth, also written as a CSV.
#' @export
catcheR_call_accuracy <- function(x,
                                  depths        = c(2, 3, 5, 10, 20, 30, 50, 100),
                                  min.truth.UMI = 100,
                                  UMI.count     = NULL,
                                  percentage    = NULL,
                                  ratio         = NULL,
                                  n.rep         = 1L,
                                  seed          = 1L,
                                  output_dir    = NULL,
                                  verbose       = TRUE) {
  if (is.character(x)) x <- readRDS(x)
  stopifnot(inherits(x, "pipcatch"))
  if (is.null(UMI.count))  UMI.count  <- x$params$UMI.count
  if (is.null(percentage)) percentage <- x$params$percentage
  if (is.null(ratio))      ratio      <- x$params$ratio
  if (is.null(output_dir)) output_dir <- x$output_dir

  ## ground truth: deep cells that resolve to exactly one UCI at full depth
  v <- .cat_validate(x$tab, UMI.count, percentage, ratio)
  deep <- v$cells$cellID[v$cells$True_UCIxcellID == 1 &
                         x$tab$UMIxcellID[match(v$cells$cellID, x$tab$cellID)] >= min.truth.UMI]
  if (length(deep) < 50)
    stop(sprintf(paste0("Only %d cells reach %d UMIs, too few to serve as ground truth. ",
                        "Lower min.truth.UMI, or accept that this library is too shallow ",
                        "for the test."), length(deep), min.truth.UMI), call. = FALSE)
  sub <- x$tab[x$tab$cellID %in% deep, c("cellID", "UCI", "UMIxUCI")]
  sp  <- split(sub, sub$cellID)
  tru <- vapply(sp, function(d) d$UCI[which.max(d$UMIxUCI)], character(1))

  set.seed(seed)
  res <- bind_rows(lapply(depths, function(n) {
    ok <- wrong <- none <- 0L
    for (rep in seq_len(n.rep)) for (i in seq_along(sp)) {
      d <- sp[[i]]
      draw <- stats::rmultinom(1, n, d$UMIxUCI / sum(d$UMIxUCI))[, 1]
      # Mirror .cat_validate: strict `>` on both criteria, then the
      # dominant-to-second ratio with catcheR's 100%-share escape clause.
      # Scoring the subsample on a looser rule than the filter makes the error
      # curve optimistic about the very threshold it is used to choose.
      pass <- which(draw > UMI.count & 100 * draw / n > percentage)
      okratio <- FALSE
      if (length(pass) == 1L) {
        srt <- sort(draw, decreasing = TRUE)
        okratio <- (100 * draw[pass] / n == 100) || length(srt) < 2L ||
          srt[2] == 0 || (srt[1] / srt[2] >= ratio)
      }
      if (length(pass) != 1L || !okratio)  none  <- none  + 1L
      else if (d$UCI[pass] == tru[[i]])    ok    <- ok    + 1L
      else                                 wrong <- wrong + 1L
    }
    tot <- ok + wrong + none
    data.frame(UMIs_in_cell = n, reference_cells = length(sp),
               correct = round(100 * ok / tot, 1),
               wrong   = round(100 * wrong / tot, 1),
               no_call = round(100 * none / tot, 1),
               error_among_calls = round(100 * wrong / max(ok + wrong, 1), 2))
  }))
  attr(res, "params") <- list(UMI.count = UMI.count, percentage = percentage,
                              ratio = ratio, min.truth.UMI = min.truth.UMI,
                              n.rep = n.rep, seed = seed)
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(res, file.path(output_dir,
                     sprintf("catcheR_call_accuracy_%s.csv", x$sample)), row.names = FALSE)
  }
  if (verbose) {
    message(sprintf("Ground truth: %s cells with >=%d UMIs and one clone at UMI.count=%s, percentage=%s",
                    format(length(sp), big.mark = ","), min.truth.UMI, UMI.count, percentage))
    print(as.data.frame(res))
  }
  res
}

# ── catcheR_threshold_sweep ─────────────────────────────────────────────────

#' Sweep the calling thresholds and report both ways a cell can fail
#'
#' A cell without a call has failed in one of two opposite ways: no UCI cleared
#' the cut, or more than one did. Collapsing them into "no call" hides what a
#' threshold change actually does - lowering `percentage` mostly moves cells
#' from the first class into the second rather than into a call.
#'
#' @param x a `pipcatch` or `parsecatch` object, or path to its `.rds`.
#' @param percentage,UMI.count,min.cell.UMI values to sweep. Each is crossed
#'   with the others.
#' @param ratio held fixed across the sweep.
#' @param output_dir,verbose as elsewhere.
#' @return a data frame with one row per combination, also written as a CSV.
#' @export
catcheR_threshold_sweep <- function(x,
                                    percentage   = c(10, 12, 15, 18, 20),
                                    UMI.count    = c(0, 1, 2),
                                    min.cell.UMI = c(0, 5, 10, 20, 50),
                                    ratio        = NULL,
                                    output_dir   = NULL,
                                    verbose      = TRUE) {
  if (is.character(x)) x <- readRDS(x)
  stopifnot(inherits(x, "pipcatch"))
  if (is.null(ratio))      ratio      <- x$params$ratio
  if (is.null(output_dir)) output_dir <- x$output_dir
  depth <- x$tab %>% group_by(cellID) %>%
    summarise(UMIxcellID = first(UMIxcellID), .groups = "drop")

  grid <- expand.grid(percentage = percentage, UMI.count = UMI.count,
                      min.cell.UMI = min.cell.UMI)
  res <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    # Call the filter itself instead of restating its rule. This used to
    # re-implement it with `>=` on the percentage - the filter uses a strict
    # `>`, as catcheR does - and without the ratio test, so the sweep
    # advertised a yield the filter would not deliver: on AB019c it reported
    # 3,676 where the filter validates 3,500.
    v <- .cat_validate(x$tab, g$UMI.count, g$percentage, ratio)
    cells_in <- depth %>% filter(UMIxcellID >= g$min.cell.UMI) %>%
      left_join(v$cells %>% select(cellID, True_UCIxcellID, ValidCell), by = "cellID") %>%
      mutate(True_UCIxcellID = tidyr::replace_na(True_UCIxcellID, 0L),
             ValidCell       = tidyr::replace_na(ValidCell, FALSE))
    data.frame(percentage = g$percentage, UMI.count = g$UMI.count, ratio = ratio,
               min.cell.UMI = g$min.cell.UMI, cells_kept = nrow(cells_in),
               one_clone  = sum(cells_in$True_UCIxcellID == 1),
               valid      = sum(cells_in$ValidCell),
               ambiguous  = sum(cells_in$True_UCIxcellID > 1),
               none_pass  = sum(cells_in$True_UCIxcellID == 0),
               pct_valid  = round(100 * mean(cells_in$ValidCell), 1))
  }))
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(res, file.path(output_dir,
                     sprintf("catcheR_threshold_sweep_%s.csv", x$sample)), row.names = FALSE)
  }
  if (verbose) {
    best <- res[which.max(res$valid), ]
    message(sprintf("Highest yield: percentage=%s UMI.count=%s cell floor=%s -> %s ValidCell",
                    best$percentage, best$UMI.count, best$min.cell.UMI,
                    format(best$valid, big.mark = ",")))
    message("Yield alone is not the criterion - check catcheR_call_accuracy() before adopting it.")
  }
  res
}
