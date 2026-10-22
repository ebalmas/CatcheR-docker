# Shared, chemistry-independent internals.
#
# Everything here is used by every chemistry the package supports. The parts
# that differ between chemistries are confined to an "adapter" (see
# [.pip_adapter()] and [.par_adapter()]), which answers three questions:
# which read carries the cell barcode, how to parse it, and whether the
# molecular index is a real UMI.

#' Classify the cassette read. The empty reference MUST be tested before the constant:
#' the no-insertion junction regenerates the constant 21 nt early (measured on
#' AB019b: 1-based position 17 in empty reads vs 38 in cassette reads - the
#' cassette length, 8 + 6 + 7), so a constant-first classifier files
#' empty-vector reads as "barcode not in reference".
#'
#' `lead` is the number of bases that may precede the reference on the cassette
#' read. PIPseq reads the cassette from base 1 (`lead = 0`); Parse carries one
#' leading base before `rc(reference)`, so the anchor and the empty-vector test
#' both have to tolerate the shift.
.cat_classify <- function(s2, pos, ANCH, RC_EMPTY, lead = 0L) {
  sim <- integer(length(s2))
  off <- if (lead > 0L) {
    a <- regexpr(substr(ANCH, 1, 12), substr(s2, 1, lead + 13), fixed = TRUE)
    ifelse(a > 0L, a - 1L, 0L)
  } else integer(length(s2))
  for (i in seq_len(nchar(ANCH)))
    sim <- sim + (substr(s2, off + i, off + i) == substr(ANCH, i, i))
  polyA <- grepl(strrep("A", 10), substr(s2, 1, 45), fixed = TRUE)
  epos <- if (lead > 0L) regexpr(RC_EMPTY, substr(s2, 1, lead + nchar(RC_EMPTY)), fixed = TRUE)
          else ifelse(startsWith(s2, RC_EMPTY), 1L, -1L)
  ifelse(epos > 0L, "empty vector (no insertion)",
  ifelse(pos >= 15,                "cassette",
  ifelse(sim >= nchar(ANCH) - 2,
         ifelse(polyA, "primer + poly-A chimera", "primer, no cassette, no poly-A"),
         ifelse(polyA, "non-construct insert (poly-A)", "non-construct insert"))))
}

#' Stream a read pair and pull out the cassette, one chunk at a time.
#'
#' Chemistry enters only through `adapter`: which read carries the cell
#' barcode (`cell_read`), how to parse it (`cell_fun`, vectorised over a
#' character vector, returning `valid` / `cell` / `umi`), and how many bases
#' may precede the reference on the cassette read (`lead`). The cassette slice
#' itself is identical everywhere - UCI then shRNA barcode, both reverse
#' complemented, immediately 5' of the vector constant.
.cat_stream <- function(f1, f2, adapter, ANCH, RC_EMPTY, CONST, UL, BL,
                        chunk, max.reads, verbose) {
  c1 <- gzfile(f1, "r"); c2 <- gzfile(f2, "r"); on.exit({ close(c1); close(c2) })
  tot <- 0L; done <- 0; k <- 0L; valid_cell <- 0L
  lead <- if (is.null(adapter$lead)) 0L else as.integer(adapter$lead)
  cls_tot <- numeric(0); cell_cls <- list(); cassette <- list(); empty <- list()
  repeat {
    l1 <- readLines(c1, n = chunk * 4L); if (!length(l1)) break
    l2 <- readLines(c2, n = chunk * 4L)
    sA <- l1[seq(2, length(l1), 4)]; sB <- l2[seq(2, length(l2), 4)]
    if (adapter$cell_read == 1L) { s1 <- sA; s2 <- sB } else { s1 <- sB; s2 <- sA }
    tot <- tot + length(s1); k <- k + 1L
    pr <- adapter$cell_fun(s1); keep <- which(pr$valid)
    # reads whose cell barcode parsed and corrected; everything downstream is
    # conditioned on this, so without the counter the read funnel has a hole
    # between 'read pairs' and the first class count.
    valid_cell <- valid_cell + length(keep)
    if (length(keep)) {
      v2  <- s2[keep]; pos <- regexpr(CONST, v2, fixed = TRUE)
      cls <- .cat_classify(v2, pos, ANCH, RC_EMPTY, lead); cell <- pr$cell[keep]
      t <- tapply(rep(1, length(cls)), cls, sum)
      u <- union(names(cls_tot), names(t))
      cls_tot <- setNames(
        replace(numeric(length(u)), match(names(cls_tot), u), cls_tot) +
        replace(numeric(length(u)), match(names(t), u), as.numeric(t)), u)
      cell_cls[[k]] <- as.data.frame(table(cell = cell, cls = cls), stringsAsFactors = FALSE) %>%
                       dplyr::filter(Freq > 0)
      j <- which(cls == "cassette")
      if (length(j)) cassette[[k]] <- data.frame(
        cellID = cell[j], index = pr$umi[keep][j],
        UCI_rc = substr(v2[j], pos[j] - UL - BL, pos[j] - BL - 1L),
        bc_rc  = substr(v2[j], pos[j] - BL,      pos[j] - 1L), stringsAsFactors = FALSE)
      e <- which(cls == "empty vector (no insertion)")
      if (length(e)) empty[[k]] <- data.frame(
        cellID = cell[e], index = pr$umi[keep][e], stringsAsFactors = FALSE)
    }
    done <- done + length(s1)
    if (verbose) message(sprintf("  ... %s read pairs", format(done, big.mark = ",")))
    if (done >= max.reads) break
  }
  list(total_pairs = tot, valid_cell_reads = valid_cell, class_counts = cls_tot,
       cell_class = bind_rows(cell_cls) %>% count(cell, cls, wt = Freq, name = "n"),
       cassette   = bind_rows(cassette),
       empty      = bind_rows(empty))
}

#' Copy the files that feed a later step into a shared scratch folder
#'
#' A `.csv` written into `Output/` is a deliverable; a `.csv` that the next
#' step reads is an input, and in catcheR's own layout those belong in
#' `scratch/`. Two of them - `cells_with_at_least_1_UCI.csv` and
#' `new_cell_names.csv` - carry neither the chemistry nor the sample in their
#' name, which is harmless inside one chemistry's folder and wrong the moment
#' two pipelines share a scratch. Copying is deliberate and the copy is
#' renamed; the originals stay where they were written.
#'
#' @param files named character vector: name = file in `from`, value = base
#'   name to use in `to` (without the extension).
#' @return the paths written, invisibly.
.cat_to_scratch <- function(files, from, to, chem, sample, verbose = TRUE) {
  dir.create(to, showWarnings = FALSE, recursive = TRUE)
  out <- character(0)
  for (f in names(files)) {
    src <- file.path(from, f)
    if (!file.exists(src)) next
    ext <- sub("^.*\\.", "", f)
    dst <- file.path(to, sprintf("%s_%s_%s.%s", files[[f]], chem, sample, ext))
    file.copy(src, dst, overwrite = TRUE)
    out <- c(out, dst)
  }
  if (verbose && length(out))
    message(sprintf("  to scratch (%s): %s", to,
                    paste(basename(out), collapse = ", ")))
  invisible(out)
}

#' Record which file played which role in a run
#'
#' Answers "where did this result come from" without trusting a folder name.
#' Roles are explicit, paths are the ones actually opened, and size plus
#' modification time make an accidental substitution visible. This is an input
#' RECORD, not a checksum: it cannot prove the biological identity of a sample.
#'
#' @param roles named character vector: name = role, value = path (a role may
#'   hold several paths, e.g. two FASTQ pairs).
#' @return a data.frame, also written as `<sample>_<chem>_input_files.csv`.
.cat_input_manifest <- function(roles, sample, chem, dir) {
  roles <- roles[!vapply(roles, is.null, logical(1))]
  rows <- do.call(rbind, lapply(names(roles), function(r) {
    p <- as.character(roles[[r]])
    data.frame(role = r, file = p,
               exists = file.exists(p),
               bytes = ifelse(file.exists(p), file.size(p), NA_real_),
               modified = as.character(ifelse(file.exists(p),
                                              format(file.mtime(p), "%Y-%m-%d %H:%M"), NA)),
               stringsAsFactors = FALSE)
  }))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(rows, file.path(dir, sprintf("%s_%s_input_files.csv", sample, chem)),
                   row.names = FALSE)
  rows
}
