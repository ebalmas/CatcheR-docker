# Parse (Evercode) internals.
#
# Parse differs from PIPseq in three ways that matter here:
#   * the cassette is on READ 1, in reverse complement, after one leading base;
#   * the cell barcode is on READ 2, as three ligation rounds separated by two
#     fixed 12 nt linkers, each round drawn from a published whitelist;
#   * read 2 opens with a real 10 nt UMI, so events can be deduplicated by
#     molecule rather than counted as reads.
# Everything downstream of extraction is the shared .cat_ layer.

#' Plate well label to Parse's well index: A1 = 1, A12 = 12, B1 = 13, E11 = 59.
.par_well_index <- function(w) {
  row <- match(toupper(substr(w, 1, 1)), LETTERS)
  (row - 1L) * 12L + as.integer(substr(w, 2, nchar(w)))
}

#' Read split-pipe's `barcode_data.csv`.
#'
#' The file is not a plain CSV: `#` comment lines carry the kit and chemistry
#' and also separate the rounds, and the single `bci,sequence,...` header sits
#' after the first of them rather than at the top, so `read.csv()` fails on it.
#' The round each row belongs to is read from the `bc_round` column, not from
#' the comment banners. Round 1 holds TWO sequences per well - one oligo-dT
#' and one random hexamer (`stype` T and R) - and both must map to the SAME
#' well, because that is how split-pipe assigns cells. Rounds 2 and 3 hold one
#' ligation barcode per well.
#'
#' @param path path to `barcode_data.csv`.
#' @return a list with `table` (the parsed rows), `maps` (one named
#'   sequence -> well-index vector per round) and `header` (the comment lines).
.par_whitelist <- function(path) {
  if (!file.exists(path))
    stop("Parse barcode file not found: ", path,
         "\n  This is split-pipe's barcode_data.csv; it is required because a ",
         "whitelist guessed from the data misses the random-hexamer barcodes.",
         call. = FALSE)
  ln  <- readLines(path, warn = FALSE)
  hdr <- ln[startsWith(ln, "#")]
  ln  <- ln[!startsWith(ln, "#") & nzchar(trimws(ln))]
  cols <- strsplit(ln[startsWith(ln, "bci,")][1], ",")[[1]]
  dat  <- ln[!startsWith(ln, "bci,")]
  m <- do.call(rbind, strsplit(dat, ",", fixed = TRUE))
  df <- as.data.frame(m, stringsAsFactors = FALSE); names(df) <- cols
  df$bci <- as.integer(df$bci); df$bc_round <- as.integer(df$bc_round)
  df$well_index <- .par_well_index(df$well)
  maps <- lapply(sort(unique(df$bc_round)), function(r) {
    d <- df[df$bc_round == r, ]
    stats::setNames(d$well_index, toupper(d$sequence))
  })
  names(maps) <- paste0("round", sort(unique(df$bc_round)))
  list(table = df, maps = maps, header = hdr)
}

#' Expand a whitelist into a correction map.
#'
#' Exact sequences always map to themselves. Neighbours within
#' `max.edit.distance` substitutions are added only when unambiguous: a variant
#' reachable from two whitelist entries with different well indices is dropped
#' rather than assigned arbitrarily. split-pipe itself corrects out to distance
#' 2, so that is available, but distance 1 is the default because the ambiguous
#' fraction grows quickly.
.par_corr_map <- function(seq2idx, max.edit.distance = 1L) {
  stopifnot(max.edit.distance %in% c(0L, 1L, 2L))
  if (max.edit.distance == 0L) return(seq2idx)
  bases <- c("A", "C", "G", "T")
  seqs  <- names(seq2idx)
  vars  <- character(0); owner <- integer(0)
  for (si in seq_along(seqs)) {
    ch <- strsplit(seqs[si], NULL)[[1]]; L <- length(ch); idx <- seq2idx[[si]]
    pos1 <- unlist(lapply(seq_len(L), function(i) {
      b <- setdiff(bases, ch[i])
      vapply(b, function(x) { z <- ch; z[i] <- x; paste(z, collapse = "") }, character(1))
    }), use.names = FALSE)
    v <- pos1
    if (max.edit.distance == 2L) {
      pr <- utils::combn(L, 2)
      v2 <- unlist(lapply(seq_len(ncol(pr)), function(k) {
        i <- pr[1, k]; j <- pr[2, k]
        unlist(lapply(setdiff(bases, ch[i]), function(a)
          vapply(setdiff(bases, ch[j]), function(b) {
            z <- ch; z[i] <- a; z[j] <- b; paste(z, collapse = "")
          }, character(1))), use.names = FALSE)
      }), use.names = FALSE)
      v <- c(v, v2)
    }
    vars <- c(vars, v); owner <- c(owner, rep(idx, length(v)))
  }
  keep <- !(vars %in% seqs)                       # exact entries win outright
  vars <- vars[keep]; owner <- owner[keep]
  first <- !duplicated(vars)
  amb   <- unique(vars[duplicated(vars) &
                       owner[match(vars, vars)] != owner])   # conflicting claims
  v <- vars[first]; o <- owner[first]
  drop <- v %in% amb
  c(seq2idx, stats::setNames(o[!drop], v[!drop]))
}

#' Parse adapter for [.cat_stream()].
#'
#' Read 2 is `[UMI][bc3][linker A][bc2][linker B][bc1][poly-T]`, with the
#' rounds running 3-2-1 from the 5' end because round 3 ligates last. The
#' returned cell ID is split-pipe's own `bc_wells` string, well indices and
#' sublibrary suffix included, so the join to `cell_metadata.csv` is a plain
#' string match with no translation step.
.par_adapter <- function(cmaps, sublibrary, linker.a, linker.b,
                         umi.length = 10L, barcode.length = 8L, lead = 1L) {
  UL <- as.integer(umi.length); BL <- as.integer(barcode.length)
  LA <- nchar(linker.a); LB <- nchar(linker.b)
  p3 <- UL + 1L;            p3e <- UL + BL
  la <- p3e + 1L;           lae <- p3e + LA
  p2 <- lae + 1L;           p2e <- lae + BL
  lb <- p2e + 1L;           lbe <- p2e + LB
  p1 <- lbe + 1L;           p1e <- lbe + BL
  list(name = "parse", cell_read = 2L, has_umi = TRUE, lead = as.integer(lead),
       offsets = c(umi = 1L, bc3 = p3, bc2 = p2, bc1 = p1),
       cell_fun = function(s) {
         okl <- substr(s, la, lae) == linker.a & substr(s, lb, lbe) == linker.b
         i1 <- unname(cmaps$round1[substr(s, p1, p1e)])
         i2 <- unname(cmaps$round2[substr(s, p2, p2e)])
         i3 <- unname(cmaps$round3[substr(s, p3, p3e)])
         ok <- okl & !is.na(i1) & !is.na(i2) & !is.na(i3)
         cell <- rep(NA_character_, length(s))
         cell[ok] <- sprintf("%02d_%02d_%02d__%s", i1[ok], i2[ok], i3[ok], sublibrary)
         list(valid = ok, cell = cell, umi = substr(s, 1L, UL),
              linker_ok = okl)
       })
}

#' Read a Parse DGE directory written by split-pipe.
#'
#' `count_matrix.mtx` is stored CELLS x GENES, the transpose of the 10X and
#' PIPseeker convention, so it is transposed on load; `Seurat::ReadMtx()`
#' cannot be used for the same reason. Genes are named from `gene_name` where
#' present and `gene_id` otherwise.
.par_read_dge <- function(mtx, genes, cells) {
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Matrix is required to read the Parse count matrix.", call. = FALSE)
  m  <- Matrix::readMM(mtx)
  gn <- utils::read.csv(genes, stringsAsFactors = FALSE)
  cm <- utils::read.csv(cells, stringsAsFactors = FALSE)
  if (nrow(m) != nrow(cm) || ncol(m) != nrow(gn))
    stop(sprintf(paste0("Parse matrix is %d x %d but cell_metadata has %d rows and ",
                        "all_genes %d rows. Expected cells x genes."),
                 nrow(m), ncol(m), nrow(cm), nrow(gn)), call. = FALSE)
  nm <- ifelse(nzchar(gn$gene_name) & !is.na(gn$gene_name), gn$gene_name, gn$gene_id)
  m <- Matrix::t(m)
  rownames(m) <- make.unique(nm); colnames(m) <- cm$bc_wells
  list(counts = m, meta = cm, genes = gn)
}
