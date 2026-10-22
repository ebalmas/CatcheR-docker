# 10x Genomics internals.
#
# The simplest of the three chemistries: read 1 is [cell barcode][UMI] at
# fixed offsets, with no linkers and no phasing, which is why the original
# `catcheR_10Xcatch()` could slice it with `cut -c 1-16` and `cut -c 17-28`.
# Read 2 carries the cassette exactly as PIPseq does, so everything from the
# reference anchor onwards is the shared .cat_ layer untouched.
#
# The one thing 10x does NOT give us offline is its barcode whitelist (the
# 737K / 3M lists ship with Cell Ranger, not with the data). Correction is
# therefore opt-in: pass `cell.whitelist`, typically the called-cell list from
# `filtered_feature_bc_matrix.h5`, and barcodes are accepted (and corrected at
# `max.edit.distance`) against it. With no whitelist every structurally valid
# 16-mer is accepted, which inflates the observed-barcode count with the
# ambient floor and with sequencing errors - see the warning in
# `catcheR_tenxcatch()`.

#' Read-1 parser for 10x v2/v3: fixed-width cell barcode then UMI.
#'
#' @param s character vector of read-1 sequences.
#' @param cmap optional named integer/character vector mapping an observed
#'   barcode to its accepted form, as built by [.par_corr_map()].
#' @return list(valid, cell, umi), vectorised over `s`.
.tenx_parse_r1 <- function(s, cb.length = 16L, umi.length = 12L, cmap = NULL) {
  need <- cb.length + umi.length
  cell <- substr(s, 1L, cb.length)
  umi  <- substr(s, cb.length + 1L, need)
  ok   <- nchar(s) >= need & !grepl("N", cell, fixed = TRUE)
  if (!is.null(cmap)) {
    hit  <- cmap[cell]
    ok   <- ok & !is.na(hit)
    cell <- ifelse(is.na(hit), cell, as.character(hit))
  }
  list(valid = ok, cell = cell, umi = umi)
}

#' Adapter describing 10x to the shared streamer.
.tenx_adapter <- function(cb.length = 16L, umi.length = 12L, cmap = NULL) {
  list(name      = "tenx",
       cell_read = 1L,    # cell barcode on read 1, cassette on read 2
       has_umi   = TRUE,  # a real 12 nt UMI: events are molecules, not reads
       lead      = 0L,    # cassette starts at the very start of read 2
       cell_fun  = function(s)
         .tenx_parse_r1(s, cb.length, umi.length, cmap))
}

#' Read the called-cell barcodes out of a Cell Ranger HDF5 matrix.
#'
#' Accepts `filtered_feature_bc_matrix.h5`, a `barcodes.tsv[.gz]`, or a plain
#' text file of barcodes. The `-1` GEM-well suffix Cell Ranger appends is
#' stripped, because the barcode this package reads off read 1 does not carry
#' it.
.tenx_called <- function(path) {
  if (is.null(path)) return(NULL)
  if (grepl("\\.h5$", path)) {
    if (!requireNamespace("hdf5r", quietly = TRUE))
      stop("Reading a Cell Ranger .h5 needs the 'hdf5r' package. Install it, or pass ",
           "barcodes.tsv.gz instead.", call. = FALSE)
    h <- hdf5r::H5File$new(path, mode = "r"); on.exit(h$close_all())
    nm <- if (h$exists("matrix/barcodes")) "matrix/barcodes" else "GRCh38/barcodes"
    bc <- h[[nm]]$read()
  } else {
    bc <- readLines(if (grepl("\\.gz$", path)) gzfile(path) else path)
  }
  unique(sub("-[0-9]+$", "", as.character(bc)))
}