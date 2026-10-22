#' Sanger sequencing QC for shRNA plasmid colony PCR
#'
#' Analyses Sanger sequencing reads (FASTA format) from colony PCR of an
#' intermediate shRNA plasmid library. For each well, the function:
#'
#' 1. Finds the MCS anchor sequence `GGCGCGCCATTTAAATGTCGAC` (RC of
#'    `GTCGACATTTAAATGGCGCGCC`) with fuzzy matching.
#' 2. Extracts the barcode (8 bp immediately before the anchor) and the
#'    UCI (6 bp before the barcode).
#' 3. Extracts the shRNA sense and antisense regions from the other side of
#'    the anchor (past the polyA tail and CTCGAG loop).
#' 4. Matches BC and shRNA against all entries in the design plate layout,
#'    using hamming-distance tolerance.
#' 5. Reads per-base quality from ABI PHD files (if provided) and reports
#'    minimum and mean Q scores in the BC and shRNA regions.
#' 6. Calls each well as PASS, CHECK, INVALID HAIRPIN, etc.
#'
#' **Read orientation.** The sequencing primer reads the bottom strand
#' 5'→3', so landmarks appear in the order:
#' `MluI → BC → UCI → GGCGCGCC·ATTTAAAT·GTCGAC → polyA → sense → CTCGAG →
#' antisense`.
#' The barcode in the read equals the **RC of the design BC column**.
#'
#' @param fasta        Path to the FASTA file of Sanger sequences (one
#'   entry per well, headers ending in the well ID e.g. `_A1`).
#' @param design_xlsx  Path to the design Excel file. Must contain a sheet
#'   named as `plate_sheet` with columns: position (well), oligo name,
#'   and the TOP OLIGO sequence.
#' @param output_dir   Output directory. Default `"Output/"`.
#' @param sample_name  Optional sample name added to the run folder name.
#' @param plate_sheet  Name of the Excel sheet containing the plate layout.
#'   Default `"Plate 01"`.
#' @param oligo_col    0-based column index of the oligo name in
#'   `plate_sheet`. Default `2`.
#' @param sequence_col 0-based column index of the TOP OLIGO sequence.
#'   Default `4`.
#' @param first_row    0-based row index of the first oligo (after headers).
#'   Default `8`.
#' @param last_row     0-based row index of the last oligo (inclusive).
#'   Default `67`.
#' @param ligation_split Column number (1-based) at which ligation condition
#'   changes. Wells in columns 1:`ligation_split` = condition A (e.g.
#'   `"1:12"`), columns above = condition B (e.g. `"1:20"`).
#'   Default `6`.
#' @param ligation_labels Character vector of length 2 naming the two
#'   ligation conditions. Default `c("1:12", "1:20")`.
#' @param anchor       The anchor sequence to search for in the read.
#'   Default `"GGCGCGCCATTTAAATGTCGAC"`.
#' @param anchor_mm    Maximum mismatches allowed when locating the anchor.
#'   Default `2`.
#' @param bc_mm_pass   Maximum BC mismatches to still call a PASS.
#'   Default `2`.
#' @param shrna_mm_pass Maximum total shRNA mismatches (sense + antisense)
#'   to still call a PASS. Default `4`.
#' @param min_read_len Minimum read length to attempt extraction.
#'   Default `150`.
#' @param abi_zip      Optional path to a ZIP file containing ABI/PHD files
#'   for per-base quality scores. PHD files must be named
#'   `<prefix>_<well>.phd.1`. Default `NULL`.
#' @param phd_prefix   Prefix used in PHD filenames inside `abi_zip`
#'   (everything before `_<well>.phd.1`). Default `NULL` (auto-detected
#'   from the first PHD file found).
#' @param sh_id_prefix Prefix to strip from oligo names to produce the
#'   short shRNA ID (e.g. `"siKDBC_"`). Default `"siKDBC_"`.
#'
#' @return A named list (invisibly) with:
#' \describe{
#'   \item{`results`}{Data frame with one row per Sanger well.}
#'   \item{`design`}{Data frame of the parsed design plate layout.}
#'   \item{`summary`}{Named integer vector of status counts.}
#'   \item{`paths`}{List of output paths.}
#' }
#'
#' @section Output files written to `output_dir/sangerQC/<run>/`:
#' \describe{
#'   \item{`csv/sanger_QC_results.csv`}{Full per-well result table.}
#'   \item{`stats/sanger_QC_summary.txt`}{Summary counts and key findings.}
#'   \item{`R_objects/sangerQC_result.rds`}{Full result list.}
#'   \item{`to_scratch/sanger_QC_results.csv`}{Copy for downstream use.}
#' }
#'
#' @section Status codes:
#' \describe{
#'   \item{`PASS`}{BC and shRNA both match perfectly (0 mm each).}
#'   \item{`PASS (Xbp BC mut)`}{shRNA perfect, BC has 1–2 mismatches.}
#'   \item{`PASS (Xbp shRNA mut)`}{BC perfect, shRNA has 1–2 total mm.}
#'   \item{`CHECK`}{BC matches (0 mm) but shRNA has > `shrna_mm_pass` mm —
#'     likely a real clone but inspect the sequence manually.}
#'   \item{`INVALID HAIRPIN`}{Sense ≠ RC(antisense) — possible sequencing
#'     artefact in the shRNA hairpin region.}
#'   \item{`shRNA OK / BC WRONG`}{shRNA perfect but BC > `bc_mm_pass` mm.}
#'   \item{`BC/shRNA discordant`}{BC and shRNA each match, but different
#'     library members.}
#'   \item{`SEQ TRUNCATED`}{Anchor not found in the read.}
#'   \item{`SEQ FAILED`}{Read shorter than `min_read_len`.}
#' }
#'
#' @section Required packages:
#' `readxl` for reading the Excel design file. Install with
#' `install.packages("readxl")`.
#'
#' @examples
#' \dontrun{
#' result <- catcheR_sangerQC(
#'   fasta        = "11109801799-1.fasta",
#'   design_xlsx  = "EB003C_DES01_shRNA_tdark_Ligation.xlsx",
#'   output_dir   = "Output/",
#'   sample_name  = "plate1_ligation",
#'   plate_sheet  = "Plate 01",
#'   abi_zip      = "11109801799-1_SCF_SEQ_ABI.zip"
#' )
#'
#' # Inspect results
#' result$results
#' result$summary
#'
#' # All PASS wells
#' result$results[grepl("^PASS", result$results$status), ]
#' }
#'
#' @seealso [catcheR_step1QC()], [catcheR_step2QC_extraction()]
#' @export
catcheR_sangerQC <- function(fasta,
                              design_xlsx,
                              output_dir       = "Output/",
                              sample_name      = NULL,
                              plate_sheet      = "Plate 01",
                              oligo_col        = 1L,
                              sequence_col     = 3L,
                              first_row        = 8L,
                              last_row         = 67L,
                              ligation_split   = 6L,
                              ligation_labels  = c("1:12", "1:20"),
                              anchor           = "GGCGCGCCATTTAAATGTCGAC",
                              anchor_mm        = 2L,
                              bc_mm_pass       = 2L,
                              shrna_mm_pass    = 4L,
                              min_read_len     = 150L,
                              abi_zip          = NULL,
                              phd_prefix       = NULL,
                              sh_id_prefix     = "siKDBC_") {

  if (!requireNamespace("readxl", quietly = TRUE))
    stop("Package 'readxl' is required. Install with: install.packages('readxl')",
         call. = FALSE)

  ptm <- proc.time()

  # ── Validate inputs ────────────────────────────────────────────────────────
  if (!file.exists(fasta))
    stop("FASTA file not found: ", fasta, call. = FALSE)
  if (!file.exists(design_xlsx))
    stop("Design Excel not found: ", design_xlsx, call. = FALSE)
  if (!is.null(abi_zip) && !file.exists(abi_zip))
    stop("ABI ZIP not found: ", abi_zip, call. = FALSE)

  dirs <- .setup_output_dirs(output_dir, "sangerQC", sample_name)

  message("=============================================================")
  message("  catcheR::catcheR_sangerQC()")
  message("=============================================================")
  message(sprintf("  FASTA       : %s", fasta))
  message(sprintf("  Design      : %s  [sheet: %s]", design_xlsx, plate_sheet))
  message(sprintf("  Anchor      : %s  (max %d mm)", anchor, anchor_mm))
  message(sprintf("  BC pass     : <= %d mm", bc_mm_pass))
  message(sprintf("  shRNA pass  : <= %d total mm", shrna_mm_pass))
  if (!is.null(abi_zip)) message(sprintf("  ABI ZIP     : %s", abi_zip))

  # ── Internal helpers ────────────────────────────────────────────────────────

  .rc_dna <- function(seq) {
    comp <- c(A="T", T="A", G="C", C="G", N="N")
    paste(rev(comp[strsplit(toupper(seq), "")[[1]]]), collapse = "")
  }

  .hamming <- function(a, b, n = NULL) {
    if (is.null(n)) n <- min(nchar(a), nchar(b))
    a_ch <- strsplit(substr(a, 1, n), "")[[1]]
    b_ch <- strsplit(substr(b, 1, n), "")[[1]]
    l <- min(length(a_ch), length(b_ch))
    if (l == 0) return(999L)
    sum(a_ch[seq_len(l)] != b_ch[seq_len(l)])
  }

  .hamming_sense <- function(found, design_sense, n = 21L) {
    mm1 <- .hamming(found, design_sense, n)
    if (nchar(design_sense) > n) {
      mm2 <- .hamming(found, substr(design_sense, 2L, nchar(design_sense)), n)
      return(min(mm1, mm2))
    }
    mm1
  }

  .find_fuzzy <- function(seq, pattern, max_mm = 2L) {
    plen <- nchar(pattern)
    slen <- nchar(seq)
    if (slen < plen) return(list(pos = -1L, mm = 999L))
    best_pos <- -1L; best_mm <- max_mm + 1L
    pat_ch <- strsplit(pattern, "")[[1]]
    for (i in seq_len(slen - plen + 1L)) {
      sub_ch <- strsplit(substr(seq, i, i + plen - 1L), "")[[1]]
      mm <- sum(pat_ch != sub_ch)
      if (mm < best_mm) { best_mm <- mm; best_pos <- i }
      if (best_mm == 0L) break
    }
    if (best_mm <= max_mm) list(pos = best_pos, mm = best_mm)
    else                   list(pos = -1L,       mm = 999L)
  }

  .mut_list <- function(found, expected, prefix = "") {
    if (nchar(found) == 0 || nchar(expected) == 0) return("")
    n   <- min(nchar(found), nchar(expected))
    fch <- strsplit(substr(found,    1, n), "")[[1]]
    ech <- strsplit(substr(expected, 1, n), "")[[1]]
    idx <- which(fch != ech)
    if (length(idx) == 0) return("")
    paste(sprintf("%sp%d:%s\u2192%s", prefix, idx, ech[idx], fch[idx]),
          collapse = ", ")
  }

  # ── Parse FASTA ─────────────────────────────────────────────────────────────
  message("\n  Parsing FASTA...")
  lines      <- readLines(fasta)
  sequences  <- list()
  cur_id     <- NULL
  cur_seq    <- character(0)

  for (ln in lines) {
    ln <- trimws(ln)
    if (nchar(ln) == 0) next
    if (startsWith(ln, ">")) {
      if (!is.null(cur_id))
        sequences[[cur_id]] <- toupper(paste(cur_seq, collapse = ""))
      # Well ID = last underscore-delimited token
      parts  <- strsplit(sub("^>", "", ln), "_")[[1]]
      cur_id <- parts[length(parts)]
      cur_seq <- character(0)
    } else {
      cur_seq <- c(cur_seq, ln)
    }
  }
  if (!is.null(cur_id))
    sequences[[cur_id]] <- toupper(paste(cur_seq, collapse = ""))

  message(sprintf("  %d sequences loaded", length(sequences)))

  # ── Parse design plate ──────────────────────────────────────────────────────
  message("  Parsing design Excel...")
  plate_raw <- readxl::read_excel(design_xlsx, sheet = plate_sheet,
                                   col_names = FALSE)

  # Rows first_row..last_row (1-based in R, so +1 for header offset)
  row_range  <- seq(first_row + 1L, last_row + 1L)
  oligo_rows <- plate_raw[row_range, , drop = FALSE]

  design <- data.frame(
    oligo_well  = character(0),
    sh_id       = character(0),
    bc_design   = character(0),
    bc_in_read  = character(0),
    sense       = character(0),
    anti        = character(0),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(oligo_rows))) {
    row  <- oligo_rows[i, ]
    well <- trimws(as.character(row[[1]]))   # col index 2 = position column
    nm   <- trimws(as.character(row[[oligo_col + 1L]])) #oligo name
    top  <- trimws(toupper(as.character(row[[sequence_col + 1L]]))) #sequence

    if (is.na(well) || well == "NA" || nchar(well) == 0) next
    if (is.na(nm)   || nm   == "NA" || nchar(nm)   == 0) next
    if (is.na(top)  || top  == "NA" || nchar(top)  == 0) next

    # BC: between GGCGCGCC+NNNNNN and ACGCGT in TOP OLIGO
    bc_m <- regmatches(top, regexpr("GGCGCGCC[ACGTN]{6}([ACGT]{8})ACGCGT", top))
    if (length(bc_m) == 0) next
    bc_design  <- substr(bc_m, 15L, 22L)   # chars 15-22 after GGCGCGCC(8)+NNNNNN(6)
    bc_in_read <- .rc_dna(bc_design)

    # shRNA: between GGATCCC and TTTTTTT in TOP OLIGO
    shrna_m <- regmatches(top, regexpr("GGATCCC([ACGT]+?)TTTTTTT", top))
    if (length(shrna_m) == 0) next
    shrna_seq <- substr(shrna_m, 8L, nchar(shrna_m) - 7L)

    if (grepl("CTCGAG", shrna_seq)) {
      parts <- strsplit(shrna_seq, "CTCGAG")[[1]]
      sense <- parts[1]; anti <- if (length(parts) > 1) parts[2] else ""
    } else {
      sense <- shrna_seq; anti <- ""
    }

    sh_id <- sub(sh_id_prefix, "", nm, fixed = TRUE)

    design <- rbind(design, data.frame(
      oligo_well = well, sh_id = sh_id,
      bc_design  = bc_design, bc_in_read = bc_in_read,
      sense = sense, anti = anti,
      stringsAsFactors = FALSE
    ))
  }

  message(sprintf("  %d oligos loaded from plate sheet", nrow(design)))

  # ── Load ABI quality (optional) ─────────────────────────────────────────────
  .load_quals <- function(well) {
    if (is.null(abi_zip)) return(integer(0))
    zcon <- tryCatch(unz(abi_zip, ""), error = function(e) NULL)
    # List files to find the PHD for this well
    all_files <- tryCatch({
      con <- unz(abi_zip, "")  # dummy — use utils::unzip to list
      utils::unzip(abi_zip, list = TRUE)$Name
    }, error = function(e) character(0))

    # Find matching PHD
    phd_file <- all_files[grepl(paste0("_", well, "\\.phd\\.1$"), all_files)]
    if (length(phd_file) == 0) return(integer(0))

    lines <- tryCatch(
      readLines(unz(abi_zip, phd_file[1])),
      error = function(e) character(0)
    )
    if (length(lines) == 0) return(integer(0))

    in_dna <- FALSE; quals <- integer(0)
    for (ln in lines) {
      if (grepl("BEGIN_DNA", ln)) { in_dna <- TRUE; next }
      if (grepl("END_DNA",   ln)) break
      if (in_dna && nchar(trimws(ln)) > 0) {
        parts <- strsplit(trimws(ln), "\\s+")[[1]]
        if (length(parts) >= 2) {
          q <- suppressWarnings(as.integer(parts[2]))
          quals <- c(quals, if (is.na(q)) 0L else q)
        }
      }
    }
    quals
  }

  .region_quality <- function(quals, start, end) {
    if (length(quals) == 0 || start < 1 || end > length(quals) || start > end)
      return(list(min = NA_integer_, mean = NA_real_))
    qr <- quals[seq(start, end)]
    list(min = min(qr), mean = round(mean(qr), 1))
  }

  # ── Process each Sanger well ─────────────────────────────────────────────
  message("  Processing Sanger wells...")

  rows_order   <- c("A","B","C","D","E","F","G","H")
  sanger_wells <- paste0(rep(rows_order, times = 12),
                         rep(1:12, each  = 8))

  results <- vector("list", length(sanger_wells))

  for (wi in seq_along(sanger_wells)) {
    swell <- sanger_wells[wi]
    col_n <- as.integer(gsub("[A-Z]", "", swell))
    lig   <- if (col_n <= ligation_split) ligation_labels[1] else ligation_labels[2]
    seq   <- sequences[[swell]]
    slen  <- if (is.null(seq)) 0L else nchar(seq)

    r <- list(
      sanger_well         = swell,
      ligation            = lig,
      seq_len             = slen,
      status              = NA_character_,
      assigned_sh         = NA_character_,
      assigned_oligo_well = NA_character_,
      bc_found            = NA_character_,
      exp_bc_in_read      = NA_character_,
      bc_mm               = NA_integer_,
      bc_mutations        = NA_character_,
      uci                 = NA_character_,
      sense_found         = NA_character_,
      sense_mm            = NA_integer_,
      sense_mutations     = NA_character_,
      anti_found          = NA_character_,
      anti_mm             = NA_integer_,
      anti_mutations      = NA_character_,
      hairpin_ok          = NA,
      min_q_bc            = NA_integer_,
      mean_q_bc           = NA_real_,
      min_q_sense         = NA_integer_,
      mean_q_sense        = NA_real_,
      anchor_mm           = NA_integer_,
      notes               = NA_character_
    )

    if (is.null(seq) || slen < min_read_len) {
      r$status <- if (is.null(seq)) "SEQ FAILED" else
        sprintf("SEQ FAILED (<150 bp, got %d)", slen)
      results[[wi]] <- r; next
    }

    # Find anchor
    anc <- .find_fuzzy(seq, anchor, anchor_mm)
    if (anc$pos < 0) {
      r$status <- "SEQ TRUNCATED"
      r$notes  <- sprintf("Anchor '%s' not found (max %d mm)", anchor, anchor_mm)
      results[[wi]] <- r; next
    }

    anc_pos <- anc$pos  # 1-based start of anchor in seq
    r$anchor_mm <- anc$mm

    # BC and UCI extraction (left of anchor)
    bc_start  <- anc_pos - 14L   # 8 bp BC + 6 bp UCI before anchor
    uci_start <- anc_pos - 6L
    if (bc_start < 1L) {
      r$status <- "SEQ TRUNCATED"
      r$notes  <- "Read starts too close to anchor to extract BC/UCI"
      results[[wi]] <- r; next
    }

    bc_found  <- substr(seq, bc_start,  bc_start + 7L)
    uci_found <- substr(seq, uci_start, uci_start + 5L)
    mlu_seq   <- substr(seq, bc_start - 6L, bc_start - 1L)
    r$bc_found <- bc_found
    r$uci      <- uci_found

    # Quality for BC/UCI region
    quals <- .load_quals(swell)
    if (length(quals) > 0) {
      q_bc  <- .region_quality(quals, bc_start, anc_pos - 1L)
      r$min_q_bc  <- q_bc$min
      r$mean_q_bc <- q_bc$mean
    }

    # Sense and antisense extraction (right of anchor)
    sal_end <- anc_pos + nchar(anchor)  # anchor ends with GTCGAC
    xho_pos <- regexpr("CTCGAG", substr(seq, sal_end, nchar(seq)))
    if (xho_pos[1] < 0) {
      sense_found <- ""; anti_found <- ""
      r$hairpin_ok <- NA
    } else {
      xho_abs  <- sal_end + xho_pos[1] - 1L  # absolute pos of CTCGAG in seq (start)
      raw_after <- substr(seq, sal_end, xho_abs - 1L)
      sense_found <- sub("^A+", "", raw_after)   # strip leading A's (polyA)
      anti_end    <- xho_abs + 5L + nchar(sense_found)
      anti_found  <- substr(seq, xho_abs + 6L, anti_end)
      r$hairpin_ok <- (nchar(sense_found) > 0 && nchar(anti_found) > 0 &&
                       sense_found == .rc_dna(anti_found))

      if (length(quals) > 0 && xho_abs > sal_end) {
        q_s <- .region_quality(quals, sal_end, xho_abs - 1L)
        r$min_q_sense  <- q_s$min
        r$mean_q_sense <- q_s$mean
      }
    }
    r$sense_found <- if (nchar(sense_found) > 0) substr(sense_found, 1, 25) else NA
    r$anti_found  <- if (nchar(anti_found)  > 0) substr(anti_found,  1, 25) else NA

    # ── Match against all designs ──────────────────────────────────────────
    # Primary: shRNA match (more signal = 42 bp vs 8 bp)
    shrna_best_idx <- NA_integer_; shrna_best_mm <- 999L
    if (nchar(sense_found) >= 10) {
      for (i in seq_len(nrow(design))) {
        mm <- .hamming_sense(sense_found, design$sense[i]) +
      .hamming(anti_found, design$anti[i], 21L)
        if (mm < shrna_best_mm) { shrna_best_mm <- mm; shrna_best_idx <- i }
      }
    }

    # Secondary: BC match (8 bp, compare bc_found to bc_in_read)
    bc_best_idx <- NA_integer_; bc_best_mm <- 999L
    for (i in seq_len(nrow(design))) {
      mm <- .hamming(bc_found, design$bc_in_read[i], 8L)
      if (mm < bc_best_mm) { bc_best_mm <- mm; bc_best_idx <- i }
    }

    # Assign from shRNA match (more reliable)
    if (!is.na(shrna_best_idx)) {
      r$assigned_sh         <- design$sh_id[shrna_best_idx]
      r$assigned_oligo_well <- design$oligo_well[shrna_best_idx]
      r$exp_bc_in_read      <- design$bc_in_read[shrna_best_idx]
      r$bc_mm  <- .hamming(bc_found, design$bc_in_read[shrna_best_idx], 8L)
      r$sense_mm <- .hamming_sense(sense_found, design$sense[shrna_best_idx])
      r$anti_mm  <- .hamming(anti_found,  design$anti[shrna_best_idx],  21L)
      r$bc_mutations    <- .mut_list(bc_found, design$bc_in_read[shrna_best_idx])
      sense_for_mut <- design$sense[shrna_best_idx]
      if (nchar(sense_for_mut) > 21L &&
          .hamming(sense_found, substr(sense_for_mut, 2L, nchar(sense_for_mut)), 21L) <
          .hamming(sense_found, sense_for_mut, 21L)) {
        sense_for_mut <- substr(sense_for_mut, 2L, nchar(sense_for_mut))
      }
      r$sense_mutations <- .mut_list(sense_found, sense_for_mut)
      r$anti_mutations  <- .mut_list(anti_found,  design$anti[shrna_best_idx])
    }

    # Check BC/shRNA consistency
    bc_sh    <- if (!is.na(bc_best_idx))    design$sh_id[bc_best_idx]    else NA
    shrna_sh <- if (!is.na(shrna_best_idx)) design$sh_id[shrna_best_idx] else NA
    consistent <- identical(bc_sh, shrna_sh)

    # ── Status ────────────────────────────────────────────────────────────
    sm   <- if (is.na(r$sense_mm)) 999L else r$sense_mm
    am   <- if (is.na(r$anti_mm))  999L else r$anti_mm
    bcm  <- if (is.na(r$bc_mm))    999L else r$bc_mm
    hp   <- isTRUE(r$hairpin_ok)
    has_shrna <- !is.na(shrna_best_idx) && shrna_best_mm < 999L

    r$status <-
      if (!isTRUE(r$hairpin_ok) && !is.na(r$hairpin_ok) && has_shrna) {
        "INVALID HAIRPIN"
      } else if (!has_shrna) {
        "NOT IN LIBRARY"
      } else if (sm == 0 && am == 0 && bcm == 0) {
        "PASS"
      } else if (sm == 0 && am == 0 && bcm == 1) {
        "PASS (1bp BC mut)"
      } else if (sm == 0 && am == 0 && bcm == 2) {
        "PASS (2bp BC mut)"
      } else if (sm <= 1 && am <= 1 && bcm == 0) {
        "PASS (1bp shRNA mut)"
      } else if (sm <= 2 && am <= 2 && bcm == 0) {
        "PASS (2bp shRNA mut)"
      } else if (sm == 0 && am == 0 && bcm > bc_mm_pass) {
        "shRNA OK / BC WRONG"
      } else if (!consistent && bcm <= bc_mm_pass && (sm + am) <= shrna_mm_pass) {
        "BC/shRNA discordant"
      } else if (bcm == 0 && (sm + am) > shrna_mm_pass) {
        "CHECK"
      } else {
        "CHECK"
      }

    # Notes
    note_parts <- character(0)
    if (!identical(mlu_seq, "ACGCGT"))
      note_parts <- c(note_parts, sprintf("MluI mutated: %s", mlu_seq))
    if (!consistent && !is.na(bc_sh) && !is.na(shrna_sh))
      note_parts <- c(note_parts,
                      sprintf("BC→%s shRNA→%s", bc_sh, shrna_sh))
    r$notes <- if (length(note_parts)) paste(note_parts, collapse = "; ") else NA

    results[[wi]] <- r
  }

  # ── Assemble data frame ──────────────────────────────────────────────────
  res_df <- do.call(rbind, lapply(results, as.data.frame,
                                   stringsAsFactors = FALSE))
  rownames(res_df) <- NULL

  # ── Summary ──────────────────────────────────────────────────────────────
  summary_counts <- sort(table(res_df$status), decreasing = TRUE)

  pass_12  <- sum(grepl("^PASS", res_df$status[res_df$ligation == ligation_labels[1]]))
  pass_20  <- sum(grepl("^PASS", res_df$status[res_df$ligation == ligation_labels[2]]))
  total_12 <- sum(res_df$ligation == ligation_labels[1])
  total_20 <- sum(res_df$ligation == ligation_labels[2])

  # ── Save outputs ──────────────────────────────────────────────────────────
  .save_csv(res_df, dirs, "sanger_QC_results")

  stats_lines <- c(
    "=============================================================",
    "  SANGER QC SUMMARY",
    "=============================================================",
    sprintf("  FASTA    : %s", fasta),
    sprintf("  Design   : %s  [%s]", design_xlsx, plate_sheet),
    sprintf("  Oligos   : %d", nrow(design)),
    sprintf("  Wells    : %d", nrow(res_df)),
    "",
    "STATUS COUNTS",
    "",
    sprintf("  %-28s %s", "Status", "Count"),
    strrep("-", 40),
    unname(mapply(function(s, n) sprintf("  %-28s %d", s, n),
                  names(summary_counts), as.integer(summary_counts))),
    "",
    "BY LIGATION",
    sprintf("  %s: %d / %d pass (%.0f%%)", ligation_labels[1],
            pass_12, total_12, 100 * pass_12 / max(total_12, 1)),
    sprintf("  %s: %d / %d pass (%.0f%%)", ligation_labels[2],
            pass_20, total_20, 100 * pass_20 / max(total_20, 1)),
    "",
    "EXTRACTION PARAMETERS",
    sprintf("  Anchor           : %s  (max %d mm)", anchor, anchor_mm),
    sprintf("  BC position      : 8 bp immediately before anchor"),
    sprintf("  UCI position     : 6 bp before BC"),
    sprintf("  BC orientation   : RC of design BC column"),
    sprintf("  Sense            : between anchor end + polyA and CTCGAG"),
    sprintf("  Antisense        : 21-22 bp after CTCGAG"),
    "============================================================="
  )
  .save_stats(stats_lines, dirs, "sanger_QC_summary")

  result <- list(
    results  = res_df,
    design   = design,
    summary  = summary_counts,
    paths    = list(
      output     = dirs$root,
      csv        = dirs$csv,
      stats      = dirs$stats,
      R_objects  = dirs$R_objects,
      to_scratch = dirs$to_scratch
    )
  )

  .save_r_objects(result, dirs, "sangerQC")
  .to_scratch(file.path(dirs$csv, "sanger_QC_results.csv"), dirs)

  elapsed <- proc.time() - ptm
  message("")
  message("=== SUMMARY ===")
  for (i in seq_along(summary_counts))
    message(sprintf("  %-28s %d", names(summary_counts)[i],
                    as.integer(summary_counts)[i]))
  message(sprintf("\n  %s: %d/%d pass  |  %s: %d/%d pass",
                  ligation_labels[1], pass_12, total_12,
                  ligation_labels[2], pass_20, total_20))
  message(sprintf("\n  Done in %.1f seconds. Output in:", elapsed["elapsed"]))
  message(sprintf("    %s", dirs$root))
  message("=============================================================")

  invisible(result)
}
