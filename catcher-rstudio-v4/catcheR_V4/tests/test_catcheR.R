# =============================================================================
# test_catcheR.R
# =============================================================================
# Run this script to verify the catcheR package installs correctly and
# all functions behave as expected without needing real data.
#
# USAGE:
#   Rscript test_catcheR.R           # from Terminal
#   source("test_catcheR.R")         # from RStudio
#
# REQUIREMENTS:
#   devtools::install("path/to/catcheR/")  # install the package first
#
# All tests use synthetic in-memory data — no FASTQ or ABI files needed.
# =============================================================================

cat("=============================================================\n")
cat("  catcheR package test suite\n")
cat("=============================================================\n\n")

# ── Helper ────────────────────────────────────────────────────────────────────
.pass <- function(msg) cat(sprintf("  [PASS] %s\n", msg))
.fail <- function(msg) { cat(sprintf("  [FAIL] %s\n", msg)); .n_fail <<- .n_fail + 1L }
.section <- function(msg) cat(sprintf("\n--- %s ---\n", msg))
.n_fail  <- 0L

# ── 1. Package loads ──────────────────────────────────────────────────────────
.section("1. Package installation")
tryCatch({
  library(catcheR)
  .pass("library(catcheR) loaded without errors")
}, error = function(e) {
  .fail(paste("library(catcheR) failed:", conditionMessage(e)))
  stop("Cannot continue without the package.", call. = FALSE)
})

# ── 2. All functions exported ─────────────────────────────────────────────────
.section("2. Exported functions")
expected_fns <- c(
  "catcheR_step1QC",
  "catcheR_step2QC_extraction",
  "catcheR_step2QC_plasmidQC",
  "catcheR_step2QC_combinedQC",
  "catcheR_step2QC_plots",
  "catcheR_10Xcatch",
  "catcheR_10XcatchQC",
  "catcheR_filtercatch",
  "catcheR_nocatch",
  "catcheR_sortcatch",
  "catcheR_sangerQC"
)
for (fn in expected_fns) {
  if (existsMethod <- tryCatch(is.function(get(fn, envir = asNamespace("catcheR"))),
                               error = function(e) FALSE)) {
    .pass(sprintf("%-35s exported", fn))
  } else {
    .fail(sprintf("%-35s NOT found in package namespace", fn))
  }
}

# ── 3. Bundled scripts present ────────────────────────────────────────────────
.section("3. Bundled scripts")
scripts <- c(
  "python/plasmid_final_corrected.py",
  "scripts/plasmid_inter.sh",
  "scripts/plasmid_inter2.R",
  "scripts/barcode_silencing_slicing.sh",
  "scripts/barcode_silencing_explorative_analysis.R",
  "scripts/barcode_silencing_cell_filtering.R",
  "scripts/barcode_silencing_empty_selection.R",
  "scripts/barcode_silencing_all_samples.R",
  "scripts/grep.sh"
)
for (s in scripts) {
  path <- system.file(s, package = "catcheR")
  if (nzchar(path) && file.exists(path)) {
    .pass(sprintf("%-50s found", s))
  } else {
    .fail(sprintf("%-50s MISSING from inst/", s))
  }
}

# ── 4. Internal helpers ───────────────────────────────────────────────────────
.section("4. Internal birthday-problem mathematics")
ns <- asNamespace("catcheR")

# .expected_distinct
ed <- get(".expected_distinct", envir = ns)(n = 10, K = 10)
if (abs(ed - 6.513) < 0.01) {
  .pass(".expected_distinct(10, 10) ≈ 6.51")
} else {
  .fail(sprintf(".expected_distinct(10,10) = %.4f, expected ≈ 6.51", ed))
}

# .pct_dup: K=1 → 0%
pd1 <- get(".pct_dup", envir = ns)(n = 5, K = 1)
if (pd1 == 0) {
  .pass(".pct_dup(n=5, K=1) = 0  (no duplicates possible with 1 clone)")
} else {
  .fail(sprintf(".pct_dup(n=5,K=1) = %.4f, expected 0", pd1))
}

# .pct_dup: large n, small K → near 0
pd2 <- get(".pct_dup", envir = ns)(n = 1000, K = 5)
if (pd2 < 1) {
  .pass(sprintf(".pct_dup(n=1000, K=5) = %.4f  (< 1%%, as expected)", pd2))
} else {
  .fail(sprintf(".pct_dup(n=1000,K=5) = %.4f, expected < 1%%", pd2))
}

# .p_at_least_m: n=1, K=1, m=1 → P=1
p1 <- get(".p_at_least_m", envir = ns)(n = 1, K = 1, m = 1)
if (abs(p1 - 1) < 1e-9) {
  .pass(".p_at_least_m(n=1, K=1, m=1) = 1  (certain)")
} else {
  .fail(sprintf(".p_at_least_m(1,1,1) = %.6f, expected 1", p1))
}

# .p_at_least_m: impossible case
p2 <- get(".p_at_least_m", envir = ns)(n = 1, K = 10, m = 2)
if (p2 == 0) {
  .pass(".p_at_least_m(n=1, K=10, m=2) = 0  (impossible, only 1 UCI)")
} else {
  .fail(sprintf(".p_at_least_m(1,10,2) = %.6f, expected 0", p2))
}

# .fmt_p
fmt_lo <- get(".fmt_p", envir = ns)(0.000001)
fmt_hi <- get(".fmt_p", envir = ns)(0.05)
if (fmt_lo == "P < 0.0001" && grepl("P = 0.0500", fmt_hi)) {
  .pass(".fmt_p() formats p-values correctly")
} else {
  .fail(sprintf(".fmt_p: got '%s' and '%s'", fmt_lo, fmt_hi))
}

# ── 5. .setup_output_dirs ─────────────────────────────────────────────────────
.section("5. Output directory setup")
tmp_base <- file.path(tempdir(), "catcheR_test_output")
dirs <- get(".setup_output_dirs", envir = ns)(tmp_base, "testfunc", "mysample")

expected_subdirs <- c("root","csv","plots","stats","R_objects","to_scratch")
for (sub in expected_subdirs) {
  if (!is.null(dirs[[sub]]) && dir.exists(dirs[[sub]])) {
    .pass(sprintf(".setup_output_dirs() created: %s/", sub))
  } else {
    .fail(sprintf(".setup_output_dirs() missing: %s/", sub))
  }
}
# Check naming convention: YYMMDD_testfunc_mysample
run_name <- basename(dirs$root)
if (grepl("^\\d{6}_testfunc_mysample$", run_name)) {
  .pass(sprintf("Run folder named correctly: %s", run_name))
} else {
  .fail(sprintf("Run folder name wrong: %s (expected YYMMDD_testfunc_mysample)", run_name))
}

# ── 6. .save_csv / .save_stats / .save_r_objects ─────────────────────────────
.section("6. File saving helpers")
test_df <- data.frame(x = 1:5, y = letters[1:5], stringsAsFactors = FALSE)
get(".save_csv", envir = ns)(test_df, dirs, "test_table")
csv_path <- file.path(dirs$csv, "test_table.csv")
if (file.exists(csv_path)) {
  df_back <- read.csv(csv_path, stringsAsFactors = FALSE)
  if (nrow(df_back) == 5 && ncol(df_back) == 2) {
    .pass(".save_csv() writes readable CSV with correct dimensions")
  } else {
    .fail(sprintf(".save_csv() CSV has wrong dims: %d×%d", nrow(df_back), ncol(df_back)))
  }
} else {
  .fail(".save_csv() did not create the file")
}

get(".save_stats", envir = ns)(c("Line 1", "Line 2", "Line 3"), dirs, "test_report")
stats_path <- file.path(dirs$stats, "test_report.txt")
if (file.exists(stats_path) && length(readLines(stats_path)) == 3) {
  .pass(".save_stats() writes correct number of lines")
} else {
  .fail(".save_stats() file missing or wrong line count")
}

test_result <- list(data = test_df, value = 42)
get(".save_r_objects", envir = ns)(test_result, dirs, "testfunc")
rds_path <- file.path(dirs$R_objects, "testfunc_result.rds")
if (file.exists(rds_path)) {
  loaded <- readRDS(rds_path)
  if (identical(loaded$value, 42L) || identical(loaded$value, 42)) {
    .pass(".save_r_objects() writes and reloads correctly")
  } else {
    .fail(".save_r_objects() reloaded value mismatch")
  }
} else {
  .fail(".save_r_objects() did not create .rds file")
}

# ── 7. .to_scratch ────────────────────────────────────────────────────────────
.section("7. to_scratch handoff")
dummy_file <- file.path(dirs$csv, "handoff_test.csv")
write.csv(data.frame(a=1), dummy_file, row.names=FALSE)
get(".to_scratch", envir = ns)(dummy_file, dirs)
scratch_copy <- file.path(dirs$to_scratch, "handoff_test.csv")
if (file.exists(scratch_copy)) {
  .pass(".to_scratch() correctly copies file to to_scratch/")
} else {
  .fail(".to_scratch() did not copy file")
}

# ── 8. .clone_stats ──────────────────────────────────────────────────────────
.section("8. Clone statistics (.clone_stats)")
clones_df <- data.frame(
  name = c(rep("GENE_1", 5), rep("UNMATCHED", 3), rep("GENE_2", 2)),
  Freq = c(500, 400, 300, 200, 50, 100, 20, 5, 600, 800),
  stringsAsFactors = FALSE
)
cs <- get(".clone_stats", envir = ns)(clones_df, DIs = 300)

if (cs$total == 10)         .pass(".clone_stats: total = 10")
  else .fail(sprintf(".clone_stats: total = %d, expected 10", cs$total))
if (cs$matched == 7)        .pass(".clone_stats: matched = 7")
  else .fail(sprintf(".clone_stats: matched = %d, expected 7", cs$matched))
if (cs$above_DIs == 4)      .pass(".clone_stats: above_DIs = 4  (Freq >= 300)")
  else .fail(sprintf(".clone_stats: above_DIs = %d, expected 4", cs$above_DIs))
if (cs$above_matched == 3)  .pass(".clone_stats: above_matched = 3")
  else .fail(sprintf(".clone_stats: above_matched = %d, expected 3", cs$above_matched))

# ── 9. .birthday_results ─────────────────────────────────────────────────────
.section("9. Birthday analysis (.birthday_results)")
shrna_df <- data.frame(
  name    = paste0("shRNA_", 1:5),
  gene    = rep("GENE1", 5),
  n_UCI   = c(8L, 10L, 6L, 12L, 9L),
  n_reads = c(1000L, 1200L, 800L, 1500L, 950L),
  pct_UCI   = rep(20, 5),
  pct_reads = rep(20, 5),
  stringsAsFactors = FALSE
)
br <- get(".birthday_results", envir = ns)(
  shrna = shrna_df, N = 5L, DIs = 300L,
  eff = 100 / 2e6,
  nucleofect_cells = 2e6,
  transfect_clones = 100L,
  transfect_cells  = 2e6
)

if (is.data.frame(br$results) && nrow(br$results) == 10) {
  .pass(".birthday_results: returns 10-row results table (m = 1..10)")
} else {
  .fail(sprintf(".birthday_results: results has %d rows, expected 10",
                nrow(br$results)))
}
if (all(c("m_target","K","cells_M","pct_dup_mean","verdict") %in% names(br$results))) {
  .pass(".birthday_results: results has all expected columns")
} else {
  .fail(".birthday_results: missing columns in results")
}
if (br$results$pct_dup_mean[1] == 0) {
  .pass(".birthday_results: m=1 has 0% duplication (K=1, always unique)")
} else {
  .fail(sprintf(".birthday_results: m=1 pct_dup = %.1f, expected 0",
                br$results$pct_dup_mean[1]))
}
if (!is.null(br$K_planned)) {
  .pass(sprintf(".birthday_results: K_planned = %.4f", br$K_planned))
} else {
  .fail(".birthday_results: K_planned is NULL")
}
if ("dup_status" %in% names(br$shrna)) {
  .pass(".birthday_results: per-shRNA table has dup_status column")
} else {
  .fail(".birthday_results: per-shRNA table missing dup_status")
}

# ── 10. catcheR_sangerQC internal helpers ────────────────────────────────────
.section("10. catcheR_sangerQC internal DNA helpers")

# We test the internal helpers by sourcing the function and calling them
# through a minimal in-memory run
fasta_tmp  <- tempfile(fileext = ".fasta")
design_tmp <- tempfile(fileext = ".xlsx")

# Write a minimal synthetic FASTA (2 wells)
# Well A1 = COA5_5 sequence from the real experiment
# Read structure: ...ACGCGT + TTGAAGCC(BC_RC) + AAACAG(UCI_RC) + GGCGCGCCATTTAAATGTCGAC + AAAAAAA + GCAGTGTTTGAAGGAAGGATA(sense) + CTCGAG + TATCCTTCCTTCAAACACTGC(anti)...
seq_a1 <- paste0(
  "CTGATCAGCGAGCTACGCGT",          # vector left
  "TTGAAGCC",                       # BC_RC (RC of GGCTTCAA)
  "AAACAG",                         # UCI_RC
  "GGCGCGCCATTTAAATGTCGAC",         # anchor
  "AAAAAAAGCAGTGTTTGAAGGAAGGATA",   # polyA + sense
  "CTCGAG",                         # loop
  "TATCCTTCCTTCAAACACTGC",          # antisense
  "GGGATCTCTATCACTGATAGGG"          # vector right
)
# Well A2 = B2M_1-like (perfect match to design)
seq_a2 <- paste0(
  "CTGATCAGCGAGCTACGCGT",
  "GGTTCAGG",                       # BC_RC of CCTGAACC
  "AATTAG",                         # UCI_RC
  "GGCGCGCCATTTAAATGTCGAC",
  "AAAAAAAGCAGCAGAGAATGGAAAGTCAA",  # polyA + sense B2M
  "CTCGAG",
  "TTGACTTTCCATTCTCTGCTG",          # antisense B2M
  "GGGATCTCTATCACTGATAGGG"
)

fasta_lines <- c(
  ">KYO_123456_A1", seq_a1,
  ">KYO_123456_A2", seq_a2
)
writeLines(fasta_lines, fasta_tmp)

# Write a minimal synthetic design Excel using openxlsx if available,
# otherwise skip the full integration test
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Plate 01")
  # Row 1 = header row (rows 1-7 skipped by first_row=8 → R index 9)
  # Row 9 (0-based row 8) = first oligo
  # Col 2=position, col 3=oligo name, col 5=sequence
  header_row <- c(NA, "Position", "Oligo Name", "5' Mods", "Sequence")
  openxlsx::writeData(wb, "Plate 01", rbind(
    # rows 1-7: filler
    matrix("", nrow=7, ncol=5),
    # row 8: headers
    header_row,
    # row 9: A1 = COA5_5
    c(NA, "A1", "siKDBC_COA5_5", NA,
      paste0("AGACTCGGATCCCGCAGTGTTTGAAGGAAGGATACTCGAGTATCCTTCCTTCAAACACTGC",
             "TTTTTTTGTCGACATTTAAATGGCGCGCCNNNNNNGGCTTCAAACGCGTGTAGCTCGCTGATCAGC")),
    # row 10: B1 = B2M_1
    c(NA, "B1", "siKDBC_B2M_1", NA,
      paste0("AGACTCGGATCCCGCAGCAGAGAATGGAAAGTCAACTCGAGTTGACTTTCCATTCTCTGCTG",
             "TTTTTTTGTCGACATTTAAATGGCGCGCCNNNNNNCCTGAACCACGCGTGTAGCTCGCTGATCAGC"))
  ), colNames = FALSE)
  openxlsx::saveWorkbook(wb, design_tmp, overwrite = TRUE)

  # Run catcheR_sangerQC with the synthetic data
  out_dir <- file.path(tempdir(), "sangerQC_test_output")
  result <- tryCatch(
    catcheR_sangerQC(
      fasta        = fasta_tmp,
      design_xlsx  = design_tmp,
      output_dir   = out_dir,
      sample_name  = "synthetic",
      plate_sheet  = "Plate 01",
      first_row    = 8L,
      last_row     = 9L,
      ligation_split = 6L,
      abi_zip      = NULL
    ),
    error = function(e) {
      .fail(paste("catcheR_sangerQC() threw error:", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(result)) {
    .pass("catcheR_sangerQC() ran without error on synthetic data")

    if (is.data.frame(result$results)) {
      .pass(sprintf("catcheR_sangerQC() returned results data frame (%d rows)",
                    nrow(result$results)))
    } else {
      .fail("catcheR_sangerQC() result$results is not a data frame")
    }

    a1_row <- result$results[result$results$sanger_well == "A1", ]
    if (nrow(a1_row) == 1 && grepl("^PASS", a1_row$status)) {
      .pass(sprintf("A1 correctly identified as COA5_5 with status: %s", a1_row$status))
    } else if (nrow(a1_row) == 1) {
      .fail(sprintf("A1 status = '%s', expected PASS (COA5_5)", a1_row$status))
    } else {
      .fail("A1 not found in results")
    }

    a2_row <- result$results[result$results$sanger_well == "A2", ]
    if (nrow(a2_row) == 1 && grepl("^PASS", a2_row$status)) {
      .pass(sprintf("A2 correctly identified as B2M_1 with status: %s", a2_row$status))
    } else if (nrow(a2_row) == 1) {
      .fail(sprintf("A2 status = '%s', expected PASS (B2M_1)", a2_row$status))
    } else {
      .fail("A2 not found in results")
    }

    # Check output files were created
    csv_out <- file.path(result$paths$csv, "sanger_QC_results.csv")
    if (file.exists(csv_out)) {
      .pass("catcheR_sangerQC() wrote sanger_QC_results.csv")
    } else {
      .fail("catcheR_sangerQC() did not write sanger_QC_results.csv")
    }

    rds_out <- file.path(result$paths$R_objects, "sangerQC_result.rds")
    if (file.exists(rds_out)) {
      reloaded <- readRDS(rds_out)
      if (!is.null(reloaded$results)) {
        .pass("catcheR_sangerQC() .rds reloads correctly")
      } else {
        .fail("catcheR_sangerQC() reloaded .rds missing $results")
      }
    } else {
      .fail("catcheR_sangerQC() did not write .rds file")
    }
  }
} else {
  cat("  [SKIP] openxlsx not installed — skipping catcheR_sangerQC integration test\n")
  cat("         Install with: install.packages('openxlsx')\n")
}

# ── 11. Input validation (error handling) ─────────────────────────────────────
.section("11. Input validation")

# catcheR_step2QC_plasmidQC: missing file → informative error
err <- tryCatch(
  catcheR_step2QC_plasmidQC(results_dir = "/nonexistent/path/"),
  error = function(e) conditionMessage(e)
)
if (grepl("distribution_all_clones", err)) {
  .pass("catcheR_step2QC_plasmidQC: missing file gives informative error")
} else {
  .fail(sprintf("catcheR_step2QC_plasmidQC: unexpected error message: %s", err))
}

# catcheR_step1QC: missing folder → informative error
err2 <- tryCatch(
  catcheR_step1QC(folder = "/nonexistent/", fastq.read1 = "x.fastq"),
  error = function(e) conditionMessage(e)
)
if (grepl("not found", err2, ignore.case = TRUE)) {
  .pass("catcheR_step1QC: missing folder gives informative error")
} else {
  .fail(sprintf("catcheR_step1QC: unexpected error: %s", err2))
}

# catcheR_sangerQC: missing FASTA → informative error
err3 <- tryCatch(
  catcheR_sangerQC(fasta = "/no/such.fasta", design_xlsx = design_tmp),
  error = function(e) conditionMessage(e)
)
if (grepl("not found|FASTA", err3, ignore.case = TRUE)) {
  .pass("catcheR_sangerQC: missing FASTA gives informative error")
} else {
  .fail(sprintf("catcheR_sangerQC: unexpected error: %s", err3))
}

# ── Final summary ─────────────────────────────────────────────────────────────
cat("\n=============================================================\n")
if (.n_fail == 0L) {
  cat("  ALL TESTS PASSED\n")
} else {
  cat(sprintf("  %d TEST(S) FAILED — review output above\n", .n_fail))
}
cat("=============================================================\n")

# Cleanup temp output
unlink(file.path(tempdir(), "catcheR_test_output"), recursive = TRUE)
unlink(file.path(tempdir(), "catcheR_test_output"), recursive = TRUE)
unlink(fasta_tmp); unlink(design_tmp)

invisible(.n_fail)
