test_that("reverse complement is an involution and correct", {
  expect_equal(catcheR:::.cat_rc("GGCGCGTTCATCTGGGGGAGCCG"), "CGGCTCCCCCAGATGAACGCGCC")
  x <- c("ACGT", "TTTT", "GATTACA")
  expect_equal(catcheR:::.cat_rc(catcheR:::.cat_rc(x)), x)
})

test_that("PIPseq R1 parser finds the 4-tier barcode at every phasing offset", {
  bc <- list(b1 = "AAAACCCC", b2 = "GGGTTT", b3 = "ACACAC", b4 = "TTTTGGGG")
  tail <- paste0("ACG", "ACG", strrep("T", 12))   # 3 nt index + ACG + oligo-dT
  build <- function(phase)
    paste0(substr("TCA", 1, phase), bc$b1, "ATG", bc$b2, "GAG", bc$b3, "TCGAG", bc$b4, tail)
  reads <- vapply(0:3, build, character(1))
  pr <- catcheR:::.pip_parse_r1(reads)
  expect_true(all(pr$valid))
  expect_equal(unique(pr$cell), paste0(bc$b1, bc$b2, bc$b3, bc$b4))
  expect_equal(unique(pr$index), "ACG")
})

test_that("linker-only matches are rejected without the downstream oligo structure", {
  bad <- paste0("AAAACCCC", "ATG", "GGGTTT", "GAG", "ACACAC", "TCGAG", "TTTTGGGG",
                "CCCCCCCCCCCCCCCCCCCC")           # linkers present, no index/ACG/dT
  pr <- catcheR:::.pip_parse_r1(bad)
  expect_true(pr$linker_only)
  expect_false(pr$valid)
})

test_that("the empty reference is classified before the downstream constant", {
  ANCH <- substr(catcheR:::.cat_rc("GGCGCGTTCATCTGGGGGAGCCG"), 1, 15)
  RC_EMPTY <- catcheR:::.cat_rc("TACGCGTTCATCTGGGGGAGCCG")
  CONST <- "ACGCGTAGCTCGCTGATCAG"
  empty <- paste0(RC_EMPTY, CONST)                 # junction regenerates CONST early
  pos <- regexpr(CONST, empty, fixed = TRUE)
  expect_equal(catcheR:::.cat_classify(empty, pos, ANCH, RC_EMPTY),
               "empty vector (no insertion)")
})

test_that("ValidCell follows catcheR's rule: single true UCI, then ratio", {
  tab <- data.frame(
    cellID            = c("a", "b", "b", "c", "c"),
    UMIxUCI           = c(5L, 4L, 1L, 3L, 3L),
    UMIxcellID        = c(5L, 5L, 5L, 6L, 6L),
    UMIpercentagexUCI = c(100, 80, 20, 50, 50))
  v <- catcheR:::.cat_validate(tab, UMI.count = 0, percentage = 15, ratio = 2)$cells
  # a: one UCI at 100%           -> valid via the pct == 100 escape clause
  # b: two UCIs both > 15%       -> True_UCIxcellID == 2, rejected
  # c: tied UCIs                 -> True == 2, rejected (and `second` must not be -Inf)
  expect_equal(v$ValidCell[v$cellID == "a"], TRUE)
  expect_equal(v$ValidCell[v$cellID == "b"], FALSE)
  expect_equal(v$ValidCell[v$cellID == "c"], FALSE)
  expect_equal(v$second[v$cellID == "c"], 3)
  # raising the threshold drops b's 1-read UCI and resolves the cell
  v1 <- catcheR:::.cat_validate(tab, UMI.count = 1, percentage = 15, ratio = 2)$cells
  expect_equal(v1$ValidCell[v1$cellID == "b"], TRUE)
})

# ── Parse (Evercode) ────────────────────────────────────────────────────────

test_that(".par_well_index follows Parse's 12-column plate numbering", {
  expect_equal(catcheR:::.par_well_index("A1"), 1)
  expect_equal(catcheR:::.par_well_index("A12"), 12)
  expect_equal(catcheR:::.par_well_index("B1"), 13)
  expect_equal(catcheR:::.par_well_index("E11"), 59)
  expect_equal(catcheR:::.par_well_index("H12"), 96)
})

test_that(".par_whitelist reads the split-pipe format and keeps both round-1 stypes", {
  f <- tempfile(fileext = ".csv")
  writeLines(c(
    "# barcode data; Kit TEST; Chemistry v4",
    "# split-pipe v1.9.0",
    "# Round 1 test",
    "bci,sequence,uid,well,plate,stype,wind,bc_round",
    "1,AAAAAAAA,p1,A1,1,R,1,1",
    "2,CCCCCCCC,p2,A2,1,R,2,1",
    "97,GGGGGGGG,p3,A1,1,T,1,1",
    "98,TTTTTTTT,p4,A2,1,T,2,1",
    "# Round 2 test",
    "1,ACACACAC,p5,A1,1,L,1,2",
    "59,GTGTTCTA,p6,E11,1,L,59,2",
    "# Round 3 test",
    "1,TGTGTGTG,p7,A1,1,L,1,3",
    "59,CACACACA,p8,E11,1,L,59,3"), f)
  wl <- catcheR:::.par_whitelist(f)
  expect_equal(nrow(wl$table), 8L)
  expect_named(wl$maps, c("round1", "round2", "round3"))
  # the dT and random-hexamer sequences of a well map to the SAME well index
  expect_equal(unname(wl$maps$round1[["AAAAAAAA"]]), unname(wl$maps$round1[["GGGGGGGG"]]))
  expect_equal(length(wl$maps$round1), 4L)        # 2 wells x 2 stypes
  expect_equal(length(unique(wl$maps$round1)), 2L)
  expect_equal(unname(wl$maps$round2[["GTGTTCTA"]]), 59)
})

test_that(".par_corr_map adds unambiguous neighbours only", {
  m <- c(AAAAAAAA = 1L, TTTTTTTT = 2L)
  c1 <- catcheR:::.par_corr_map(m, max.edit.distance = 1L)
  expect_equal(length(c1), 2L + 2L * 24L)          # exact + 8 positions x 3 bases each
  expect_equal(unname(c1[["AAAAAAAC"]]), 1L)
  expect_equal(unname(catcheR:::.par_corr_map(m, max.edit.distance = 0L)), unname(m))
  # a variant reachable from two different wells is dropped, not assigned
  amb <- c(AAAAAAAA = 1L, AAAAAAAT = 2L)
  c2 <- catcheR:::.par_corr_map(amb, max.edit.distance = 1L)
  expect_false("AAAAAAAG" %in% names(c2))
})

test_that(".par_adapter emits split-pipe bc_wells IDs and reads the UMI", {
  wl <- list(round1 = c(AAAAAAAA = 1L), round2 = c(CCCCCCCC = 8L), round3 = c(GGGGGGGG = 59L))
  LA <- "ATGAGGGGTCAG"; LB <- "TCCAACCACCTC"
  ad <- catcheR:::.par_adapter(wl, "s1", LA, LB, umi.length = 10L, barcode.length = 8L)
  expect_equal(ad$cell_read, 2L)
  expect_true(ad$has_umi)
  good <- paste0("ACGTACGTAC", "GGGGGGGG", LA, "CCCCCCCC", LB, "AAAAAAAA", "TTTTTTTTTT")
  bad  <- paste0("ACGTACGTAC", "GGGGGGGG", "ATGAGGGGTCAA", "CCCCCCCC", LB, "AAAAAAAA", "TTTT")
  r <- ad$cell_fun(c(good, bad))
  expect_equal(r$valid, c(TRUE, FALSE))
  expect_equal(r$cell[1], "01_08_59__s1")
  expect_equal(r$umi[1], "ACGTACGTAC")
})

test_that(".cat_classify tolerates a leading base on the cassette read", {
  REF <- "GGCGCGTTCATCTGGGGGAGCCG"; EMPTY <- "TACGCGTTCATCTGGGGGAGCCG"
  RCREF <- catcheR:::.cat_rc(REF); RC_EMPTY <- catcheR:::.cat_rc(EMPTY)
  ANCH <- substr(RCREF, 1, nchar(RCREF) - 8)
  CONST <- "ACGCGTAGCTCGCTGATCAG"
  read <- paste0("N", RCREF, "ACGTAC", "GGGGTTTT", CONST)
  pos  <- regexpr(CONST, read, fixed = TRUE)
  expect_equal(catcheR:::.cat_classify(read, pos, ANCH, RC_EMPTY, lead = 1L), "cassette")
  emptyread <- paste0("N", RC_EMPTY, "TTTTTTTTTT")
  expect_equal(catcheR:::.cat_classify(emptyread, -1L, ANCH, RC_EMPTY, lead = 1L),
               "empty vector (no insertion)")
})
