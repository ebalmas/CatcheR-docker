test_that("birthday math is correct", {
  # E[distinct] with n=10, K=10 should be close to 10*(1-(9/10)^10) = 6.51
  expect_equal(round(catcheR:::.expected_distinct(10, 10), 2), 6.51)

  # With K=1, % dup should be 0
  expect_equal(catcheR:::.pct_dup(5, 1), 0)

  # With n > K, all cells should be unique
  expect_lt(catcheR:::.pct_dup(1000, 5), 1)

  # P(at least 1 distinct from 1 UCI in 1 draw) = 1
  expect_equal(catcheR:::.p_at_least_m(1, 1, 1), 1)

  # P(at least 2 distinct from 1 UCI) = 0 (impossible)
  expect_equal(catcheR:::.p_at_least_m(1, 10, 2), 0)
})

test_that("clone_stats returns correct structure", {
  # Create minimal test data
  df <- data.frame(
    name = c("GENE_1", "GENE_1", "UNMATCHED"),
    Freq = c(500, 100, 50),
    stringsAsFactors = FALSE
  )
  s <- catcheR:::.clone_stats(df, DIs = 300)
  expect_equal(s$total,         3)
  expect_equal(s$matched,       2)
  expect_equal(s$above_DIs,     1)
  expect_equal(s$above_matched, 1)
})

test_that("fmt_p formats p-values correctly", {
  expect_equal(catcheR:::.fmt_p(0.00001), "P < 0.0001")
  expect_match(catcheR:::.fmt_p(0.05), "P = 0.0500")
})
