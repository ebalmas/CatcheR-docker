# =============================================================================
# R/utils.R
# Internal helper functions — not exported to the user namespace.
# Prefixed with . to signal they are private.
# =============================================================================

# ── Formatting ────────────────────────────────────────────────────────────────

.fmt_p <- function(p) {
  if (p < 0.0001) "P < 0.0001" else sprintf("P = %.4f", p)
}

# ── Structured output directory setup ────────────────────────────────────────
# Creates the date-stamped, themed folder structure for every function run.
# Returns a list of paths: $root, $csv, $plots, $stats, $R_objects, $to_scratch

.setup_output_dirs <- function(output_base, func_name, sample_name = NULL) {
  date_stamp <- format(Sys.Date(), "%y%m%d")
  run_name   <- if (!is.null(sample_name) && nzchar(sample_name))
    paste0(date_stamp, "_", func_name, "_", sample_name)
  else
    paste0(date_stamp, "_", func_name)

  root <- file.path(output_base, func_name, run_name)

  dirs <- list(
    root       = root,
    csv        = file.path(root, "csv"),
    plots      = file.path(root, "plots"),
    stats      = file.path(root, "stats"),
    R_objects  = file.path(root, "R_objects"),
    to_scratch = file.path(root, "to_scratch")
  )

  for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

  message(sprintf("  Output folder: %s", root))
  dirs
}

# ── Save a plot to the plots/ subfolder ───────────────────────────────────────
.save_plot <- function(p, dirs, name, w = 12, h = 7) {
  base <- file.path(dirs$plots, name)
  ggplot2::ggsave(paste0(base, ".pdf"), plot = p, width = w, height = h,
                  limitsize = FALSE)
  ggplot2::ggsave(paste0(base, ".jpg"), plot = p, width = w, height = h,
                  limitsize = FALSE)
  invisible(base)
}

# ── Save a data frame to the csv/ subfolder ───────────────────────────────────
.save_csv <- function(df, dirs, name) {
  path <- file.path(dirs$csv, paste0(name, ".csv"))
  utils::write.csv(df, path, row.names = FALSE)
  invisible(path)
}

# ── Save text report to the stats/ subfolder ─────────────────────────────────
.save_stats <- function(lines, dirs, name) {
  path <- file.path(dirs$stats, paste0(name, ".txt"))
  writeLines(lines, path)
  invisible(path)
}

# ── Auto-save all R objects and return the result list ───────────────────────
# Saves the full result list as a single .rds in R_objects/
# Individual named objects (data frames, plots) saved separately.
.save_r_objects <- function(result, dirs, func_name) {
  # Full result list
  saveRDS(result, file.path(dirs$R_objects,
                              paste0(func_name, "_result.rds")))
  # Individual tables
  for (nm in names(result)) {
    obj <- result[[nm]]
    if (is.data.frame(obj)) {
      saveRDS(obj, file.path(dirs$R_objects, paste0(nm, ".rds")))
    } else if (is.list(obj) && !inherits(obj, "gg")) {
      for (nm2 in names(obj)) {
        obj2 <- obj[[nm2]]
        if (is.data.frame(obj2) || inherits(obj2, "gg")) {
          saveRDS(obj2, file.path(dirs$R_objects,
                                   paste0(nm, "_", nm2, ".rds")))
        }
      }
    }
  }
  message(sprintf("  R objects saved to: %s/", basename(dirs$R_objects)))
  invisible(dirs$R_objects)
}

# ── Copy file(s) to the to_scratch/ folder ───────────────────────────────────
.to_scratch <- function(paths, dirs) {
  for (p in paths) {
    if (file.exists(p)) {
      file.copy(p, file.path(dirs$to_scratch, basename(p)),
                overwrite = TRUE)
    }
  }
  message(sprintf("  Intermediate files for next step -> %s/",
                  basename(dirs$to_scratch)))
  invisible(dirs$to_scratch)
}


.load_clones <- function(results_dir) {
  f <- file.path(results_dir, "distribution_all_clones.csv")
  if (!file.exists(f)) {
    stop(
      "Cannot find distribution_all_clones.csv in: ", results_dir,
      "\nHave you run catcheR_step2QC_extraction() first?",
      call. = FALSE
    )
  }
  read.csv(f, stringsAsFactors = FALSE)
}

.apply_colours <- function(p, results_dir) {
  col_file <- file.path(results_dir, "colors.csv")
  if (file.exists(col_file)) {
    cols <- read.csv(col_file, header = FALSE,
                     col.names = c("name", "colors"))
    p <- p + ggplot2::scale_fill_manual(
      values = stats::setNames(cols$colors, cols$name))
  }
  p
}

# ── Script path patching ──────────────────────────────────────────────────────
# The shell scripts bundled from the original Docker image contain hardcoded
# paths like /home/barcode_silencing_explorative_analysis.R. This function
# creates a temporary patched copy with the correct paths to the bundled R
# scripts, so the package works on any server without Docker.

.patch_sh_script <- function(sh_path, folder) {
  lines <- readLines(sh_path)

  scripts_dir <- system.file("scripts", package = "catcheR")
  grep_sh     <- file.path(scripts_dir, "grep.sh")

  # Map: old Docker path -> new bundled path
  replacements <- list(
    "/home/plasmid_inter2.R"                       =
      file.path(scripts_dir, "plasmid_inter2.R"),
    "/home/barcode_silencing_explorative_analysis.R" =
      file.path(scripts_dir, "barcode_silencing_explorative_analysis.R"),
    "/home/barcode_silencing_cell_filtering.R"     =
      file.path(scripts_dir, "barcode_silencing_cell_filtering.R"),
    "/home/barcode_silencing_empty_selection.R"    =
      file.path(scripts_dir, "barcode_silencing_empty_selection.R"),
    "/home/barcode_silencing_all_samples.R"        =
      file.path(scripts_dir, "barcode_silencing_all_samples.R"),
    "/home/grep.sh"                                =
      grep_sh,
    # Also patch chmod /data/scratch lines — not relevant outside Docker
    "chmod 777 /data/scratch"                      = "true"
  )

  for (old in names(replacements)) {
    lines <- gsub(old, replacements[[old]], lines, fixed = TRUE)
  }

  # Write patched script to a temp file
  tmp <- tempfile(fileext = ".sh")
  writeLines(lines, tmp)
  Sys.chmod(tmp, mode = "0755")
  tmp
}


.base_theme <- function() {
  ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line     = ggplot2::element_line(linewidth = 0.4),
      axis.ticks    = ggplot2::element_line(linewidth = 0.4),
      plot.margin   = ggplot2::margin(8, 20, 8, 8),
      plot.subtitle = ggplot2::element_text(colour = "grey50", size = 9),
      plot.caption  = ggplot2::element_text(colour = "grey50", size = 9)
    )
}

# ── Birthday problem mathematics ──────────────────────────────────────────────

.expected_distinct <- function(n, K) {
  if (n <= 0 || K <= 0) return(0)
  n * (1 - (1 - 1 / n)^K)
}

.pct_dup <- function(n, K) {
  if (K <= 1) return(0)
  max(0, (K - .expected_distinct(n, K)) / K * 100)
}

.p_exactly_j <- function(n, K, j) {
  if (j > min(n, K) || j < 0) return(0)
  s <- sum(sapply(0:j, function(i) (-1)^i * choose(j, i) * (j - i)^K))
  choose(n, j) * s / n^K
}

.p_at_least_m <- function(n, K, m) {
  if (m > n || K < m) return(0)
  p <- sum(sapply(m:min(n, K), function(j) .p_exactly_j(n, K, j)))
  max(0, min(1, p))
}

# ── shRNA aggregation (shared between QC functions) ───────────────────────────

.shrna_agg <- function(clones_above_dis) {
  clones_above_dis %>%
    dplyr::group_by(name, gene) %>%
    dplyr::summarise(
      n_UCI   = dplyr::n_distinct(UCI),
      n_reads = sum(Freq),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      pct_UCI   = n_UCI   / sum(n_UCI)   * 100,
      pct_reads = n_reads / sum(n_reads) * 100
    )
}

# ── Section 1: clone-level summary stats ─────────────────────────────────────

.clone_stats <- function(clones, DIs) {
  total         <- nrow(clones)
  matched       <- sum(clones$name != "UNMATCHED")
  below_DIs     <- sum(clones$Freq <  DIs)
  above_DIs     <- sum(clones$Freq >= DIs)
  above_matched <- sum(clones$Freq >= DIs & clones$name != "UNMATCHED")

  breaks_b <- c(1, 2, 11, 51, 101, 501, 1001, 5001, 10001, Inf)
  labels_b <- c("1", "2-10", "11-50", "51-100", "101-500",
                "501-1,000", "1,001-5,000", "5,001-10,000", ">10,000")
  clones$bucket <- cut(clones$Freq, breaks = breaks_b, labels = labels_b,
                       right = FALSE, include.lowest = TRUE)
  bc <- clones %>%
    dplyr::count(bucket, name = "n_clones") %>%
    dplyr::mutate(pct = round(n_clones / total * 100, 1))

  list(total = total, matched = matched, unmatched = total - matched,
       below_DIs = below_DIs, above_DIs = above_DIs,
       above_matched = above_matched, bc = bc, clones = clones)
}

# ── Birthday analysis: per-target results table ───────────────────────────────

.birthday_results <- function(shrna, N, DIs, eff,
                               nucleofect_cells, transfect_clones,
                               transfect_cells) {
  m_targets <- 1:10
  K_planned <- (nucleofect_cells * eff) / N

  results <- dplyr::bind_rows(lapply(m_targets, function(m) {
    K_m      <- m
    pd       <- .pct_dup(mean(shrna$n_UCI), K_m)
    p_ok     <- .p_at_least_m(round(mean(shrna$n_UCI)), K_m, m)
    n_lib_ok <- sum(shrna$n_UCI >= m)
    verdict  <- dplyr::case_when(
      pd > 20 & mean(shrna$n_UCI) < m ~ "LOW UCI + HIGH DUP",
      pd > 20                          ~ "HIGH DUP (>20%)",
      mean(shrna$n_UCI) < m            ~ "LOW UCI in library",
      TRUE                             ~ "OK"
    )
    data.frame(
      m_target      = m,
      K             = K_m,
      cells_M       = round(K_m * N / eff / 1e6, 1),
      pct_dup_mean  = round(pd, 1),
      p_recover_mean = round(p_ok * 100, 1),
      n_shrnas_ok   = n_lib_ok,
      pct_shrnas_ok = round(n_lib_ok / N * 100, 1),
      verdict       = verdict
    )
  }))

  shrna <- shrna %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pct_dup_planned   = round(.pct_dup(n_UCI, K_planned), 1),
      exp_distinct_plan = round(.expected_distinct(n_UCI, K_planned), 2),
      dup_status = dplyr::case_when(
        pct_dup_planned <= 20 ~ "GOOD  (<=20% dup)",
        pct_dup_planned <= 40 ~ "FAIR  (<=40% dup)",
        TRUE                  ~ "POOR  (>40% dup)"
      )
    ) %>%
    dplyr::ungroup()

  list(results = results, shrna = shrna, K_planned = K_planned,
       n_good = sum(shrna$pct_dup_planned <= 20),
       n_fair = sum(shrna$pct_dup_planned > 20 & shrna$pct_dup_planned <= 40),
       n_poor = sum(shrna$pct_dup_planned > 40))
}
