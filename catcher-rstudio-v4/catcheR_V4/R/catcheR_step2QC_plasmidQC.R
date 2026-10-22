#' Run single-library plasmid QC analysis
#'
#' Reads `distribution_all_clones.csv` from the results folder produced by
#' [catcheR_step2QC_extraction()] and runs five analysis sections:
#'
#' 1. **Clone-level statistics** — total clones, matched %, above DIs %.
#' 2. **TXT report** — detailed read-count breakdown saved to disk.
#' 3. **Clone frequency histogram** — log-scale plot distinguishing noise
#'    from real clones.
#' 4. **shRNA representation 3-panel** — % UCIs / raw UCI count /
#'    % reads per shRNA.
#' 5. **Birthday-problem analysis** — for each target number of distinct
#'    UCIs per shRNA (m = 1..10): duplication rate, cells to nucleofect,
#'    and probability of recovery.
#'
#' @param results_dir      Path to the `results/` folder from [catcheR_step2QC_extraction()].
#'   Must contain `distribution_all_clones.csv`. Default `"results/"`.
#' @param DIs              DIs threshold — same value used in [catcheR_step2QC_extraction()].
#'   Clones with fewer reads are treated as artefacts. Default `300`.
#' @param transfect_clones Number of clones recovered in your **test**
#'   nucleofection experiment. Used to estimate efficiency. Default `100`.
#' @param transfect_cells  Number of cells nucleofected in that test.
#'   Default `2000000`.
#' @param nucleofect_cells Number of cells you **plan** to nucleofect in the
#'   actual experiment. Determines expected clones per shRNA (K).
#'   Default `2000000`.
#'
#' @details
#' **K is never a user input.** It is derived as:
#' \deqn{K = \frac{\text{nucleofect\_cells} \times \text{efficiency}}{N_{\text{shRNAs}}}}
#' where efficiency = `transfect_clones / transfect_cells`.
#'
#' **Birthday problem model.** Each recovered clone independently received
#' one plasmid at random from the pool. If a shRNA has *n* distinct UCIs in
#' the library and *K* clones are recovered for it, some will share a UCI
#' (duplicates). The expected number of distinct UCIs is:
#' \deqn{E[\text{distinct}] = n \left(1 - \left(1 - \frac{1}{n}\right)^K\right)}
#'
#' @return A named list (returned invisibly) containing:
#' \describe{
#'   \item{`stats`}{Clone-level summary statistics (list).}
#'   \item{`shrna`}{Per-shRNA data frame with UCI counts and birthday metrics.}
#'   \item{`birthday`}{Birthday analysis results including `results` table
#'     and per-shRNA `shrna` data frame.}
#' }
#'
#' @section Output files written to `results_dir`:
#' \describe{
#'   \item{`DIs_threshold_stats_DIs<N>.txt`}{Clone count breakdown.}
#'   \item{`clone_frequency_dist_DIs<N>.pdf/.jpg`}{Log-scale histogram.}
#'   \item{`shrna_representation_3panel_DIs<N>.pdf/.jpg`}{3-panel figure.}
#'   \item{`UCI_birthday_analysis_DIs<N>.txt`}{Birthday analysis report.}
#'   \item{`calc_cells_and_dup_DIs<N>.pdf/.jpg`}{Cells and dup rate per m.}
#'   \item{`calc_pct_dup_per_shrna_DIs<N>.pdf/.jpg`}{% dup per shRNA.}
#'   \item{`calc_recovery_prob_DIs<N>.pdf/.jpg`}{Recovery probability.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- catcheR_step2QC_plasmidQC(
#'   results_dir      = "CATCHER1_new/results/",
#'   DIs              = 300,
#'   transfect_clones = 100,
#'   transfect_cells  = 2000000,
#'   nucleofect_cells = 2000000
#' )
#'
#' # Inspect birthday results programmatically
#' result$birthday$results
#' result$shrna
#' }
#'
#' @seealso [catcheR_step2QC_extraction()], [catcheR_step2QC_combinedQC()], [catcheR_step2QC_plots()]
#' @export
catcheR_step2QC_plasmidQC <- function(results_dir      = "results/",
                                       output_dir       = "Output/",
                                       sample_name      = NULL,
                                       DIs              = 300,
                                       transfect_clones = 100,
                                       transfect_cells  = 2000000,
                                       nucleofect_cells = 2000000) {

  eff <- transfect_clones / transfect_cells

  # ── Setup output directories ──────────────────────────────────────────────
  dirs <- .setup_output_dirs(output_dir, "step2QC_plasmidQC", sample_name)

  message("=============================================================")
  message("  catcheR::catcheR_step2QC_plasmidQC()")
  message("=============================================================")
  message(sprintf("  Input (results_dir): %s", results_dir))
  message(sprintf("  DIs                : %d", DIs))
  message(sprintf("  Efficiency         : %d per %s cells",
                  transfect_clones, format(transfect_cells, big.mark = ",")))
  message(sprintf("  Planned nucleofection: %s cells",
                  format(nucleofect_cells, big.mark = ",")))

  # ── Load clones ───────────────────────────────────────────────────────────
  clones <- .load_clones(results_dir)

  # ── Section 1: clone-level statistics ─────────────────────────────────────
  message("\n-- Section 1: clone-level statistics --")
  s1 <- .clone_stats(clones, DIs)
  message(sprintf("  Total: %s | Above DIs: %s (%.1f%%) | Matched: %s (%.1f%%)",
                  format(s1$total,     big.mark = ","),
                  format(s1$above_DIs, big.mark = ","),
                  s1$above_DIs / s1$total * 100,
                  format(s1$matched,   big.mark = ","),
                  s1$matched   / s1$total * 100))

  # ── Section 2: stats TXT ──────────────────────────────────────────────────
  message("\n-- Section 2: stats report --")
  notes_map <- c("1" = "almost certainly artefacts", "2-10" = "likely noise",
                 "11-50" = "", "51-100" = "", "101-500" = "", "501-1,000" = "")
  stats_lines <- c(
    "=============================================================",
    "  DIs THRESHOLD STATISTICS",
    "=============================================================",
    sprintf("  DIs threshold  : %s reads", format(DIs, big.mark = ",")),
    "-------------------------------------------------------------",
    "",
    "OVERVIEW",
    sprintf("  Total clones         : %s", format(s1$total,         big.mark=",")),
    sprintf("  Matched barcodes     : %s  (%.1f%%)",
            format(s1$matched,   big.mark=","), s1$matched   / s1$total * 100),
    sprintf("  Unmatched            : %s  (%.1f%%)",
            format(s1$unmatched, big.mark=","), s1$unmatched / s1$total * 100),
    sprintf("  Below DIs (< %s)    : %s  (%.1f%%)  -- filtered",
            format(DIs, big.mark=","), format(s1$below_DIs, big.mark=","),
            s1$below_DIs / s1$total * 100),
    sprintf("  Above DIs (>= %s)   : %s  (%.1f%%)  -- kept",
            format(DIs, big.mark=","), format(s1$above_DIs, big.mark=","),
            s1$above_DIs / s1$total * 100),
    sprintf("    of which matched   : %s  (%.1f%%)",
            format(s1$above_matched, big.mark=","),
            s1$above_matched / s1$above_DIs * 100),
    "",
    "-------------------------------------------------------------",
    "READ COUNT BREAKDOWN",
    "",
    sprintf("  %-20s  %10s  %7s  %s", "Range", "Clones", "%", "Notes")
  )
  for (i in seq_len(nrow(s1$bc))) {
    lbl  <- as.character(s1$bc$bucket[i])
    note <- ifelse(lbl %in% names(notes_map), notes_map[[lbl]], "")
    stats_lines <- c(stats_lines,
      sprintf("  %-20s  %10s  %6.1f%%  %s",
              lbl, format(s1$bc$n_clones[i], big.mark=","),
              s1$bc$pct[i], note))
  }
  stats_lines <- c(stats_lines, "=============================================================")
  .save_stats(stats_lines, dirs, paste0("DIs_threshold_stats_DIs", DIs))

  # ── Section 3: histogram ──────────────────────────────────────────────────
  message("\n-- Section 3: clone frequency histogram --")
  cl3 <- s1$clones %>%
    dplyr::mutate(match_status = dplyr::if_else(name=="UNMATCHED","Unmatched","Matched"))
  p_hist <- ggplot2::ggplot(cl3, ggplot2::aes(x=Freq, fill=match_status)) +
    ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(count)),
                            bins=60, colour=NA, alpha=0.9) +
    ggplot2::geom_vline(xintercept=DIs, colour="#D85A30",
                        linewidth=0.8, linetype="dashed") +
    ggplot2::annotate("text", x=DIs*1.15, y=Inf,
                      label=paste0("DIs=",format(DIs,big.mark=",")),
                      hjust=0, vjust=1.5, size=3.5, colour="#D85A30") +
    ggplot2::scale_x_log10(labels=scales::label_comma(), breaks=10^(0:6)) +
    ggplot2::scale_y_log10(labels=scales::label_comma()) +
    ggplot2::scale_fill_manual(
      values=c("Matched"="#378ADD","Unmatched"="#B4B2A9"), name=NULL) +
    ggplot2::labs(x="Reads per clone", y="Number of clones (log scale)",
                  title=paste0("Clone frequency distribution  (DIs=",DIs,")")) +
    ggplot2::theme_minimal(base_size=13) +
    ggplot2::theme(legend.position="top")
  .save_plot(p_hist, dirs, paste0("clone_frequency_dist_DIs",DIs))

  # ── Section 4: 3-panel ────────────────────────────────────────────────────
  message("\n-- Section 4: shRNA representation 3-panel --")
  above <- s1$clones %>% dplyr::filter(Freq >= DIs, name != "UNMATCHED")
  shrna <- .shrna_agg(above)
  N     <- nrow(shrna)

  sw_uci   <- stats::shapiro.test(shrna$pct_UCI)
  sw_reads <- stats::shapiro.test(shrna$pct_reads)
  sw_n     <- stats::shapiro.test(shrna$n_UCI)
  wx_uci   <- stats::wilcox.test(shrna$pct_UCI,   mu=100/N)
  wx_reads <- stats::wilcox.test(shrna$pct_reads, mu=100/N)
  fill_col <- "#BA7517"
  bt       <- .base_theme()

  make_hp <- function(vals, tn, xl, tt, np, mp, bw=0.5) {
    df <- data.frame(x=vals)
    ggplot2::ggplot(df, ggplot2::aes(x=x)) +
      ggplot2::geom_histogram(
        ggplot2::aes(y=ggplot2::after_stat(count)/tn*100),
        binwidth=bw, fill=fill_col, colour=NA, boundary=0) +
      ggplot2::geom_vline(xintercept=stats::median(vals),
                          linetype="dashed", linewidth=0.6) +
      ggplot2::annotate("text", x=max(vals)*0.6, y=Inf,
                        vjust=1.5, hjust=0, size=3, colour="grey30",
                        lineheight=1.3,
                        label=paste0("Normality\n",np,"\n\nMedian\n",mp)) +
      ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(0,0.02))) +
      ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,0.06))) +
      ggplot2::labs(x=xl, y="% of shRNAs", title=tt) + bt
  }
  p_A <- make_hp(shrna$pct_UCI,   N, "% of UCIs per shRNA",
                  "A", .fmt_p(sw_uci$p.value),   .fmt_p(wx_uci$p.value))
  p_B <- ggplot2::ggplot(shrna, ggplot2::aes(x=n_UCI)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y=ggplot2::after_stat(count)/N*100),
      binwidth=2, fill=fill_col, colour=NA, boundary=0.5) +
    ggplot2::geom_vline(xintercept=stats::median(shrna$n_UCI),
                        linetype="dashed", linewidth=0.6) +
    ggplot2::scale_x_continuous(
      breaks=seq(0,max(shrna$n_UCI)+2,by=4),
      expand=ggplot2::expansion(mult=c(0,0.02))) +
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,0.06))) +
    ggplot2::labs(x="Number of UCIs per shRNA",y="% of shRNAs",title="B") + bt
  p_C <- make_hp(shrna$pct_reads, N, "% of reads per shRNA",
                  "C", .fmt_p(sw_reads$p.value), .fmt_p(wx_reads$p.value))

  p_3panel <- NULL
  if (requireNamespace("patchwork", quietly=TRUE)) {
    p_3panel <- (p_A | p_B | p_C) +
      patchwork::plot_annotation(
        title    = paste0("shRNA library representation  (DIs=",DIs,")"),
        subtitle = sprintf("%d shRNAs | %d UCIs | mean=%.1f | median=%.0f",
                           N, sum(shrna$n_UCI), mean(shrna$n_UCI),
                           stats::median(shrna$n_UCI)),
        theme = ggplot2::theme(
          plot.title    = ggplot2::element_text(size=13, face="plain"),
          plot.subtitle = ggplot2::element_text(size=10, colour="grey50")))
    .save_plot(p_3panel, dirs,
               paste0("shrna_representation_3panel_DIs",DIs), w=18, h=5)
  }

  # Save per-shRNA table
  .save_csv(shrna, dirs, paste0("shrna_summary_DIs",DIs))

  # ── Section 5: birthday analysis ─────────────────────────────────────────
  message("\n-- Section 5: birthday analysis --")
  b <- .birthday_results(shrna, N, DIs, eff,
                          nucleofect_cells, transfect_clones, transfect_cells)

  # Birthday TXT
  b_lines <- c(
    "=============================================================",
    "  UCI UNIQUENESS ANALYSIS - BIRTHDAY PROBLEM",
    "=============================================================",
    sprintf("  Efficiency        : %d per %s cells",
            transfect_clones, format(transfect_cells, big.mark=",")),
    sprintf("  Planned nucleofection : %s cells",
            format(nucleofect_cells, big.mark=",")),
    sprintf("  Expected K/shRNA  : %.2f clones", b$K_planned),
    sprintf("  Library           : %d shRNAs, mean=%.1f UCIs", N, mean(shrna$n_UCI)),
    "",
    "RESULTS PER TARGET m",
    "",
    sprintf("  %-4s  %-6s  %-10s  %-8s  %-12s  %-13s  %s",
            "m","K","cells(M)","% dup","P(>=m UCI)","% shRNAs OK","verdict"),
    strrep("-",72)
  )
  for (i in seq_len(nrow(b$results))) {
    r <- b$results[i,]
    b_lines <- c(b_lines,
      sprintf("  %-4d  %-6d  %-10.1f  %-7.1f%%  %-11.1f%%  %-12.1f%%  %s",
              r$m_target, r$K, r$cells_M, r$pct_dup_mean,
              r$p_recover_mean, r$pct_shrnas_ok, r$verdict))
  }
  .save_stats(b_lines, dirs, paste0("UCI_birthday_analysis_DIs",DIs))

  # Birthday plots
  rp <- b$results %>%
    dplyr::mutate(
      m_label  = paste0("m=", m_target),
      bar_fill = dplyr::case_when(
        pct_dup_mean <= 20 & pct_shrnas_ok >= 80 ~ "OK (<20% dup, lib OK)",
        pct_dup_mean <= 20                        ~ "Library insufficient",
        TRUE                                       ~ ">20% duplicates"))

  p_cells <- ggplot2::ggplot(
    rp, ggplot2::aes(x=factor(m_label,levels=paste0("m=",1:10)),
                     y=cells_M, fill=bar_fill)) +
    ggplot2::geom_col(width=0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label=sprintf("%.1f%%\ndup\n%d%%\nrecov.",
                                  pct_dup_mean, round(p_recover_mean))),
      vjust=-0.15, size=2.7, lineheight=1.1, colour="grey20") +
    ggplot2::geom_hline(yintercept=nucleofect_cells/1e6,
                        linetype="dotted", colour="black", linewidth=0.9) +
    ggplot2::scale_fill_manual(
      values=c("OK (<20% dup, lib OK)"="#1D9E75",
               "Library insufficient"="#BA7517",
               ">20% duplicates"="#D85A30"), name=NULL) +
    ggplot2::scale_y_continuous(
      labels=scales::label_number(suffix="M"),
      expand=ggplot2::expansion(mult=c(0,0.20))) +
    ggplot2::labs(x="Target distinct UCIs per shRNA (m)",
                  y="Cells to nucleofect (millions)",
                  title="Cells needed and duplication rate per target m") +
    .base_theme() + ggplot2::theme(legend.position="top")
  .save_plot(p_cells, dirs,
             paste0("calc_cells_and_dup_DIs",DIs), w=12, h=7)

  so <- b$shrna %>% dplyr::arrange(n_UCI) %>%
    dplyr::mutate(name=factor(name,levels=name))
  p_dup <- ggplot2::ggplot(
    so, ggplot2::aes(x=name, y=pct_dup_planned, fill=dup_status)) +
    ggplot2::geom_col(width=0.85) +
    ggplot2::geom_hline(yintercept=20, linetype="dashed",
                        colour="black", linewidth=0.7) +
    ggplot2::scale_fill_manual(
      values=c("GOOD  (<=20% dup)"="#1D9E75",
               "FAIR  (<=40% dup)"="#BA7517",
               "POOR  (>40% dup)" ="#D85A30"), name=NULL) +
    ggplot2::scale_y_continuous(
      labels=scales::label_number(suffix="%"),
      limits=c(0, max(b$shrna$pct_dup_planned)*1.15),
      expand=ggplot2::expansion(mult=c(0,0))) +
    ggplot2::labs(x="shRNA (ordered by n UCIs)", y="% duplicate clones",
                  title="% duplicate clones per shRNA at planned nucleofection") +
    .base_theme() +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5,size=6),
                   legend.position="top")
  .save_plot(p_dup, dirs,
             paste0("calc_pct_dup_per_shrna_DIs",DIs), w=16, h=6)

  rec_df <- dplyr::bind_rows(lapply(2:5, function(m) {
    b$shrna %>% dplyr::mutate(
      m_target=m, m_label=paste0("m=",m),
      p_rec=sapply(n_UCI, function(n) .p_at_least_m(n,m,m)*100),
      lib_ok=n_UCI>=m)
  }))
  p_rec <- ggplot2::ggplot(
    rec_df %>% dplyr::mutate(
      name=factor(name, levels=b$shrna%>%dplyr::arrange(n_UCI)%>%dplyr::pull(name))),
    ggplot2::aes(x=name, y=p_rec, fill=lib_ok)) +
    ggplot2::geom_col(width=0.85) +
    ggplot2::geom_hline(yintercept=80, linetype="dashed",
                        colour="black", linewidth=0.6) +
    ggplot2::facet_wrap(~m_label, ncol=2) +
    ggplot2::scale_fill_manual(
      values=c("TRUE"="#1D9E75","FALSE"="#D85A30"),
      labels=c("TRUE"="n >= m","FALSE"="n < m (impossible)"), name=NULL) +
    ggplot2::scale_y_continuous(
      labels=scales::label_number(suffix="%"),
      limits=c(0,105), expand=ggplot2::expansion(mult=c(0,0))) +
    ggplot2::labs(x="shRNA (ordered by UCI count)",
                  y="P(recover >= m distinct UCIs) %",
                  title="Probability of recovering >= m distinct UCIs per shRNA") +
    ggplot2::theme_classic(base_size=11) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5,size=4),
                   legend.position="top",
                   strip.background=ggplot2::element_rect(fill="grey95",colour=NA))
  .save_plot(p_rec, dirs,
             paste0("calc_recovery_prob_DIs",DIs), w=14, h=10)

  # Save birthday tables
  .save_csv(b$results, dirs, paste0("birthday_results_DIs",DIs))
  .save_csv(b$shrna,   dirs, paste0("birthday_per_shrna_DIs",DIs))

  # ── Assemble return object ────────────────────────────────────────────────
  result <- list(
    # Data frames
    clones       = s1$clones,
    clone_stats  = s1$bc,
    shrna        = shrna,
    birthday     = b$results,
    shrna_birthday = b$shrna,
    # Plots
    plots = list(
      histogram  = p_hist,
      panel_A    = p_A,
      panel_B    = p_B,
      panel_C    = p_C,
      threepanel = p_3panel,
      cells_dup  = p_cells,
      pct_dup    = p_dup,
      recovery   = p_rec
    ),
    # Summary values
    summary = list(
      N              = N,
      K_planned      = b$K_planned,
      n_good         = b$n_good,
      n_fair         = b$n_fair,
      n_poor         = b$n_poor,
      mean_uci       = mean(shrna$n_UCI),
      median_uci     = stats::median(shrna$n_UCI)
    ),
    # Paths
    paths = list(
      output     = dirs$root,
      csv        = dirs$csv,
      plots      = dirs$plots,
      stats      = dirs$stats,
      R_objects  = dirs$R_objects,
      to_scratch = dirs$to_scratch
    )
  )

  # ── Save R objects ────────────────────────────────────────────────────────
  .save_r_objects(result, dirs, "step2QC_plasmidQC")

  # ── to_scratch: copy files the next step needs ────────────────────────────
  .to_scratch(
    c(file.path(results_dir, "distribution_all_clones.csv"),
      file.path(dirs$csv, paste0("shrna_summary_DIs",DIs,".csv"))),
    dirs
  )

  message("\n=============================================================")
  message("  Done. Output in:")
  message(sprintf("    %s", dirs$root))
  message("=============================================================")

  invisible(result)
}
