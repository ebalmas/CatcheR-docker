#' Run combined two-library QC analysis
#'
#' Compares two independent preparations of the same shRNA plasmid library
#' (technical replicates) by merging their `distribution_all_clones.csv`
#' files, computing reproducibility metrics, and producing a go/no-go
#' recommendation for proceeding with the screen.
#'
#' @param lib1_dir         Path to library 1 `results/` folder.
#' @param lib2_dir         Path to library 2 `results/` folder.
#' @param out_dir          Output directory. Created if absent.
#'   Default `"combined_results/"`.
#' @param lib1_name        Display label for library 1. Default `"CATCHER1"`.
#' @param lib2_name        Display label for library 2. Default `"CATCHER2"`.
#' @param DIs              DIs threshold. Default `300`.
#' @param transfect_clones Clones from test nucleofection. Default `100`.
#' @param transfect_cells  Cells in test nucleofection. Default `2000000`.
#' @param nucleofect_cells Cells planned for actual experiment.
#'   Default `2000000`.
#'
#' @section Go/no-go criteria:
#' \describe{
#'   \item{Reproducibility}{Spearman rho of UCI counts >= 0.80.}
#'   \item{shRNA overlap}{>= 80% of shRNAs present in both libraries.}
#' }
#'
#' @return A named list (invisibly) containing:
#' \describe{
#'   \item{`comp`}{Per-shRNA comparison data frame (lib1 vs lib2).}
#'   \item{`merged`}{Merged library data frame.}
#'   \item{`birthday`}{Birthday analysis on merged library.}
#'   \item{`repro_pass`}{Logical — reproducibility criterion met.}
#'   \item{`balance_pass`}{Logical — overlap criterion met.}
#'   \item{`overall_go`}{Logical — both criteria met.}
#' }
#'
#' @section Output files written to `out_dir`:
#' \describe{
#'   \item{`combined_report.txt`}{Full go/no-go report.}
#'   \item{`combined_uci_correlation.pdf/.jpg`}{UCI count scatter plot.}
#'   \item{`combined_read_correlation.pdf/.jpg`}{Read % scatter plot.}
#'   \item{`combined_uci_distribution_overlay.pdf/.jpg`}{Distribution overlay.}
#'   \item{`combined_merged_3panel.pdf/.jpg`}{3-panel for merged library.}
#'   \item{`calc_combined_cells_and_dup.pdf/.jpg`}{Birthday: cells and dup %.}
#'   \item{`calc_combined_pct_dup_per_shrna.pdf/.jpg`}{Birthday: per shRNA.}
#'   \item{`calc_combined_recovery_prob.pdf/.jpg`}{Birthday: recovery prob.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- catcheR_step2QC_combinedQC(
#'   lib1_dir         = "../CATCHER1_new/results/",
#'   lib2_dir         = "../CATCHER2_new/results/",
#'   out_dir          = "combined_results/",
#'   DIs              = 300,
#'   transfect_clones = 100,
#'   transfect_cells  = 2000000,
#'   nucleofect_cells = 2000000
#' )
#' result$overall_go   # TRUE = proceed, FALSE = no-go
#' }
#'
#' @seealso [catcheR_step2QC_extraction()], [catcheR_step2QC_plasmidQC()], [catcheR_step2QC_plots()]
#' @export
catcheR_step2QC_combinedQC <- function(lib1_dir         = "CATCHER1/results/",
                                        lib2_dir         = "CATCHER2/results/",
                                        output_dir       = "Output/",
                                        sample_name      = NULL,
                                        lib1_name        = "CATCHER1",
                                        lib2_name        = "CATCHER2",
                                        DIs              = 300,
                                        transfect_clones = 100,
                                        transfect_cells  = 2000000,
                                        nucleofect_cells = 2000000) {

  eff  <- transfect_clones / transfect_cells
  dirs <- .setup_output_dirs(output_dir, "step2QC_combinedQC", sample_name)

  message("=============================================================")
  message("  catcheR::catcheR_step2QC_combinedQC()")
  message("=============================================================")
  message(sprintf("  Library 1 : %s  (%s)", lib1_name, lib1_dir))
  message(sprintf("  Library 2 : %s  (%s)", lib2_name, lib2_dir))
  message(sprintf("  Output    : %s",        dirs$root))
  message(sprintf("  DIs=%d | eff=%d/%s | planned=%s cells",
                  DIs, transfect_clones,
                  format(transfect_cells,  big.mark = ","),
                  format(nucleofect_cells, big.mark = ",")))

  # ── Load ──────────────────────────────────────────────────────────────────
  raw1 <- .load_clones(lib1_dir); raw1$library <- lib1_name
  raw2 <- .load_clones(lib2_dir); raw2$library <- lib2_name
  above1 <- raw1 %>% dplyr::filter(Freq >= DIs, name != "UNMATCHED")
  above2 <- raw2 %>% dplyr::filter(Freq >= DIs, name != "UNMATCHED")
  s1 <- .shrna_agg(above1)
  s2 <- .shrna_agg(above2)

  # ── Overlap ───────────────────────────────────────────────────────────────
  message("\n-- Overlap --")
  shared_shrnas <- length(intersect(s1$name, s2$name))
  ucis1         <- paste(above1$barcode, above1$UCI, sep = "_")
  ucis2         <- paste(above2$barcode, above2$UCI, sep = "_")
  shared_ucis   <- length(intersect(ucis1, ucis2))
  union_ucis    <- length(union(ucis1, ucis2))
  balance_pass  <- shared_shrnas / max(nrow(s1), nrow(s2)) >= 0.80
  message(sprintf("  shRNAs shared: %d/%d (%.1f%%) [%s]",
                  shared_shrnas, max(nrow(s1), nrow(s2)),
                  shared_shrnas / max(nrow(s1), nrow(s2)) * 100,
                  if (balance_pass) "PASS" else "FAIL"))
  message(sprintf("  UCIs shared  : %d/%d (%.1f%% of union)",
                  shared_ucis, union_ucis, shared_ucis / union_ucis * 100))

  # ── Correlation ───────────────────────────────────────────────────────────
  message("\n-- Correlation --")
  comp <- dplyr::full_join(
    s1 %>% dplyr::select(name, gene,
                          n_UCI_1 = n_UCI, pct_UCI_1 = pct_UCI,
                          n_reads_1 = n_reads, pct_reads_1 = pct_reads),
    s2 %>% dplyr::select(name, gene,
                          n_UCI_2 = n_UCI, pct_UCI_2 = pct_UCI,
                          n_reads_2 = n_reads, pct_reads_2 = pct_reads),
    by = c("name","gene")) %>%
    dplyr::mutate(
      dplyr::across(where(is.numeric), ~ tidyr::replace_na(., 0)),
      in_both = n_UCI_1 > 0 & n_UCI_2 > 0)

  cor_uci_s   <- stats::cor(comp$n_UCI_1,     comp$n_UCI_2,     method = "spearman")
  cor_uci_p   <- stats::cor(comp$n_UCI_1,     comp$n_UCI_2,     method = "pearson")
  cor_reads_s <- stats::cor(comp$pct_reads_1, comp$pct_reads_2, method = "spearman")
  cor_reads_p <- stats::cor(comp$pct_reads_1, comp$pct_reads_2, method = "pearson")
  repro_pass  <- cor_uci_s >= 0.80
  message(sprintf("  UCI Spearman rho=%.3f  Pearson r=%.3f  [%s]",
                  cor_uci_s, cor_uci_p, if (repro_pass) "PASS" else "FAIL"))

  # Correlation plots
  p_uci_cor <- ggplot2::ggplot(comp,
      ggplot2::aes_string(x = "n_UCI_1", y = "n_UCI_2", colour = "in_both")) +
    ggplot2::geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey60") +
    ggplot2::geom_point(alpha=0.75,size=2.5) +
    ggplot2::annotate("text",x=-Inf,y=Inf,hjust=-0.1,vjust=1.4,size=3.5,
                      colour="grey30",
                      label=sprintf("Pearson r=%.3f\nSpearman rho=%.3f",
                                    cor_uci_p, cor_uci_s)) +
    ggplot2::scale_colour_manual(
      values=c("TRUE"="#378ADD","FALSE"="#D85A30"),
      labels=c("TRUE"="in both","FALSE"="one library only"), name=NULL) +
    ggplot2::labs(x=paste("UCIs/shRNA -",lib1_name),
                  y=paste("UCIs/shRNA -",lib2_name),
                  title="UCI count reproducibility") +
    .base_theme() + ggplot2::theme(legend.position="top")
  .save_plot(p_uci_cor, dirs, "combined_uci_correlation")

  p_reads_cor <- ggplot2::ggplot(comp,
      ggplot2::aes_string(x="pct_reads_1",y="pct_reads_2",colour="in_both")) +
    ggplot2::geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey60") +
    ggplot2::geom_point(alpha=0.75,size=2.5) +
    ggplot2::annotate("text",x=-Inf,y=Inf,hjust=-0.1,vjust=1.4,size=3.5,
                      colour="grey30",
                      label=sprintf("Pearson r=%.3f\nSpearman rho=%.3f",
                                    cor_reads_p, cor_reads_s)) +
    ggplot2::scale_colour_manual(
      values=c("TRUE"="#378ADD","FALSE"="#D85A30"),
      labels=c("TRUE"="in both","FALSE"="one library only"), name=NULL) +
    ggplot2::labs(x=paste("% reads -",lib1_name),
                  y=paste("% reads -",lib2_name),
                  title="Read representation reproducibility") +
    .base_theme() + ggplot2::theme(legend.position="top")
  .save_plot(p_reads_cor, dirs, "combined_read_correlation")

  # ── Merged library ────────────────────────────────────────────────────────
  message("\n-- Merged library --")
  merged <- dplyr::bind_rows(
    above1 %>% dplyr::mutate(uid = paste(barcode, UCI, sep="_")),
    above2 %>% dplyr::mutate(uid = paste(barcode, UCI, sep="_"))) %>%
    dplyr::group_by(name, gene) %>%
    dplyr::summarise(n_UCI         = dplyr::n_distinct(uid),
                     n_reads_total = sum(Freq), .groups="drop") %>%
    dplyr::mutate(pct_UCI   = n_UCI         / sum(n_UCI)         * 100,
                  pct_reads = n_reads_total / sum(n_reads_total) * 100)
  N_m <- nrow(merged)
  message(sprintf("  Merged: %d shRNAs | %d UCIs | mean=%.1f median=%.0f",
                  N_m, sum(merged$n_UCI),
                  mean(merged$n_UCI), stats::median(merged$n_UCI)))

  # Distribution overlay
  both_long <- dplyr::bind_rows(
    s1     %>% dplyr::mutate(library = lib1_name),
    s2     %>% dplyr::mutate(library = lib2_name),
    merged %>% dplyr::mutate(library = "MERGED (union)"))
  p_overlay <- ggplot2::ggplot(both_long,
      ggplot2::aes(x=pct_UCI, fill=library)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y=ggplot2::after_stat(count)/ggplot2::after_stat(ncount)*100),
      binwidth=0.5, position="identity", alpha=0.55, colour=NA, boundary=0) +
    ggplot2::geom_vline(
      data=both_long %>%
        dplyr::group_by(library) %>%
        dplyr::summarise(med=stats::median(pct_UCI),.groups="drop"),
      ggplot2::aes(xintercept=med,colour=library),
      linetype="dashed", linewidth=0.8) +
    ggplot2::scale_fill_manual(values=c("#378ADD","#BA7517","#1D9E75"),name=NULL)+
    ggplot2::scale_colour_manual(values=c("#378ADD","#BA7517","#1D9E75"),name=NULL)+
    ggplot2::labs(x="% of UCIs per shRNA", y="% of shRNAs (normalised)",
                  title="UCI distribution - individual libraries and merged") +
    ggplot2::theme_classic(base_size=13)+ggplot2::theme(legend.position="top")
  .save_plot(p_overlay, dirs, "combined_uci_distribution_overlay")

  # 3-panel for merged
  sw1 <- stats::shapiro.test(merged$pct_UCI)
  sw2 <- stats::shapiro.test(merged$pct_reads)
  sw3 <- stats::shapiro.test(merged$n_UCI)
  wx1 <- stats::wilcox.test(merged$pct_UCI,   mu=100/N_m)
  wx2 <- stats::wilcox.test(merged$pct_reads, mu=100/N_m)
  fc  <- "#1D9E75"; bt <- .base_theme()

  mp <- function(vals, tn, xl, tt, np, mp_p, bw=0.5) {
    df <- data.frame(x=vals)
    ggplot2::ggplot(df,ggplot2::aes(x=x)) +
      ggplot2::geom_histogram(
        ggplot2::aes(y=ggplot2::after_stat(count)/tn*100),
        binwidth=bw,fill=fc,colour=NA,boundary=0) +
      ggplot2::geom_vline(xintercept=stats::median(vals),
                          linetype="dashed",linewidth=0.6) +
      ggplot2::annotate("text",x=max(vals)*0.6,y=Inf,vjust=1.5,hjust=0,
                        size=3,colour="grey30",lineheight=1.3,
                        label=paste0("Normality\n",np,"\n\nMedian\n",mp_p)) +
      ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(0,0.02))) +
      ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,0.06))) +
      ggplot2::labs(x=xl,y="% of shRNAs",title=paste(tt,"- merged")) + bt
  }
  pA <- mp(merged$pct_UCI,  N_m,"% of UCIs per shRNA",  "A",
           .fmt_p(sw1$p.value),.fmt_p(wx1$p.value))
  pB <- ggplot2::ggplot(merged,ggplot2::aes(x=n_UCI))+
    ggplot2::geom_histogram(
      ggplot2::aes(y=ggplot2::after_stat(count)/N_m*100),
      binwidth=2,fill=fc,colour=NA,boundary=0.5)+
    ggplot2::geom_vline(xintercept=stats::median(merged$n_UCI),
                        linetype="dashed",linewidth=0.6)+
    ggplot2::scale_x_continuous(
      breaks=seq(0,max(merged$n_UCI)+2,by=4),
      expand=ggplot2::expansion(mult=c(0,0.02)))+
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,0.06)))+
    ggplot2::labs(x="Number of UCIs per shRNA",y="% of shRNAs",
                  title="B - merged")+bt
  pC <- mp(merged$pct_reads,N_m,"% of reads per shRNA", "C",
           .fmt_p(sw2$p.value),.fmt_p(wx2$p.value))

  p_3panel <- NULL
  if (requireNamespace("patchwork",quietly=TRUE)) {
    p_3panel <- (pA | pB | pC) +
      patchwork::plot_annotation(
        title=sprintf("Merged library representation  (DIs=%d)",DIs),
        subtitle=sprintf("%d shRNAs | %d UCIs | mean=%.1f | median=%.0f",
                         N_m, sum(merged$n_UCI),
                         mean(merged$n_UCI), stats::median(merged$n_UCI)),
        theme=ggplot2::theme(
          plot.title=ggplot2::element_text(size=13,face="plain"),
          plot.subtitle=ggplot2::element_text(size=10,colour="grey50")))
    .save_plot(p_3panel, dirs, "combined_merged_3panel", w=18, h=5)
  }

  # ── Birthday on merged ────────────────────────────────────────────────────
  message("\n-- Birthday analysis (merged) --")
  b_m <- .birthday_results(merged, N_m, DIs, eff,
                             nucleofect_cells, transfect_clones, transfect_cells)

  # Birthday plots
  rp_m <- b_m$results %>%
    dplyr::mutate(
      m_label  = paste0("m=",m_target),
      bar_fill = dplyr::case_when(
        pct_dup_mean<=20 & pct_shrnas_ok>=80 ~ "OK (<20% dup)",
        pct_dup_mean<=20                      ~ "Library insufficient",
        TRUE                                   ~ ">20% duplicates"))

  p_bc <- ggplot2::ggplot(rp_m,
      ggplot2::aes(x=factor(m_label,levels=paste0("m=",1:10)),
                   y=cells_M,fill=bar_fill)) +
    ggplot2::geom_col(width=0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label=sprintf("%.1f%%\ndup\n%d%%\nrecov.",
                                  pct_dup_mean,round(p_recover_mean))),
      vjust=-0.15,size=2.7,lineheight=1.1,colour="grey20") +
    ggplot2::geom_hline(yintercept=nucleofect_cells/1e6,
                        linetype="dotted",colour="black",linewidth=0.9) +
    ggplot2::scale_fill_manual(
      values=c("OK (<20% dup)"="#1D9E75","Library insufficient"="#BA7517",
               ">20% duplicates"="#D85A30"),name=NULL) +
    ggplot2::scale_y_continuous(labels=scales::label_number(suffix="M"),
                                expand=ggplot2::expansion(mult=c(0,0.20))) +
    ggplot2::labs(x="Target distinct UCIs per shRNA (m)",
                  y="Cells to nucleofect (millions)",
                  title="Merged library - cells needed and duplication rate") +
    .base_theme()+ggplot2::theme(legend.position="top")
  .save_plot(p_bc, dirs, "calc_combined_cells_and_dup")

  so_m <- b_m$shrna %>% dplyr::arrange(n_UCI) %>%
    dplyr::mutate(name=factor(name,levels=name))
  p_dm <- ggplot2::ggplot(so_m,
      ggplot2::aes(x=name,y=pct_dup_planned,fill=dup_status)) +
    ggplot2::geom_col(width=0.85) +
    ggplot2::geom_hline(yintercept=20,linetype="dashed",
                        colour="black",linewidth=0.7) +
    ggplot2::scale_fill_manual(
      values=c("GOOD  (<=20% dup)"="#1D9E75","FAIR  (<=40% dup)"="#BA7517",
               "POOR  (>40% dup)" ="#D85A30"),name=NULL) +
    ggplot2::scale_y_continuous(
      labels=scales::label_number(suffix="%"),
      limits=c(0,max(b_m$shrna$pct_dup_planned)*1.15),
      expand=ggplot2::expansion(mult=c(0,0))) +
    ggplot2::labs(x="shRNA (ordered by merged UCI count)",
                  y="% duplicate clones",
                  title="% dup per shRNA at planned nucleofection (merged)") +
    .base_theme() +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,
                                                      vjust=0.5,size=6),
                   legend.position="top")
  .save_plot(p_dm, dirs, "calc_combined_pct_dup_per_shrna", w=16, h=6)

  rec_m <- dplyr::bind_rows(lapply(2:5, function(m) {
    b_m$shrna %>% dplyr::mutate(
      m_target=m, m_label=paste0("m=",m),
      p_rec=sapply(n_UCI,function(n).p_at_least_m(n,m,m)*100),
      lib_ok=n_UCI>=m)
  }))
  p_rm <- ggplot2::ggplot(
    rec_m %>% dplyr::mutate(
      name=factor(name,
                  levels=b_m$shrna%>%dplyr::arrange(n_UCI)%>%dplyr::pull(name))),
    ggplot2::aes(x=name,y=p_rec,fill=lib_ok)) +
    ggplot2::geom_col(width=0.85) +
    ggplot2::geom_hline(yintercept=80,linetype="dashed",
                        colour="black",linewidth=0.6) +
    ggplot2::facet_wrap(~m_label,ncol=2) +
    ggplot2::scale_fill_manual(
      values=c("TRUE"="#1D9E75","FALSE"="#D85A30"),
      labels=c("TRUE"="n>=m","FALSE"="n<m (impossible)"),name=NULL) +
    ggplot2::scale_y_continuous(labels=scales::label_number(suffix="%"),
                                limits=c(0,105),
                                expand=ggplot2::expansion(mult=c(0,0))) +
    ggplot2::labs(x="shRNA",y="P(recover >= m distinct UCIs) %",
                  title="Recovery probability per shRNA (merged library)") +
    ggplot2::theme_classic(base_size=11) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,
                                                      vjust=0.5,size=4),
                   legend.position="top",
                   strip.background=ggplot2::element_rect(fill="grey95",colour=NA))
  .save_plot(p_rm, dirs, "calc_combined_recovery_prob", w=14, h=10)

  # ── Save CSVs ──────────────────────────────────────────────────────────────
  .save_csv(comp,       dirs, "comparison_lib1_vs_lib2")
  .save_csv(merged,     dirs, "merged_library_summary")
  .save_csv(b_m$results,dirs, paste0("birthday_results_DIs",DIs))
  .save_csv(b_m$shrna,  dirs, paste0("birthday_per_shrna_DIs",DIs))

  # ── Go/no-go report ────────────────────────────────────────────────────────
  overall_go <- repro_pass & balance_pass
  r <- b_m$results
  report_lines <- c(
    "=============================================================",
    "  COMBINED LIBRARY QC REPORT",
    "=============================================================",
    sprintf("  Library 1  : %s  (%s)", lib1_name, lib1_dir),
    sprintf("  Library 2  : %s  (%s)", lib2_name, lib2_dir),
    sprintf("  DIs=%d | eff=%d/%s | planned=%s cells",
            DIs, transfect_clones,
            format(transfect_cells,big.mark=","),
            format(nucleofect_cells,big.mark=",")),
    "-------------------------------------------------------------",
    "",
    "REPRODUCIBILITY",
    sprintf("  UCI Spearman rho : %.3f  [%s >= 0.80]",
            cor_uci_s, if(repro_pass)"PASS" else "FAIL"),
    sprintf("  UCI Pearson  r   : %.3f", cor_uci_p),
    "",
    "OVERLAP",
    sprintf("  shRNAs shared    : %d/%d (%.1f%%)  [%s >= 80%%]",
            shared_shrnas, max(nrow(s1),nrow(s2)),
            shared_shrnas/max(nrow(s1),nrow(s2))*100,
            if(balance_pass)"PASS" else "FAIL"),
    sprintf("  UCIs shared      : %d/%d (%.1f%% of union)",
            shared_ucis, union_ucis, shared_ucis/union_ucis*100),
    "",
    "MERGED LIBRARY",
    sprintf("  %d shRNAs | %d UCIs | mean=%.1f median=%.0f range=%d-%d",
            N_m, sum(merged$n_UCI), mean(merged$n_UCI),
            stats::median(merged$n_UCI),
            min(merged$n_UCI), max(merged$n_UCI)),
    "",
    "BIRTHDAY ANALYSIS (merged)",
    sprintf("  %-4s %-8s %-8s %-12s %s","m","cells(M)","% dup","P(>=m UCI)","verdict")
  )
  for (i in seq_len(nrow(r)))
    report_lines <- c(report_lines,
      sprintf("  %-4d %-8.1f %-7.1f%% %-11.1f%%  %s",
              r$m_target[i],r$cells_M[i],r$pct_dup_mean[i],
              r$p_recover_mean[i],r$verdict[i]))
  report_lines <- c(report_lines,
    "",
    "-------------------------------------------------------------",
    "DECISION",
    "",
    sprintf("  [%s] Reproducibility (Spearman rho >= 0.80): %.3f",
            if(repro_pass)"PASS" else "FAIL", cor_uci_s),
    sprintf("  [%s] shRNA overlap  (>= 80%% in both)       : %.1f%%",
            if(balance_pass)"PASS" else "FAIL",
            shared_shrnas/max(nrow(s1),nrow(s2))*100),
    "",
    sprintf("  OVERALL: %s",
            if(overall_go)">>> GO - proceed with the screen <<<" else
              ">>> NO-GO - address issues above <<<"),
    "=============================================================")
  .save_stats(report_lines, dirs, "combined_report")

  message("\n  RECOMMENDATION: ", if(overall_go)"GO" else "NO-GO")

  # ── Assemble result ────────────────────────────────────────────────────────
  result <- list(
    # Data frames
    comp          = comp,
    merged        = merged,
    shrna_lib1    = s1,
    shrna_lib2    = s2,
    birthday      = b_m$results,
    shrna_birthday = b_m$shrna,
    # Plots
    plots = list(
      uci_correlation    = p_uci_cor,
      reads_correlation  = p_reads_cor,
      distribution_overlay = p_overlay,
      threepanel         = p_3panel,
      cells_dup          = p_bc,
      pct_dup            = p_dm,
      recovery           = p_rm
    ),
    # Decision
    repro_pass   = repro_pass,
    balance_pass = balance_pass,
    overall_go   = overall_go,
    # Summary
    summary = list(
      cor_uci_spearman  = cor_uci_s,
      cor_uci_pearson   = cor_uci_p,
      cor_reads_spearman = cor_reads_s,
      shared_shrnas     = shared_shrnas,
      shared_ucis       = shared_ucis,
      N_merged          = N_m,
      mean_uci_merged   = mean(merged$n_UCI),
      K_planned         = b_m$K_planned
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

  .save_r_objects(result, dirs, "step2QC_combinedQC")

  # to_scratch: merged summary needed if running further QC
  .to_scratch(
    c(file.path(dirs$csv, "merged_library_summary.csv"),
      file.path(dirs$csv, paste0("birthday_results_DIs",DIs,".csv")),
      file.path(dirs$stats, "combined_report.txt")),
    dirs
  )

  message(sprintf("\n  Output in: %s", dirs$root))
  message("=============================================================")

  invisible(result)
}
