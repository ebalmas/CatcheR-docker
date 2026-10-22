#' Generate visualisation plots from Python pipeline outputs
#'
#' Produces publication-ready QC and distribution plots by reading the
#' CSV/TSV files created by [catcheR_step2QC_extraction()].
#' CSV/TSV files created by [run_extraction()].
#'
#' @param results_dir Path to the `results/` folder. Default `"results/"`.
#' @param DIs         DIs threshold. Default `300`.
#' @param clones_path Optional path to a `clones_of_interest.csv` file.
#'   If `NULL` (default), the function looks for
#'   `results_dir/clones_of_interest.csv`; if absent, those plots are skipped.
#'
#' @return Invisibly returns `NULL`. All output is written to disk.
#'
#' @section Required packages:
#' Requires `data.table` and `forcats` in addition to the standard imports.
#' Install with `install.packages(c("data.table","forcats"))`.
#'
#' @section Output files written to `results_dir`:
#' \describe{
#'   \item{`clone_percentage_all.pdf/.jpg`}{Histogram: % reads per clone.}
#'   \item{`clone_distribution_above_DIs.pdf/.jpg`}{Bar: reads per clone.}
#'   \item{`density_above_DIs.pdf/.jpg`}{Density + histogram.}
#'   \item{`pie_chart_clone.pdf/.jpg`}{Pie by shRNA clone.}
#'   \item{`pie_chart_gene.pdf/.jpg`}{Pie aggregated by gene.}
#'   \item{`barcode_reads_barplot.pdf/.jpg`}{Reads per barcode.}
#'   \item{`gene_reads_barplot.pdf/.jpg`}{Reads per gene.}
#'   \item{`clones_of_interest_barplot.pdf/.jpg`}{If COI file present.}
#'   \item{`clones_of_interest_pie.pdf/.jpg`}{If COI file present.}
#' }
#'
#' @section Custom colours:
#' If `results_dir/colors.csv` exists (2-column CSV: `name,hex_color`),
#' its colour mapping is applied to all shRNA-coloured plots.
#'
#' @examples
#' \dontrun{
#' catcheR_step2QC_plots("results/", DIs = 300)
#' }
#'
#' @seealso [catcheR_step2QC_extraction()], [catcheR_step2QC_plasmidQC()]
#' @export
catcheR_step2QC_plots <- function(results_dir = "results/",
                               DIs         = 300,
                               clones_path = NULL) {

  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' required: install.packages('data.table')",
         call. = FALSE)
  if (!requireNamespace("forcats", quietly = TRUE))
    stop("Package 'forcats' required: install.packages('forcats')",
         call. = FALSE)

  message("=============================================================")
  message("  catcheR::catcheR_step2QC_plots()")
  message("=============================================================")
  message(sprintf("  Results dir  : %s", results_dir))
  message(sprintf("  DIs threshold: %d", DIs))

  # Internal helpers scoped to this function
  sp <- function(p, name, w = 20, h = 15) {
    base <- file.path(results_dir, name)
    ggplot2::ggsave(paste0(base, ".pdf"), plot=p, width=w, height=h, limitsize=FALSE)
    ggplot2::ggsave(paste0(base, ".jpg"), plot=p, width=w, height=h, limitsize=FALSE)
    message("  Saved: ", name, ".pdf/.jpg")
  }
  ac <- function(p) .apply_colours(p, results_dir)

  # Load data
  clone_t <- data.table::fread(
    file.path(results_dir, "distribution_all_clones.csv"), data.table = FALSE)
  stopifnot(all(c("clone","barcode","UCI","name","gene","Freq") %in% names(clone_t)))

  pct_file <- file.path(results_dir, "percentages.csv")
  clone_t_per <- if (file.exists(pct_file)) {
    data.table::fread(pct_file, data.table = FALSE)
  } else {
    clone_t %>%
      dplyr::filter(Freq > DIs, name != "UNMATCHED") %>%
      dplyr::group_by(barcode, name) %>%
      dplyr::summarise(Freq = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(percentage = Freq / sum(Freq) * 100)
  }
  above_DIs <- dplyr::filter(clone_t, Freq > DIs, name != "UNMATCHED")

  # 1. Percentage histogram
  message("  Plotting: clone_percentage_all")
  sp(ggplot2::ggplot(clone_t_per, ggplot2::aes(x = percentage)) +
    ggplot2::geom_histogram(bins=30, fill="#4682B4", colour="white") +
    ggplot2::labs(x="Percentage of assigned reads (%)", y="Number of clones",
                  title=paste0("Clone representation - DIs=", DIs)) +
    ggplot2::theme_minimal(base_size=14),
    "clone_percentage_all")

  # 2. Clone distribution bar
  message("  Plotting: clone_distribution_above_DIs")
  sp(ac(ggplot2::ggplot(above_DIs) +
    ggplot2::geom_col(
      ggplot2::aes(x=forcats::fct_reorder(clone, Freq, .desc=TRUE),
                   y=Freq, fill=name), width=0.85) +
    ggplot2::labs(x="Clone (barcode_UCI)", y="Read count", fill="shRNA",
                  title=paste0("Clone distribution - Freq > ", DIs)) +
    ggplot2::theme_minimal(base_size=12) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5,size=6))),
  "clone_distribution_above_DIs")

  # 3. Density
  message("  Plotting: density_above_DIs")
  sp(ac(ggplot2::ggplot(above_DIs, ggplot2::aes(x=Freq)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y=ggplot2::after_stat(density), fill=name),
      binwidth=max(1, diff(range(above_DIs$Freq))/40),
      colour="white", alpha=0.8) +
    ggplot2::geom_density(colour="black", linewidth=0.7) +
    ggplot2::labs(x="Reads per clone", y="Density", fill="shRNA",
                  title=paste0("Reads per clone - Freq > ", DIs)) +
    ggplot2::theme_minimal(base_size=14)),
  "density_above_DIs")

  # 4. Pie — clone
  message("  Plotting: pie_chart_clone")
  sp(ac(ggplot2::ggplot(above_DIs, ggplot2::aes(x="", y=Freq, fill=name)) +
    ggplot2::geom_bar(width=1, stat="identity") +
    ggplot2::coord_polar("y", start=0) +
    ggplot2::labs(fill="shRNA",
                  title=paste0("Read distribution by clone - DIs=", DIs)) +
    ggplot2::theme_minimal(base_size=14) +
    ggplot2::theme(axis.text=ggplot2::element_blank(),
                   axis.title=ggplot2::element_blank(),
                   panel.grid=ggplot2::element_blank())),
  "pie_chart_clone", w=15, h=15)

  # 5. Pie — gene
  message("  Plotting: pie_chart_gene")
  gene_agg <- above_DIs %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(Freq = sum(Freq), .groups = "drop")
  sp(ggplot2::ggplot(gene_agg, ggplot2::aes(x="", y=Freq, fill=gene)) +
    ggplot2::geom_bar(width=1, stat="identity") +
    ggplot2::coord_polar("y", start=0) +
    ggplot2::labs(fill="Gene",
                  title=paste0("Read distribution by gene - DIs=", DIs)) +
    ggplot2::theme_minimal(base_size=14) +
    ggplot2::theme(axis.text=ggplot2::element_blank(),
                   axis.title=ggplot2::element_blank(),
                   panel.grid=ggplot2::element_blank()),
  "pie_chart_gene", w=15, h=15)

  # 6. Barcode bar
  bc_f <- file.path(results_dir, "counts_per_barcode.tsv")
  if (file.exists(bc_f)) {
    message("  Plotting: barcode_reads_barplot")
    cbc <- utils::read.delim(bc_f, stringsAsFactors=FALSE) %>%
      dplyr::filter(reads > 0) %>% dplyr::arrange(dplyr::desc(reads))
    sp(ggplot2::ggplot(cbc,
        ggplot2::aes(x=forcats::fct_reorder(barcode,reads,.desc=TRUE),
                     y=reads, fill=gene)) +
      ggplot2::geom_col(width=0.8) +
      ggplot2::geom_text(ggplot2::aes(label=unique_UCIs), vjust=-0.3, size=2.5) +
      ggplot2::labs(x="Barcode (shRNA)", y="Total reads", fill="Gene",
                    title="Reads per barcode (label = unique UCIs)") +
      ggplot2::theme_minimal(base_size=12) +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,
                                                        vjust=0.5,size=7)),
    "barcode_reads_barplot")
  }

  # 7. Gene bar
  gf <- file.path(results_dir, "counts_per_gene.tsv")
  if (file.exists(gf)) {
    message("  Plotting: gene_reads_barplot")
    cg <- utils::read.delim(gf, stringsAsFactors=FALSE)
    sp(ggplot2::ggplot(cg,
        ggplot2::aes(x=forcats::fct_reorder(gene,reads,.desc=TRUE),
                     y=reads, fill=gene)) +
      ggplot2::geom_col(width=0.8, show.legend=FALSE) +
      ggplot2::geom_text(ggplot2::aes(label=unique_UCIs), vjust=-0.3, size=3) +
      ggplot2::labs(x="Gene", y="Total reads",
                    title="Reads per gene (label = unique UCIs)") +
      ggplot2::theme_minimal(base_size=14) +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1)),
    "gene_reads_barplot", w=14, h=10)
  }

  # 8. Clones of interest (optional)
  coi_file <- file.path(results_dir, "clones_of_interest.csv")
  if (!is.null(clones_path) && file.exists(clones_path)) coi_file <- clones_path
  if (file.exists(coi_file)) {
    message("  Plotting: clones_of_interest")
    coi <- data.table::fread(coi_file, data.table=FALSE)
    sp(ac(ggplot2::ggplot(coi) +
      ggplot2::geom_col(
        ggplot2::aes(x=forcats::fct_reorder(clone,Freq,.desc=TRUE),
                     y=Freq, fill=name), width=0.85) +
      ggplot2::labs(x="Clone", y="Read count", fill="shRNA",
                    title="Clones of interest - all reads") +
      ggplot2::theme_minimal(base_size=12) +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,
                                                        vjust=0.5,size=7))),
    "clones_of_interest_barplot")
    sp(ac(ggplot2::ggplot(coi, ggplot2::aes(x="",y=Freq,fill=name)) +
      ggplot2::geom_bar(width=1, stat="identity") +
      ggplot2::coord_polar("y", start=0) +
      ggplot2::labs(fill="shRNA", title="Clones of interest - read distribution") +
      ggplot2::theme_minimal(base_size=14) +
      ggplot2::theme(axis.text=ggplot2::element_blank(),
                     axis.title=ggplot2::element_blank(),
                     panel.grid=ggplot2::element_blank())),
    "clones_of_interest_pie", w=15, h=15)
  } else {
    message("  No clones-of-interest file found -- skipping those plots.")
  }

  message("\n=== Done. All plots written to: ", results_dir, " ===")
  invisible(NULL)
}
