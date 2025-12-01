#' Perform Infection-Level GSEA Analysis with Optional Timepoint Filtering
#'
#' This function performs differential expression and Hallmark GSEA analysis 
#' between infection levels within a Seurat object. It supports filtering by 
#' specific timepoints, ranking genes by \code{avg_log2FC} or adjusted p-values, 
#' and produces publication-ready NES barplots for each comparison.
#'
#' @param obj A \code{Seurat} object containing expression data and metadata.
#' @param group_col Character. Column name in \code{obj@meta.data} specifying 
#'   infection grouping (e.g. \code{"infection_level"}).
#' @param ref_level Character or numeric. Reference level to compare against.
#' @param min.pct Numeric. Minimum fraction of cells expressing a gene in 
#'   either group (passed to \code{Seurat::FindMarkers}).
#' @param logfc.threshold Numeric. Minimum log2 fold-change threshold 
#'   for differential testing (default = 0).
#' @param ranks_from Character. Method to rank genes for GSEA; one of 
#'   \code{"avg_log2FC"} or \code{"pval_adj_signed"} (default = \code{"avg_log2FC"}).
#' @param top_up Integer. Number of top upregulated pathways to display in each plot.
#' @param top_down Integer. Number of top downregulated pathways to display.
#' @param col_up Character. Fill color for upregulated pathways (default = \code{"#D74C4C"}).
#' @param col_down Character. Fill color for downregulated pathways (default = \code{"#3B6FB3"}).
#' @param drop_prefix Logical. Whether to remove the "HALLMARK_" prefix from pathway names.
#' @param save_plot Logical. Whether to save GSEA plots as PDF files.
#' @param outdir Character. Directory to save output plots (default = \code{"gsea_plots"}).
#' @param time_col Character. Column name specifying timepoints in \code{obj@meta.data}
#'   (default = \code{"timepoint"}).
#' @param time_include Character or numeric vector. Optional list of timepoints to include.
#'   If \code{NULL} (default), all timepoints are analyzed.
#'
#' @details
#' The function runs \code{Seurat::FindMarkers()} for each infection level compared to 
#' the reference level, then ranks genes for GSEA using either log2 fold-change 
#' or signed -log10 adjusted p-values. GSEA is performed using \pkg{clusterProfiler} 
#' with Hallmark gene sets. The function outputs a list of differential expression 
#' tables, GSEA results, and ggplot2 barplots of normalized enrichment scores (NES).
#'
#' When \code{time_include} is provided, only cells belonging to the specified 
#' timepoints are used for analysis. This allows users to focus on specific 
#' infection stages or experimental time windows.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{DE}: Differential expression results for each infection level.
#'   \item \code{GSEA}: GSEA results (\code{gseaResult} objects) for each level.
#'   \item \code{plots}: ggplot2 barplots of NES for up- and down-regulated pathways.
#' }
#'
#' @examples
#' \dontrun{
#' # Run analysis across all timepoints
#' res_all <- InfectionDynamic(
#'   obj = data,
#'   group_col = "infection_level",
#'   ref_level = 0,
#'   gmt_path = "h.all.v2025.1.Hs.symbols.gmt"
#' )
#'
#' # Restrict analysis to specific timepoints
#' res_tp <- InfectionDynamic(
#'   obj = data,
#'   group_col = "infection_level",
#'   ref_level = 0,
#'   time_col = "timepoint",
#'   time_include = c("24hpi", "48hpi"),
#'   save_plot = TRUE,
#'   outdir = "GSEA_TP_24_48"
#' )
#' }
#'
#' @importFrom Seurat FindMarkers Idents
#' @importFrom clusterProfiler GSEA read.gmt
#' @importFrom msigdbr msigdbr
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual coord_flip geom_hline geom_text labs theme_classic theme element_text ggsave
#' @importFrom dplyr select distinct arrange desc
#' @importFrom stats setNames
#' @export
InfectionDynamic <- function(
    obj,
    group_col = "infection_level",
    ref_level = 0,
    gmt_path = NULL,                 
    min.pct = 0.1,
    logfc.threshold = 0,
    ranks_from = c("avg_log2FC","pval_adj_signed"),
    top_up = 10,
    top_down = 10,
    col_up = "#D74C4C",
    col_down = "#3B6FB3",
    drop_prefix = TRUE,
    save_plot = FALSE,
    outdir = "gsea_plots",
    time_col = "timepoint",
    time_include = NULL
){
  ranks_from <- match.arg(ranks_from)

  # --- basic checks ---
  if (!inherits(obj, "Seurat")) stop("'obj' must be a Seurat object.")
  if (!group_col %in% colnames(obj[[]])) {
    stop(sprintf("Column '%s' not found in meta.data.", group_col))
  }
  if (!is.null(time_include)) {
    if (!time_col %in% colnames(obj[[]])) {
      stop(sprintf("You set time_include, but column '%s' was not found in meta.data.", time_col))
    }
    md <- obj[[]]
    ti_all  <- as.character(md[[time_col]])
    ti_want <- as.character(time_include)
    cells_use <- rownames(md)[ti_all %in% ti_want]
    if (length(cells_use) == 0) {
      stop(sprintf("No cells matched timepoints: %s (in column '%s').",
                   paste(unique(ti_want), collapse = ", "), time_col))
    }
    obj <- subset(obj, cells = cells_use)
  }

  # --- load Hallmark TERM2GENE ---
  if (is.null(gmt_path)) {
    gmt_path <- system.file(
      "h.all.v2025.1.Hs.symbols.gmt", 
      package = "SingleViral"
    )
  }

  if (!is.null(gmt_path) && nzchar(gmt_path) && file.exists(gmt_path)) {
    message("✅ Using GMT file from: ", gmt_path)
    hallmark <- clusterProfiler::read.gmt(gmt_path)
    colnames(hallmark) <- c("gs_name","gene_symbol")
  } else {
    message("⚠️ GMT file not found, falling back to msigdbr(H) ...")
    hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") |>
      dplyr::select(gs_name, gene_symbol) |>
      dplyr::distinct()
  }

  # --- levels to compare vs reference ---
  all_lvls <- sort(unique(obj[[group_col]][,1]))
  target_lvls <- setdiff(all_lvls, ref_level)
  if (length(target_lvls) == 0) stop("No other levels found besides the reference level.")

  if (save_plot && !dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  Seurat::Idents(obj) <- obj[[group_col]][,1]
  plots <- list(); de_list <- list(); gsea_list <- list()

  tp_tag <- if (is.null(time_include)) "" else {
    paste0("_TP_", paste(unique(as.character(time_include)), collapse = "-"))
  }

  for (lvl in target_lvls) {
    message(sprintf("\n==> Level %s vs %s | Running FindMarkers...", lvl, ref_level))
    de <- Seurat::FindMarkers(
      object          = obj,
      ident.1         = as.character(lvl),
      ident.2         = as.character(ref_level),
      min.pct         = min.pct,
      logfc.threshold = logfc.threshold,
      only.pos        = FALSE,
      test.use        = "wilcox"
    )
    de$gene <- rownames(de)

    if (ranks_from == "avg_log2FC") {
      ranks <- de$avg_log2FC; names(ranks) <- de$gene
    } else {
      s  <- sign(de$avg_log2FC)
      pa <- pmax(de$p_val_adj, 1e-300)
      ranks <- -log10(pa) * s; names(ranks) <- de$gene
    }
    ranks <- ranks[is.finite(ranks)]
    ranks <- ranks[!is.na(names(ranks)) & names(ranks) != ""]
    ranks <- sort(ranks, decreasing = TRUE)

    message("   Running GSEA (Hallmark)...")
    gsea_res <- tryCatch({
      clusterProfiler::GSEA(
        geneList      = ranks,
        TERM2GENE     = hallmark,
        verbose       = FALSE,
        pvalueCutoff  = 1,
        pAdjustMethod = "BH",
        minGSSize     = 10,
        eps           = 0
      )
    }, error = function(e){
      warning(paste("GSEA error at level", lvl, ":", e$message)); NULL
    })

    de_list[[as.character(lvl)]]   <- de
    gsea_list[[as.character(lvl)]] <- gsea_res

    if (!is.null(gsea_res) && inherits(gsea_res, "gseaResult")) {
      df <- tryCatch(as.data.frame(gsea_res), error = function(e) NULL)
      if (!is.null(df) && nrow(df) > 0) {
        if (!"Description" %in% colnames(df) && "ID" %in% colnames(df)) df$Description <- df$ID
        if (drop_prefix) df$Description <- gsub("^HALLMARK_", "", df$Description)

        up_df   <- df[df$NES > 0, , drop = FALSE] |> dplyr::arrange(dplyr::desc(NES)) |> head(top_up)
        down_df <- df[df$NES < 0, , drop = FALSE] |> dplyr::arrange(NES) |> head(top_down)

        plot_df <- rbind(
          transform(up_df,   direction = "Up"),
          transform(down_df, direction = "Down")
        )
        plot_df$Description <- factor(
          plot_df$Description,
          levels = rev(c(up_df$Description, down_df$Description))
        )

        p <- ggplot2::ggplot(
          plot_df,
          ggplot2::aes(x = Description, y = NES,
                       fill = factor(direction, levels = c("Up","Down")))
        ) +
          ggplot2::geom_col(width = 0.7) +
          ggplot2::scale_fill_manual(values = c(Up = col_up, Down = col_down)) +
          ggplot2::coord_flip() +
          ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
          ggplot2::geom_text(
            ggplot2::aes(label = round(NES, 2)),
            hjust = ifelse(plot_df$NES > 0, -0.1, 1.1),
            size = 3.2, fontface = "bold"
          ) +
          ggplot2::labs(
            title = paste0("Level ", lvl, " vs ", ref_level, " | Hallmark GSEA"),
            x = "Pathway",
            y = "Normalized Enrichment Score (NES)"
          ) +
          ggplot2::theme_classic(base_size = 12) +
          ggplot2::theme(
            plot.title   = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
            axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
            axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
            axis.text.x  = ggplot2::element_text(size = 12, face = "bold"),
            axis.text.y  = ggplot2::element_text(size = 12, face = "bold"),
            legend.text  = ggplot2::element_text(size = 12),
            legend.title = ggplot2::element_blank()
          )

        if (save_plot) {
          if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
          fn <- file.path(
            outdir,
            sprintf("GSEA_Level_%s_vs_%s_TopUp%d_TopDown%d%s.pdf",
                    lvl, ref_level, top_up, top_down, tp_tag)
          )
          ggplot2::ggsave(fn, p, width = 8, height = 6, useDingbats = FALSE)
        }
        plots[[as.character(lvl)]] <- p
      } else {
        plots[[as.character(lvl)]] <- NULL
        message(sprintf("⚠️ Level %s: empty or insufficient GSEA results.", lvl))
      }
    } else {
      plots[[as.character(lvl)]] <- NULL
      message(sprintf("⚠️ Level %s: GSEA object is NULL.", lvl))
    }
  }

  list(DE = de_list, GSEA = gsea_list, plots = plots)
}
