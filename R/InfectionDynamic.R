#' Perform Infection-Level GSEA Analysis with Species-Specific Hallmark GMT Support
#'
#' This function performs differential expression (DE) and Hallmark GSEA analysis
#' between infection levels within a Seurat object. It supports filtering by
#' selected infection levels and/or specific groups (e.g. timepoints) and/or cell
#' types, ranks genes by \code{avg_log2FC} or signed adjusted p-values, and produces
#' publication-ready NES barplots for each comparison.
#'
#' For \code{species = "human"} or \code{"mouse"}, the function will automatically
#' use species-specific Hallmark GMT files bundled in the \pkg{SingleViral} package.
#' For \code{species = "other"}, you must provide a valid \code{gmt_path}.
#'
#' Robustness: when a comparison group has too few cells (e.g. <3), \code{Seurat::FindMarkers}
#' will error. This function will \strong{skip} such levels and continue running,
#' returning \code{NULL} entries for the failed level while keeping all successful results.
#'
#' @param obj A \code{Seurat} object containing expression data and metadata.
#' @param level_col Character. Column name in \code{obj@meta.data} specifying infection levels
#'   (e.g. \code{"infection_level"}).
#' @param ref_level Character or numeric. Reference level to compare against.
#'   Must be present after optional filtering.
#' @param level_include Character or numeric vector. Optional list of infection levels to include.
#'   If provided, only cells whose \code{level_col} is in \code{level_include} are kept.
#'   If \code{NULL} (default), all infection levels are used.
#' @param species Character. Species used for GSEA, one of \code{"human"}, \code{"mouse"}, or \code{"other"}.
#' @param gmt_path Optional character path to a GMT file. If \code{NULL}, built-in GMT files are used
#'   for \code{species="human"} or \code{"mouse"}. For \code{species="other"}, this is required.
#' @param min.pct Numeric. Minimum fraction of cells expressing a gene in either group
#'   (passed to \code{Seurat::FindMarkers}).
#' @param logfc.threshold Numeric. Minimum log2 fold-change threshold for DE testing (default = 0).
#' @param min_cells_per_group Integer. Minimum number of cells required in each comparison group
#'   for running DE. Levels that do not meet this threshold will be skipped (default = 3).
#' @param ranks_from Character. Ranking method for GSEA; one of \code{"avg_log2FC"} or
#'   \code{"pval_adj_signed"} (signed -log10 adjusted p-values).
#' @param top_up Integer. Number of top upregulated pathways (NES > 0) to display.
#' @param top_down Integer. Number of top downregulated pathways (NES < 0) to display.
#' @param col_up Character. Fill color for upregulated pathways (default = \code{"#D74C4C"}).
#' @param col_down Character. Fill color for downregulated pathways (default = \code{"#3B6FB3"}).
#' @param drop_prefix Logical. Whether to remove the \code{"HALLMARK_"} prefix from pathway names.
#' @param save_plot Logical. Whether to save plots as PDF files.
#' @param outdir Character. Output directory for plots (default = \code{"gsea_plots"}).
#' @param group_col Character. Column name specifying groups in \code{obj@meta.data}
#'   (default = \code{"timepoint"}).
#' @param group_include Character or numeric vector. Optional list of groups to include.
#' @param cell_type Character. Column name specifying cell type (default = \code{"manual_label"}).
#' @param celltype_include Character or numeric vector. Optional list of cell types to include.
#' @param title_prefix Character. Optional title prefix added to plots.
#'
#' @details
#' Filtering order is: \code{level_include} -> \code{group_include} -> \code{celltype_include}.
#' When multiple filters are provided, the intersection of the selected cells is used.
#'
#' For \code{species="human"} and \code{"mouse"}, Hallmark GMT files are located via
#' \code{system.file()} from the \pkg{SingleViral} package. If the GMT cannot be found,
#' the function falls back to \pkg{msigdbr} Hallmark gene sets.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{DE}: Differential expression results (data.frames) for each infection level;
#'         levels that were skipped will be \code{NULL}.
#'   \item \code{GSEA}: GSEA results (\code{gseaResult} objects) for each level; may be \code{NULL}.
#'   \item \code{plots}: ggplot2 barplots of NES for each level; may be \code{NULL}.
#' }
#'
#' @examples
#' \dontrun{
#' res <- InfectionDynamic(
#'   obj            = data,
#'   level_col      = "infection_level",
#'   ref_level      = 0,
#'   level_include  = c(0,1,2,3,4),
#'   species        = "human",
#'   group_col      = "timepoint",
#'   group_include  = c("24hpi","48hpi"),
#'   cell_type      = "manual_label",
#'   celltype_include = c("Lymphoid cells"),
#'   min_cells_per_group = 3,
#'   save_plot      = TRUE,
#'   outdir         = "gsea_plots_subset"
#' )
#' }
#'
#' @importFrom Seurat FindMarkers Idents
#' @importFrom clusterProfiler GSEA read.gmt
#' @importFrom msigdbr msigdbr
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual coord_flip geom_hline geom_text labs
#'   theme_classic theme element_text ggsave
#' @importFrom dplyr select distinct arrange desc
#' @export
InfectionDynamic <- function(
    obj,
    level_col        = "infection_level",
    ref_level        = 0,
    level_include    = NULL,
    species          = c("human", "mouse", "other"),
    gmt_path         = NULL,
    min.pct          = 0.1,
    logfc.threshold  = 0,
    min_cells_per_group = 3,
    ranks_from       = c("avg_log2FC","pval_adj_signed"),
    top_up           = 10,
    top_down         = 10,
    col_up           = "#D74C4C",
    col_down         = "#3B6FB3",
    drop_prefix      = TRUE,
    save_plot        = FALSE,
    outdir           = "gsea_plots",
    group_col        = "timepoint",
    group_include    = NULL,
    cell_type        = "manual_label",
    celltype_include = NULL,
    title_prefix     = ""
){
  ranks_from <- match.arg(ranks_from)
  species    <- match.arg(species)
  
  # --- basic checks ---
  if (!inherits(obj, "Seurat")) stop("'obj' must be a Seurat object.")
  if (!level_col %in% colnames(obj[[]])) {
    stop(sprintf("Column '%s' not found in meta.data.", level_col))
  }
  
  # --- optional infection-level filtering ---
  if (!is.null(level_include)) {
    md0 <- obj[[]]
    lvl_all  <- as.character(md0[[level_col]])
    lvl_want <- as.character(level_include)
    
    cells_use_lvl <- rownames(md0)[lvl_all %in% lvl_want]
    if (length(cells_use_lvl) == 0) {
      stop(sprintf("No cells matched levels: %s (in column '%s').",
                   paste(unique(lvl_want), collapse = ", "), level_col))
    }
    obj <- subset(obj, cells = cells_use_lvl)
    
    if (!as.character(ref_level) %in% as.character(unique(obj[[level_col]][, 1]))) {
      stop(sprintf("ref_level '%s' is not present after level_include filtering.", as.character(ref_level)))
    }
  }
  
  # --- optional group/timepoint filtering ---
  if (!is.null(group_include)) {
    if (!group_col %in% colnames(obj[[]])) {
      stop(sprintf("You set group_include, but column '%s' was not found in meta.data.", group_col))
    }
    md <- obj[[]]
    ti_all  <- as.character(md[[group_col]])
    ti_want <- as.character(group_include)
    cells_use_tp <- rownames(md)[ti_all %in% ti_want]
    if (length(cells_use_tp) == 0) {
      stop(sprintf("No cells matched groups: %s (in column '%s').",
                   paste(unique(ti_want), collapse = ", "), group_col))
    }
    obj <- subset(obj, cells = cells_use_tp)
  }
  
  # --- optional cell-type filtering ---
  if (!is.null(celltype_include)) {
    if (!cell_type %in% colnames(obj[[]])) {
      stop(sprintf("You set celltype_include, but column '%s' was not found in meta.data.", cell_type))
    }
    md_ct <- obj[[]]
    ct_all  <- as.character(md_ct[[cell_type]])
    ct_want <- as.character(celltype_include)
    cells_use_ct <- rownames(md_ct)[ct_all %in% ct_want]
    if (length(cells_use_ct) == 0) {
      stop(sprintf("No cells matched cell types: %s (in column '%s').",
                   paste(unique(ct_want), collapse = ", "), cell_type))
    }
    obj <- subset(obj, cells = cells_use_ct)
  }
  
  # --- choose GMT path based on species / gmt_path ---
  if (is.null(gmt_path) || !nzchar(gmt_path)) {
    if (species == "human") {
      gmt_path <- system.file("h.all.v2025.1.Hs.symbols.gmt", package = "SingleViral")
    } else if (species == "mouse") {
      gmt_path <- system.file("mh.all.v2025.1.Mm.symbols.gmt", package = "SingleViral")
    } else {
      stop("For species = 'other', you must provide a valid 'gmt_path'.")
    }
  }
  
  # --- load Hallmark TERM2GENE ---
  hallmark <- NULL
  if (!is.null(gmt_path) && nzchar(gmt_path) && file.exists(gmt_path)) {
    message("✅ Using GMT file from: ", gmt_path)
    hallmark <- clusterProfiler::read.gmt(gmt_path)
    colnames(hallmark) <- c("gs_name","gene_symbol")
  } else {
    if (species %in% c("human","mouse")) {
      sp_name <- if (species == "human") "Homo sapiens" else "Mus musculus"
      message("⚠️ GMT not found at '", gmt_path, "'. Falling back to msigdbr Hallmark for ", sp_name, " ...")
      hallmark <- msigdbr::msigdbr(species = sp_name, category = "H") |>
        dplyr::select(gs_name, gene_symbol) |>
        dplyr::distinct()
    } else {
      stop("GMT file not found at '", gmt_path, "'. Please provide a valid GMT file.")
    }
  }
  
  # --- levels to compare vs reference ---
  ref_chr   <- as.character(ref_level)
  all_lvls  <- sort(unique(as.character(obj[[level_col]][, 1])))
  target_lvls <- setdiff(all_lvls, ref_chr)
  if (length(target_lvls) == 0) stop("No other levels found besides the reference level.")
  
  if (save_plot && !dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # IMPORTANT: use character idents to avoid numeric/continuous surprises
  Seurat::Idents(obj) <- as.character(obj[[level_col]][, 1])
  
  plots <- list()
  de_list <- list()
  gsea_list <- list()
  
  # build filename tag
  tp_tag <- ""
  if (!is.null(level_include)) {
    tp_tag <- paste0(tp_tag, "_LV_", paste(unique(as.character(level_include)), collapse = "-"))
  }
  if (!is.null(group_include)) {
    tp_tag <- paste0(tp_tag, "_TP_", paste(unique(as.character(group_include)), collapse = "-"))
  }
  if (!is.null(celltype_include)) {
    tp_tag <- paste0(tp_tag, "_CT_", paste(unique(as.character(celltype_include)), collapse = "-"))
  }
  
  for (lvl in target_lvls) {
    
    # ---- cell count check (skip instead of error) ----
    md_lvl <- obj[[]][, level_col, drop = TRUE]
    n_lvl  <- sum(as.character(md_lvl) == as.character(lvl))
    n_ref  <- sum(as.character(md_lvl) == ref_chr)
    
    if (n_lvl < min_cells_per_group || n_ref < min_cells_per_group) {
      warning(sprintf(
        "⚠️ Skipping Level %s vs %s: insufficient cells (lvl=%d, ref=%d, min=%d).",
        as.character(lvl), ref_chr, n_lvl, n_ref, min_cells_per_group
      ))
      de_list[[as.character(lvl)]]   <- NULL
      gsea_list[[as.character(lvl)]] <- NULL
      plots[[as.character(lvl)]]     <- NULL
      next
    }
    
    message(sprintf("\n==> Level %s vs %s | Running FindMarkers...", as.character(lvl), ref_chr))
    
    de <- tryCatch({
      Seurat::FindMarkers(
        object          = obj,
        ident.1         = as.character(lvl),
        ident.2         = ref_chr,
        min.pct         = min.pct,
        logfc.threshold = logfc.threshold,
        only.pos        = FALSE,
        test.use        = "wilcox"
      )
    }, error = function(e) {
      warning(sprintf("❌ FindMarkers failed at Level %s vs %s: %s",
                      as.character(lvl), ref_chr, e$message))
      NULL
    })
    
    if (is.null(de) || nrow(de) == 0) {
      warning(sprintf("⚠️ Level %s vs %s: DE result is empty. Skipping GSEA/plot.",
                      as.character(lvl), ref_chr))
      de_list[[as.character(lvl)]]   <- de
      gsea_list[[as.character(lvl)]] <- NULL
      plots[[as.character(lvl)]]     <- NULL
      next
    }
    
    de$gene <- rownames(de)
    
    # --- build ranks ---
    if (ranks_from == "avg_log2FC") {
      ranks <- de$avg_log2FC
      names(ranks) <- de$gene
    } else {
      s  <- sign(de$avg_log2FC)
      pa <- pmax(de$p_val_adj, 1e-300)
      ranks <- -log10(pa) * s
      names(ranks) <- de$gene
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
    }, error = function(e) {
      warning(paste("GSEA error at level", as.character(lvl), ":", e$message))
      NULL
    })
    
    de_list[[as.character(lvl)]]   <- de
    gsea_list[[as.character(lvl)]] <- gsea_res
    
    # --- make plot ---
    if (!is.null(gsea_res) && inherits(gsea_res, "gseaResult")) {
      df <- tryCatch(as.data.frame(gsea_res), error = function(e) NULL)
      
      if (!is.null(df) && nrow(df) > 0) {
        if (!"Description" %in% colnames(df) && "ID" %in% colnames(df)) {
          df$Description <- df$ID
        }
        if (drop_prefix) df$Description <- gsub("^HALLMARK_", "", df$Description)
        
        # 去掉 NA NES
        df <- df[is.finite(df$NES), , drop = FALSE]
        
        up_df <- df[df$NES > 0, , drop = FALSE] |>
          dplyr::arrange(dplyr::desc(NES)) |>
          utils::head(top_up)
        
        down_df <- df[df$NES < 0, , drop = FALSE] |>
          dplyr::arrange(NES) |>
          utils::head(top_down)
        
        plot_df_list <- list()
        if (nrow(up_df) > 0)   plot_df_list <- c(plot_df_list, list(transform(up_df, direction = "Up")))
        if (nrow(down_df) > 0) plot_df_list <- c(plot_df_list, list(transform(down_df, direction = "Down")))
        
        if (length(plot_df_list) == 0) {
          warning(sprintf("⚠️ Level %s: no pathways to plot after filtering (NES all NA or no +/−).", lvl))
          plots[[as.character(lvl)]] <- NULL
          next
        }
        plot_df <- do.call(rbind, plot_df_list)
        
        plot_df$Description <- factor(
          plot_df$Description,
          levels = rev(c(up_df$Description, down_df$Description))
        )
        
        p <- ggplot2::ggplot(
          plot_df,
          ggplot2::aes(
            x = Description,
            y = NES,
            fill = factor(direction, levels = c("Up", "Down"))
          )
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
            title = paste0(title_prefix, " Level ", as.character(lvl), " vs ", ref_chr, " | Hallmark GSEA"),
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
                    as.character(lvl), ref_chr, top_up, top_down, tp_tag)
          )
          ggplot2::ggsave(fn, p, width = 8, height = 6, useDingbats = FALSE)
        }
        
        plots[[as.character(lvl)]] <- p
      } else {
        plots[[as.character(lvl)]] <- NULL
        message(sprintf("⚠️ Level %s: empty or insufficient GSEA results.", as.character(lvl)))
      }
    } else {
      plots[[as.character(lvl)]] <- NULL
      message(sprintf("⚠️ Level %s: GSEA object is NULL.", as.character(lvl)))
    }
  }
  
  list(DE = de_list, GSEA = gsea_list, plots = plots)
}