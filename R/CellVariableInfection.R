#' Summarize infection patterns and visualize viral load across cell groups
#'
#' @description
#' This function summarizes infection-related features from a Seurat object
#' and provides multiple publication-ready visualizations. Specifically, it:
#' \itemize{
#'   \item Computes the proportion of infected cells per group
#'   \item Computes the mean viral load per group
#'   \item (Optional) Visualizes stacked infection-level distributions
#'   \item (Optional) Generates UMAP visualizations split by cell types:
#'     \itemize{
#'       \item DimPlot colored by infection level
#'       \item FeaturePlot showing viral load distribution
#'     }
#' }
#'
#' Infection is defined as cells with viral load greater than \code{zero_tol}.
#'
#' @param obj A \strong{Seurat} object containing single-cell expression data
#'   and associated metadata.
#' @param group_col Character. Metadata column defining the main grouping variable
#'   (e.g. \code{"manual_label"}, \code{"timepoint"}).
#' @param score_col Character. Metadata column storing viral load scores
#'   (default: \code{"viral_load"}).
#' @param level_col Character. Optional metadata column defining discrete infection
#'   levels (e.g. 0/1/2). If not present, infection-level stacked barplots are skipped.
#' @param zero_tol Numeric. Cells with \code{score_col > zero_tol} are considered infected.
#' @param order_groups Optional character vector specifying the display order of groups.
#' @param palette Optional named character vector specifying colors for groups.
#' @param lvl_colors Optional named character vector specifying colors for infection levels.
#' @param save_dir Optional character. Output directory to save plots as PDF files.
#'   If \code{NULL}, plots are not saved.
#' @param show_stack_count_label Logical. Whether to display cell-count labels on
#'   stacked count barplots (default: \code{TRUE}).
#' @param show_stack_percent_label Logical. Whether to display percentage labels on
#'   stacked proportion barplots (default: \code{TRUE}).
#'
#' @param make_umap_level Logical. Whether to generate a UMAP \code{DimPlot}
#'   colored by \code{level_col} and split by \code{cell_type} (default: \code{TRUE}).
#' @param make_umap_viral Logical. Whether to generate a UMAP \code{FeaturePlot}
#'   showing \code{score_col} and split by \code{cell_type} (default: \code{TRUE}).
#' @param reduction Character. Dimensional reduction to use for UMAP plots
#'   (default: \code{"umap"}).
#' @param cell_type Character. Metadata column used for splitting UMAP panels
#'   (default: \code{"manual_label"}).
#' @param split_ncol Integer. Number of columns used when arranging split \code{DimPlot}
#'   panels (default: \code{6}).
#' @param viral_ncol Integer. Number of columns used when arranging split viral-load
#'   \code{FeaturePlot} panels (default: \code{3}).
#'
#' @details
#' The function performs all calculations using metadata stored in
#' \code{obj@meta.data}. Grouping variables are internally coerced to character
#' or factor types to avoid unintended behavior.
#'
#' When \code{make_umap_level = TRUE}, the function generates a split UMAP
#' \code{DimPlot} colored by infection levels.
#'
#' When \code{make_umap_viral = TRUE}, the function generates a split UMAP
#' \code{FeaturePlot} showing viral load values.
#'
#' All generated plots can be accessed from the returned list and optionally
#' saved to disk if \code{save_dir} is provided.
#'
#' @return
#' A named list with two elements:
#' \describe{
#'   \item{tables}{A list of data frames summarizing infection statistics,
#'     including infection proportions, mean viral load, and infection-level
#'     distributions (if available).}
#'   \item{plots}{A list of \code{ggplot} or \code{patchwork} objects, including:
#'     \itemize{
#'       \item infection proportion barplot
#'       \item mean viral load barplot
#'       \item stacked infection-level barplots (optional)
#'       \item UMAP split DimPlot for infection levels (optional)
#'       \item UMAP split FeaturePlot for viral load (optional)
#'     }}
#' }
#'
#' @import Seurat
#' @import ggplot2
#' @import dplyr
#' @importFrom scales percent_format
#' @importFrom patchwork wrap_plots
#'
#' @examples
#' \dontrun{
#' res <- CellVariableInfection(
#'   obj        = seurat_obj,
#'   group_col = "manual_label",
#'   lvl_colors = lvl_colors,
#'   save_dir   = "CellVariableInfection_plots"
#' )
#'
#' res$plots$umap_viral_load_split
#' }
#'
#' @export
CellVariableInfection <- function(
    obj,
    group_col   = "manual_label",
    score_col   = "viral_load",
    level_col   = "infection_level",
    zero_tol    = 0,
    order_groups = NULL,
    palette     = NULL,
    lvl_colors  = NULL,
    save_dir    = NULL,
    show_stack_count_label   = TRUE,
    show_stack_percent_label = TRUE,
    make_umap_level = TRUE,
    make_umap_viral = TRUE,
    reduction   = "umap",
    cell_type   = "manual_label",
    split_ncol  = 6,
    viral_ncol  = 3
){
  if (!inherits(obj, "Seurat")) stop("'obj' must be a Seurat object.")
  meta <- obj@meta.data
  
  ## ---- basic column checks ----
  for (cn in c(group_col, score_col)) {
    if (!cn %in% colnames(meta)) {
      stop(sprintf("Column '%s' not found in meta.data.", cn))
    }
  }
  has_level <- level_col %in% colnames(meta)
  
  ## ---- clean column types ----
  g_raw <- unlist(meta[[group_col]])
  s_raw <- unlist(meta[[score_col]])
  meta[[group_col]] <- as.character(g_raw)
  meta[[score_col]] <- suppressWarnings(as.numeric(s_raw))
  
  if (any(!is.finite(meta[[score_col]]))) {
    stop(sprintf("Non-numeric or NA values detected in '%s'.", score_col))
  }
  
  if (has_level) {
    l_raw <- unlist(meta[[level_col]])
    if (all(grepl("^\\d+$", as.character(l_raw)))) {
      meta[[level_col]] <- as.integer(l_raw)
    } else {
      meta[[level_col]] <- as.character(l_raw)
    }
  }
  
  ## ---- group order ----
  if (!is.null(order_groups)) {
    lvl_groups <- order_groups
  } else if (is.factor(obj@meta.data[[group_col]])) {
    lvl_groups <- levels(obj@meta.data[[group_col]])
  } else {
    lvl_groups <- unique(meta[[group_col]])
  }
  meta[[group_col]] <- factor(meta[[group_col]], levels = lvl_groups)
  
  ## =============== ① Infection ratio per group ===============
  infection_summary <- meta %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(
      total          = dplyr::n(),
      infected       = sum(.data[[score_col]] > zero_tol, na.rm = TRUE),
      infected_ratio = infected / total,
      .groups = "drop"
    ) %>%
    dplyr::rename(group = !!group_col) %>%
    dplyr::mutate(group = factor(group, levels = lvl_groups))
  
  p_prop <- ggplot2::ggplot(
    infection_summary,
    ggplot2::aes(x = group, y = infected_ratio, fill = group)
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", infected_ratio * 100)),
      vjust = -0.3, size = 3.2
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, NA)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::ggtitle("Proportion of infected cells within each group") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      plot.title  = ggplot2::element_text(hjust = 0.5, size = 20),
      legend.position = "none"
    ) +
    ggplot2::xlab("Group") + ggplot2::ylab("Infected cell proportion (%)")
  if (!is.null(palette)) p_prop <- p_prop + ggplot2::scale_fill_manual(values = palette)
  
  ## =============== ② Mean viral load per group ===============
  mean_viral <- meta %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(
      mean_viral_load = mean(.data[[score_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(group = !!group_col) %>%
    dplyr::mutate(group = factor(group, levels = lvl_groups))
  
  p_mean <- ggplot2::ggplot(
    mean_viral,
    ggplot2::aes(x = group, y = mean_viral_load, fill = group)
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", mean_viral_load)),
      vjust = -0.3, size = 3.2
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::ggtitle("Mean viral load per group") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      plot.title  = ggplot2::element_text(hjust = 0.5, size = 20),
      legend.position = "none"
    ) +
    ggplot2::xlab("Group") + ggplot2::ylab("Mean viral load")
  if (!is.null(palette)) p_mean <- p_mean + ggplot2::scale_fill_manual(values = palette)
  
  ## =============== ③ Infection level stacked barplots ===============
  if (has_level) {
    stacked_counts <- meta %>%
      dplyr::group_by(.data[[group_col]], .data[[level_col]]) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::rename(group = !!group_col, level = !!level_col) %>%
      dplyr::mutate(
        group   = factor(group, levels = lvl_groups),
        n_log   = log2(n + 1)
      )
    
    if (all(grepl("^\\d+$", as.character(stacked_counts$level)))) {
      stacked_counts$level <- factor(
        as.integer(as.character(stacked_counts$level)),
        levels = sort(unique(as.integer(as.character(stacked_counts$level))))
      )
    } else {
      stacked_counts$level <- factor(
        as.character(stacked_counts$level),
        levels = sort(unique(as.character(stacked_counts$level)))
      )
    }
    
    p_stack_count <- ggplot2::ggplot(
      stacked_counts,
      ggplot2::aes(x = group, y = n_log, fill = level)
    ) +
      ggplot2::geom_col(width = 0.8)
    
    if (show_stack_count_label) {
      p_stack_count <- p_stack_count +
        ggplot2::geom_text(
          ggplot2::aes(label = n),
          position = ggplot2::position_stack(vjust = 0.5),
          size = 3.0, color = "white"
        )
    }
    
    p_stack_count <- p_stack_count +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::ggtitle("Counts of infection levels within each group") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 14),
        plot.title  = ggplot2::element_text(hjust = 0.5, size = 20)
      ) +
      ggplot2::xlab("Group") +
      ggplot2::ylab("log2(Cell counts + 1)") +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Infection level"))
    if (!is.null(lvl_colors)) p_stack_count <- p_stack_count + ggplot2::scale_fill_manual(values = lvl_colors)
    
    stacked_percent <- stacked_counts %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(frac = n / sum(n)) %>%
      dplyr::ungroup()
    
    p_stack_percent <- ggplot2::ggplot(
      stacked_percent,
      ggplot2::aes(x = group, y = frac, fill = level)
    ) +
      ggplot2::geom_col(width = 0.8)
    
    if (show_stack_percent_label) {
      p_stack_percent <- p_stack_percent +
        ggplot2::geom_text(
          ggplot2::aes(label = scales::percent(frac, accuracy = 1)),
          position = ggplot2::position_stack(vjust = 0.5),
          size = 3.0, color = "white"
        )
    }
    
    p_stack_percent <- p_stack_percent +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        limits = c(0, 1)
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::ggtitle("Proportion of infection levels within each group") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 14),
        plot.title  = ggplot2::element_text(hjust = 0.5, size = 20)
      ) +
      ggplot2::xlab("Group") + ggplot2::ylab("Proportion") +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Infection level"))
    if (!is.null(lvl_colors)) p_stack_percent <- p_stack_percent + ggplot2::scale_fill_manual(values = lvl_colors)
    
  } else {
    stacked_counts   <- NULL
    stacked_percent  <- NULL
    p_stack_count    <- NULL
    p_stack_percent  <- NULL
  }
  
  ## =============== ④ EXTRA: UMAP plots ===============
  p_umap_level <- NULL
  p_umap_viral <- NULL
  
  # ④.1 DimPlot: infection_level split by manual_label
  if (make_umap_level) {
    if (!has_level) {
      warning("make_umap_level = TRUE, but level_col '", level_col, "' not found. Skip DimPlot.")
    } else {
      if (!cell_type %in% colnames(obj[[]])) {
        warning("make_umap_level = TRUE, but cell_type column '", cell_type, "' not found. Skip DimPlot.")
      } else {
        p_umap_level <- Seurat::DimPlot(
          obj,
          reduction = reduction,
          group.by  = level_col,
          label     = FALSE,
          split.by  = cell_type,
          cols      = lvl_colors,
          ncol      = split_ncol
        ) + ggplot2::ggtitle("Infection level")
      }
    }
  }
  
  # ④.2 FeaturePlot: viral_load split by manual_label, wrap into 3 columns
  if (make_umap_viral) {
    if (!cell_type %in% colnames(obj[[]])) {
      warning("make_umap_viral = TRUE, but cell_type column '", cell_type, "' not found. Skip FeaturePlot.")
    } else {
      p_list <- Seurat::FeaturePlot(
        obj,
        features = score_col,
        cols     = c("grey", "red"),
        split.by = cell_type,
        combine  = FALSE
      )
      p_list <- lapply(p_list, function(p){
        p + ggplot2::theme(
          axis.title.x = ggplot2::element_text(),
          axis.text.x  = ggplot2::element_text(),
          axis.ticks.x = ggplot2::element_line()
        )
      })
      p_umap_viral <- patchwork::wrap_plots(p_list, ncol = viral_ncol)
    }
  }
  
  ## =============== Save outputs ===============
  if (!is.null(save_dir)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    
    ggplot2::ggsave(file.path(save_dir, "Infected_proportion_by_group.pdf"),
                    p_prop, width = 8, height = 5, useDingbats = FALSE)
    ggplot2::ggsave(file.path(save_dir, "Mean_viral_load_by_group.pdf"),
                    p_mean, width = 8, height = 5, useDingbats = FALSE)
    
    if (!is.null(p_stack_count)) {
      ggplot2::ggsave(file.path(save_dir, "Infection_levels_stacked_counts.pdf"),
                      p_stack_count, width = 8, height = 5, useDingbats = FALSE)
      ggplot2::ggsave(file.path(save_dir, "Infection_levels_stacked_percent.pdf"),
                      p_stack_percent, width = 8, height = 5, useDingbats = FALSE)
    }
    
    if (!is.null(p_umap_level)) {
      ggplot2::ggsave(file.path(save_dir, "UMAP_split_infection_level.pdf"),
                      p_umap_level, width = 14, height = 6, useDingbats = FALSE)
    }
    if (!is.null(p_umap_viral)) {
      ggplot2::ggsave(file.path(save_dir, "UMAP_split_viral_load.pdf"),
                      p_umap_viral, width = 14, height = 10, useDingbats = FALSE)
    }
  }
  
  list(
    tables = list(
      infection_summary        = infection_summary,
      mean_viral_load          = mean_viral,
      infection_level_counts   = stacked_counts,
      infection_level_percent  = if (has_level) stacked_percent else NULL
    ),
    plots = list(
      infection_proportion         = p_prop,
      mean_viral_load              = p_mean,
      infection_level_counts_plot  = p_stack_count,
      infection_level_percent_plot = p_stack_percent,
      umap_infection_level_split   = p_umap_level,
      umap_viral_load_split        = p_umap_viral
    )
  )
}