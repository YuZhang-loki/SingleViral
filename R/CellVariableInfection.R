#' Cell-type-wise infection summary and visualization for Seurat objects
#'
#' This function summarizes infection/viral-load related metrics across cell groups
#' (e.g., annotated cell types) and generates a set of publication-ready plots,
#' including infected proportion, mean viral load, viral-load violin with infection
#' level coloring (plus mean point and per-group infection-level pie), stacked bars
#' of infection levels, and optional UMAP visualizations.
#'
#' @param obj A \code{\link[Seurat]{Seurat}} object.
#' @param group_col Character. Column name in \code{obj@meta.data} used as the group
#'   (e.g., cell type annotation). Default: \code{"manual_label"}.
#' @param score_col Character. Column name in \code{obj@meta.data} used as the viral
#'   load / infection score. Must be numeric (or coercible to numeric). Default: \code{"viral_load"}.
#' @param level_col Character. Column name in \code{obj@meta.data} used as infection
#'   level categories (e.g., 0–4). If absent, level-dependent plots are skipped.
#'   Default: \code{"infection_level"}.
#' Mean viral load (bar chart and violin mean points) excludes infection level 0.
#' Groups with no nonzero-level cells have NA mean.
#' @param zero_tol Numeric. Threshold above which a cell is considered infected
#'   when computing infected proportion. Default: \code{0}.
#' @param order_groups Optional character vector specifying the display/order of groups.
#'   If provided, the final group levels will be intersected with groups present after filtering.
#'   Default: \code{NULL}.
#' @param celltype_include Optional character vector. If provided, only cells whose
#'   \code{group_col} value is in \code{celltype_include} will be kept for all downstream
#'   analyses and plots. Default: \code{NULL} (use all groups).
#' @param palette Optional named vector. Manual fill colors for groups in barplots
#'   (\code{p_prop} and \code{p_mean}). Names should match group labels. Default: \code{NULL}.
#' @param lvl_colors Optional named vector. Manual colors for infection levels used in points,
#'   pies, and stacked bars. Names should match level labels (e.g., \code{"0"}, \code{"1"}, ...).
#'   Default: \code{NULL}.
#' @param save_dir Optional character. If not \code{NULL}, figures are saved as PDF files
#'   into this directory (created if needed). Default: \code{NULL}.
#' @param show_stack_count_label Logical. Whether to show count labels on stacked count bars.
#'   Default: \code{TRUE}.
#' @param show_stack_percent_label Logical. Whether to show percent labels on stacked percent bars.
#'   Default: \code{TRUE}.
#' @param make_umap_level Logical. Whether to generate split UMAP colored by infection level.
#'   Default: \code{TRUE}.
#' @param make_umap_viral Logical. Whether to generate split UMAP of viral load (FeaturePlot).
#'   Default: \code{TRUE}.
#' @param reduction Character. Dimensional reduction name used for UMAP plots. Default: \code{"umap"}.
#' @param cell_type Character. Column name used to split UMAP plots (via \code{split.by}).
#'   Default: \code{"manual_label"}.
#' @param split_ncol Integer. Number of columns for split UMAP plots. Default: \code{6}.
#' @param viral_ncol Integer. Number of columns when wrapping FeaturePlot panels. Default: \code{3}.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{tables}: a list of summary tables (\code{infection_summary}, \code{mean_viral_load},
#'         \code{infection_level_counts}, \code{infection_level_percent}).
#'   \item \code{plots}: a list of ggplot/patchwork objects.
#' }
#'
#' @examples
#' \dontrun{
#' res <- CellVariableInfection(
#'   obj,
#'   group_col = "manual_label",
#'   score_col = "viral_load",
#'   level_col = "infection_level",
#'   celltype_include = c("CD4+Tcells","CD8+Tcells","Tfh","Th17"),
#'   lvl_colors = c("0"="#7FB3D5","1"="turquoise","2"="#FFFF99","3"="#fdb462","4"="#D74C4C"),
#'   save_dir = "infection_summary_plots",
#'   show_stack_count_label = FALSE,
#'   show_stack_percent_label = FALSE
#' )
#' res$plots$mean_viral_load_violin
#' }
#'
#' @importFrom Seurat DimPlot FeaturePlot
#' @importFrom ggplot2 ggplot aes geom_col geom_text scale_y_continuous theme_classic ggtitle theme xlab ylab
#' @importFrom ggplot2 geom_violin geom_point position_jitter geom_polygon geom_blank coord_cartesian
#' @importFrom ggplot2 scale_fill_manual scale_shape_manual guides guide_legend element_text expansion
#' @importFrom dplyr group_by summarise n rename mutate transmute filter count ungroup arrange lag cumsum do rowwise
#' @importFrom scales percent_format percent
#' @importFrom patchwork wrap_plots
#'
#' @export
CellVariableInfection <- function(
    obj,
    group_col   = "manual_label",
    score_col   = "viral_load",
    level_col   = "infection_level",
    zero_tol    = 0,
    order_groups = NULL,
    celltype_include = NULL,   
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
  
  ## ===================== ✅ NEW: filter selected cell types =====================
  if (!is.null(celltype_include)) {
    celltype_include <- as.character(celltype_include)
    present_types <- sort(unique(meta[[group_col]]))
    missing_types <- setdiff(celltype_include, present_types)
    if (length(missing_types) > 0) {
      warning("celltype_include contains types not found in '", group_col, "': ",
              paste(missing_types, collapse = ", "))
    }
    
    keep_cells <- rownames(meta)[meta[[group_col]] %in% celltype_include]
    if (length(keep_cells) == 0) {
      stop("After applying celltype_include, 0 cells remain. Check your type names in ", group_col, ".")
    }
    
    
    obj  <- obj[, keep_cells, drop = FALSE]
    meta <- obj@meta.data
    meta[[group_col]] <- as.character(meta[[group_col]])
    meta[[score_col]] <- suppressWarnings(as.numeric(meta[[score_col]]))
  }
  ## =========================================================================
  
  ## ---- group order ----
  if (!is.null(order_groups)) {
    lvl_groups <- as.character(order_groups)
  } else if (is.factor(obj@meta.data[[group_col]])) {
    lvl_groups <- levels(obj@meta.data[[group_col]])
  } else {
    lvl_groups <- unique(meta[[group_col]])
  }
  
  
  lvl_groups <- lvl_groups[lvl_groups %in% unique(meta[[group_col]])]
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
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::ggtitle("Proportion of infected cells within each group") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      plot.title  = ggplot2::element_text(hjust = 0.5, size = 20),
      legend.position = "none"
    ) +
    ggplot2::xlab("Group") + ggplot2::ylab("Infected cell proportion (%)")
  if (!is.null(palette)) p_prop <- p_prop + ggplot2::scale_fill_manual(values = palette)
  
  ## =============== ② Mean viral load per group (exclude level 0) ===============
  mean_meta <- if (has_level) {
    meta %>% dplyr::filter(!is.na(.data[[level_col]]),
                           as.character(.data[[level_col]]) != "0")
  } else {
    warning("Infection level column missing; using viral load > zero_tol for mean calculation.")
    meta %>% dplyr::filter(.data[[score_col]] > zero_tol)
  }
  mean_viral <- meta %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(
      mean_viral_load = {
        v <- mean_meta[[score_col]][as.character(mean_meta[[group_col]]) ==
                                      as.character(dplyr::first(.data[[group_col]]))]
        if (length(v)) mean(v, na.rm = TRUE) else NA_real_
      },
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
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::ggtitle("Mean viral load per group") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      plot.title  = ggplot2::element_text(hjust = 0.5, size = 20),
      legend.position = "none"
    ) +
    ggplot2::xlab("Group") + ggplot2::ylab("Mean viral load")
  if (!is.null(palette)) p_mean <- p_mean + ggplot2::scale_fill_manual(values = palette)
  
  ## =============== ⑤ Viral load violin per group (points colored by infection level) ===============
  p_mean_violin <- NULL
  {
    df_violin <- meta %>%
      dplyr::transmute(
        group      = .data[[group_col]],
        viral_load = .data[[score_col]],
        level      = if (has_level) .data[[level_col]] else NA
      ) %>%
      dplyr::mutate(group = factor(group, levels = lvl_groups))
    
    if (has_level) {
      if (!is.null(lvl_colors)) {
        lvl_keys <- names(lvl_colors)
        df_violin$level <- factor(as.character(df_violin$level), levels = lvl_keys)
      } else {
        df_violin$level <- factor(as.character(df_violin$level),
                                  levels = sort(unique(as.character(df_violin$level))))
      }
    }
    
    
    df_point <- df_violin %>%
      dplyr::filter(is.finite(viral_load), viral_load >= 0)
    
    ## ---- mean point data (per group) ----
    mean_pt_df <- df_violin %>%
      dplyr::filter(if (has_level) as.character(level) != "0" & !is.na(level)
                    else viral_load > zero_tol) %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(mean_v = mean(viral_load, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(group = factor(group, levels = lvl_groups))
    
    ## ---- build pies: infection level proportions, ALL at the same height ----
    pie_df <- NULL
    y0_common <- NULL
    ry <- NULL
    if (has_level) {
      pie_prop <- df_violin %>%
        dplyr::filter(!is.na(level)) %>%
        dplyr::count(group, level, name = "n") %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(frac = n / sum(n)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(group, level) %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(
          theta0 = 2 * pi * dplyr::lag(cumsum(frac), default = 0),
          theta1 = 2 * pi * cumsum(frac)
        ) %>%
        dplyr::ungroup()
      
      y_rng  <- range(df_violin$viral_load, na.rm = TRUE)
      y_span <- diff(y_rng) + 1e-9
      
      ry <- 0.10 * y_span
      rx <- 0.22
      y0_common <- max(df_violin$viral_load, na.rm = TRUE) + 0.22 * y_span + ry * 1.3
      
      pie_prop <- pie_prop %>%
        dplyr::mutate(
          y0 = y0_common,
          x0 = as.numeric(factor(group, levels = lvl_groups))
        )
      
      n_arc <- 60
      pie_df <- pie_prop %>%
        dplyr::rowwise() %>%
        dplyr::do({
          th <- seq(.$theta0, .$theta1, length.out = n_arc)
          data.frame(
            group = .$group,
            level = .$level,
            x = c(.$x0, .$x0 + rx * cos(th), .$x0),
            y = c(.$y0, .$y0 + ry * sin(th), .$y0)
          )
        }) %>%
        dplyr::ungroup()
    }
    
    ## ---- y upper (avoid clipping) ----
    y_upper <- max(df_violin$viral_load, na.rm = TRUE)
    if (!is.null(y0_common) && !is.null(ry)) y_upper <- max(y_upper, y0_common + ry * 2)
    
    p_mean_violin <- ggplot2::ggplot(df_violin, ggplot2::aes(x = group, y = viral_load)) +
      ggplot2::geom_violin(
        ggplot2::aes(group = group),
        fill  = "grey85",
        color = NA,
        alpha = 0.7,
        trim  = TRUE,
        scale = "width"
      ) +
      (if (has_level) {
        ggplot2::geom_point(
          data = df_point,
          ggplot2::aes(x = group, y = viral_load, fill = level),
          position = ggplot2::position_jitter(width = 0.18, height = 0),
          shape = 21,
          color = "grey30",
          stroke = 0.15,
          size  = 1.15,
          alpha = 0.60
        )
      } else {
        ggplot2::geom_point(
          data = df_point,
          ggplot2::aes(x = group, y = viral_load),
          position = ggplot2::position_jitter(width = 0.18, height = 0),
          size  = 1.15,
          alpha = 0.60,
          color = "grey30"
        )
      }) +
      ggplot2::geom_point(
        data = mean_pt_df,
        ggplot2::aes(x = group, y = mean_v, shape = "Mean viral load"),
        inherit.aes = FALSE,
        color = "red",
        alpha = 0.75,
        size = 2.6
      ) +
      (if (has_level && !is.null(pie_df)) {
        ggplot2::geom_polygon(
          data = pie_df,
          ggplot2::aes(x = x, y = y, group = interaction(group, level), fill = level),
          inherit.aes = FALSE,
          alpha = 0.55,
          color = NA
        )
      } else {
        ggplot2::geom_blank()
      }) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::ggtitle("Analysis of cellular susceptibility to infection") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = ggplot2::element_text(size = 14),
        plot.title  = ggplot2::element_text(hjust = 0.5, size = 18)
      ) +
      ggplot2::xlab("Group") + ggplot2::ylab("Viral load") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
      ggplot2::coord_cartesian(ylim = c(0, y_upper), clip = "off") +
      ggplot2::scale_shape_manual(values = c("Mean viral load" = 16))
    
    if (has_level && !is.null(lvl_colors)) {
      p_mean_violin <- p_mean_violin +
        ggplot2::scale_fill_manual(values = lvl_colors, na.translate = FALSE) +
        ggplot2::guides(
          fill = ggplot2::guide_legend(
            title = "Infection level",
            override.aes = list(shape = 22, size = 5, alpha = 1, color = NA)
          ),
          shape = ggplot2::guide_legend(
            title = NULL,
            override.aes = list(color = "red", size = 3, alpha = 1)
          )
        )
    } else {
      p_mean_violin <- p_mean_violin + ggplot2::guides(fill = "none", shape = "none")
    }
  }
  
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
      ggplot2::theme_classic(base_size = 12) +
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
      ggplot2::theme_classic(base_size = 12) +
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
    
    if (!is.null(p_mean_violin)) {
      ggplot2::ggsave(file.path(save_dir, "Viral_load_violin_by_group_colored_by_level.pdf"),
                      p_mean_violin, width = 10, height = 5.5, useDingbats = FALSE)
    }
    
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
      mean_viral_load_violin       = p_mean_violin,
      infection_level_counts_plot  = p_stack_count,
      infection_level_percent_plot = p_stack_percent,
      umap_infection_level_split   = p_umap_level,
      umap_viral_load_split        = p_umap_viral
    )
  )
}
