#' Summarize infection across groups and visualize proportions/means/stacked levels
#'
#' @description
#' Given a Seurat object with viral-load and (optionally) infection levels,
#' this function (i) computes per-group infection proportion (viral_load > zero_tol),
#' (ii) computes per-group mean viral load, and (iii) if an infection level column
#' is present, draws stacked barplots for counts and percentages of levels per group.
#' Plots can be saved to disk. Colors and group order are customizable.
#'
#' @param obj A \strong{Seurat} object.
#' @param group_col Character. Column in \code{obj@meta.data} defining groups
#'   (e.g., \code{"manual_label"} or \code{"timepoint"}).
#' @param score_col Character. Column in \code{obj@meta.data} storing viral-load scores.
#' @param level_col Character. Optional column storing infection levels (0/1/2/…).
#'   If absent, stacked level plots are skipped.
#' @param zero_tol Numeric. Cells with \code{viral_load > zero_tol} are considered infected.
#' @param order_groups Optional character vector giving an explicit display order of groups.
#'   If \code{NULL}, uses factor levels in \code{obj@meta.data[[group_col]]} when available,
#'   otherwise the unique order found in the data.
#' @param palette Optional named character vector of colors for \code{group_col} bars.
#'   Names should match group names. If \code{NULL}, ggplot default palette is used.
#' @param lvl_colors Optional named character vector of colors for \code{level_col} (stacked bars).
#'   Names should match infection-level labels (or their numeric values as strings).
#' @param save_dir Optional output directory to save PDF plots. If \code{NULL}, plots are not saved.
#'
#' @return A list with two sublists:
#' \itemize{
#'   \item \code{tables}: \code{infection_summary}, \code{mean_viral_load},
#'         \code{infection_level_counts}, \code{infection_level_percent} (latter two may be \code{NULL})
#'   \item \code{plots}: \code{infection_proportion}, \code{mean_viral_load},
#'         \code{infection_level_counts_plot}, \code{infection_level_percent_plot} (latter two may be \code{NULL})
#' }
#'
#' @examples
#' \dontrun{
#' res <- CellVariableInfection(
#'   obj          = data,
#'   group_col    = "manual_label",
#'   score_col    = "viral_load",
#'   level_col    = "infection_level",
#'   zero_tol     = 0,
#'   order_groups = c("CD4+Tcells","CD8+Tcells","B_cells"),
#'   palette      = c("CD4+Tcells"="#8FBC8F","CD8+Tcells"="#87CEFA","B_cells"="#FFB347"),
#'   lvl_colors   = c("0"="#cfd8dc","1"="#90caf9","2"="#64b5f6","3"="#42a5f5","4"="#1e88e5"),
#'   save_dir     = "results/summary_by_group"
#' )
#' res$plots$infection_proportion
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_text scale_y_continuous
#'   theme_minimal ggtitle theme element_text xlab ylab guides guide_legend
#'   coord_flip ggsave
#' @importFrom dplyr group_by summarise rename mutate ungroup %>%
#' @importFrom scales percent percent_format
#' @export
CellVariableInfection <- function(
    obj,
    group_col   = "manual_label",   # "manual_label" or "timepoint"
    score_col   = "viral_load",     # viral load column
    level_col   = "infection_level",# infection level column (0/1/2/...)
    zero_tol    = 0,                # > zero_tol is considered infected
    order_groups = NULL,            # explicit order for groups
    palette     = NULL,             # optional group colors
    lvl_colors  = NULL,             # optional infection level colors
    save_dir    = NULL              # output directory for saving plots
){
  if (!inherits(obj, "Seurat")) stop("'obj' must be a Seurat object.")
  meta <- obj@meta.data

  # basic column checks
  for (cn in c(group_col, score_col)) {
    if (!cn %in% colnames(meta)) {
      stop(sprintf("Column '%s' not found in meta.data.", cn))
    }
  }
  has_level <- level_col %in% colnames(meta)

  # Ensure clean column types
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

  # ------- Determine group order -------
  if (!is.null(order_groups)) {
    lvl_groups <- order_groups
  } else if (is.factor(obj@meta.data[[group_col]])) {
    lvl_groups <- levels(obj@meta.data[[group_col]])
  } else {
    lvl_groups <- unique(meta[[group_col]])
  }
  meta[[group_col]] <- factor(meta[[group_col]], levels = lvl_groups)

  # =============== ① Infection ratio per group (viral_load > zero_tol) ===============
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

  p_prop <- ggplot2::ggplot(infection_summary, ggplot2::aes(x = group, y = infected_ratio, fill = group)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", infected_ratio * 100)),
                       vjust = -0.3, size = 3.2) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
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

  # =============== ② Mean viral load per group ===============
  mean_viral <- meta %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(mean_viral_load = mean(.data[[score_col]], na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::rename(group = !!group_col) %>%
    dplyr::mutate(group = factor(group, levels = lvl_groups))

  p_mean <- ggplot2::ggplot(mean_viral, ggplot2::aes(x = group, y = mean_viral_load, fill = group)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", mean_viral_load)),
                       vjust = -0.3, size = 3.2) +
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

  # =============== ③ Infection level stacked barplots: count + percentage ===============
  if (has_level) {
    stacked_counts <- meta %>%
      dplyr::group_by(.data[[group_col]], .data[[level_col]]) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::rename(group = !!group_col, level = !!level_col)

    stacked_counts$group <- factor(stacked_counts$group, levels = lvl_groups)

    if (all(grepl("^\\d+$", as.character(stacked_counts$level)))) {
      stacked_counts$level <- factor(as.integer(as.character(stacked_counts$level)),
                                     levels = sort(unique(as.integer(as.character(stacked_counts$level)))))
    } else {
      stacked_counts$level <- factor(as.character(stacked_counts$level),
                                     levels = sort(unique(as.character(stacked_counts$level))))
    }

    # —— Count stacked barplot
    p_stack_count <- ggplot2::ggplot(stacked_counts, ggplot2::aes(
      x = group, y = n, fill = level)) +
      ggplot2::geom_col(width = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = n),
                         position = ggplot2::position_stack(vjust = 0.5), size = 3.0, color = "white") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::ggtitle("Counts of infection levels within each group") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 14),
        plot.title  = ggplot2::element_text(hjust = 0.5, size = 20)
      ) +
      ggplot2::xlab("Group") + ggplot2::ylab("Cell counts") +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Infection level"))
    if (!is.null(lvl_colors)) p_stack_count <- p_stack_count + ggplot2::scale_fill_manual(values = lvl_colors)

    # —— Percentage stacked barplot
    stacked_percent <- stacked_counts %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(frac = n / sum(n)) %>%
      dplyr::ungroup()

    p_stack_percent <- ggplot2::ggplot(stacked_percent, ggplot2::aes(
      x = group, y = frac, fill = level)) +
      ggplot2::geom_col(width = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = scales::percent(frac, accuracy = 1)),
                         position = ggplot2::position_stack(vjust = 0.5), size = 3.0, color = "white") +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
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
    p_stack_count    <- NULL
    p_stack_percent  <- NULL
    stacked_percent  <- NULL
  }

  # =============== Save outputs ===============
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
      infection_level_percent_plot = p_stack_percent
    )
  )
}
