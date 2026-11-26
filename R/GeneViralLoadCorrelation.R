#' Gene–viral load correlation and infection-level trend plot (selected cell groups)
#'
#' @description
#' For one or multiple genes, compute correlation between gene expression and viral load
#' within selected cells, and visualize gene expression dynamics along infection levels.
#' The plot annotates correlation coefficient (r) and p-value.
#'
#' If show_cell_scatter = TRUE, overlays a standard violin distribution at each infection level
#' (no extra packages needed).
#'
#' @param obj A Seurat object.
#' @param genes Character vector. One or multiple gene names (must match features).
#' @param viral_col Character. Meta-data column for viral load (default "viral_load").
#' @param level_col Character. Meta-data column for infection level (default "infection_level").
#' @param group_col Optional character. A meta-data column used to subset groups (e.g. "manual_label").
#' @param group_include Optional character vector. Groups to keep if group_col is provided.
#' @param assay Character. Assay name (default "RNA").
#' @param slot Character. Slot for expression (default "data").
#' @param method Character. Correlation method: "spearman" or "pearson".
#' @param use_infected_only Logical. If TRUE, correlation is computed using cells with infection_level > 0.
#' @param min_cells Integer. Minimum cells required after filtering.
#' @param summarize_fun Function. Summarization per infection level (default mean).
#' @param show_se Logical. Add SE ribbon (default TRUE).
#' @param point_size Numeric. Point size for level means.
#' @param line_width Numeric. Line width for mean trend.
#' @param violin_alpha Numeric. Alpha for violin fill (default 0.25).
#' @param show_cell_scatter Logical. If TRUE, overlay violin distribution (default FALSE).
#' @param violin_width Numeric. Width of violins (default 0.85).
#' @param show_cell_points Logical. If TRUE, overlay per-cell jitter points (default FALSE).
#' @param cell_point_alpha Numeric. Alpha for per-cell points (default 0.12).
#' @param cell_point_size Numeric. Size for per-cell points (default 0.35).
#' @param jitter_width Numeric. Jitter width for per-cell points (default 0.12).
#' @param save_dir Optional directory to save outputs (NULL = not save).
#' @param save_prefix Prefix for saved files.
#' @param cor_label_size Numeric. Size for correlation label text (default 4.6).
#' @param cor_label_vjust Numeric. Vertical adjustment for correlation label (default 0).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{cor_table}: data.frame of per-gene correlation results
#'   \item \code{plot}: ggplot object
#'   \item \code{data_long}: long table used for plotting
#' }
#'
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom stats cor cor.test
#' @export
GeneViralLoadCorrelation <- function(
    obj,
    genes,
    viral_col = "viral_load",
    level_col = "infection_level",
    group_col = NULL,
    group_include = NULL,
    
    assay = "RNA",
    
    # ==========================================================
    # Seurat v5 使用 layer
    # 例如：
    # layer = "data"
    # layer = "counts"
    # ==========================================================
    layer = "data",
    
    method = c("spearman", "pearson"),
    
    use_infected_only = FALSE,
    
    min_cells = 50,
    
    summarize_fun = mean,
    
    # ==========================================================
    # Mean ± SE
    # ==========================================================
    show_se = TRUE,
    
    # mean point
    point_size = 2.2,
    
    # mean line
    line_width = 1.2,
    
    # ==========================================================
    # Violin
    # ==========================================================
    violin_alpha = 0.25,
    
    show_cell_scatter = TRUE,
    
    violin_width = 0.85,
    
    # ==========================================================
    # 单细胞散点
    # ==========================================================
    show_cell_points = TRUE,
    
    cell_point_alpha = 0.15,
    
    cell_point_size = 0.1,
    
    jitter_width = 0.12,
    
    # ==========================================================
    # 保存
    # ==========================================================
    save_dir = NULL,
    
    save_prefix = "Gene_ViralLoad_Correlation",
    
    # ==========================================================
    # correlation label
    # ==========================================================
    cor_label_size = 4.6,
    
    cor_label_vjust = 0,
    
    # ==========================================================
    # 基因名字体大小
    #
    # 控制 facet strip 中：
    # NCL
    # SLC7A5
    # ...
    # ==========================================================
    gene_title_size = 18,
    
    # ==========================================================
    # 整体字体
    # ==========================================================
    base_size = 12,
    
    # ==========================================================
    # 最上方 group title 字体
    # ==========================================================
    plot_title_size = 10,
    
    # ==========================================================
    # 坐标轴标题字体
    # ==========================================================
    axis_title_size = 12,
    
    # ==========================================================
    # 坐标轴刻度字体
    # ==========================================================
    axis_text_size = 11
    
) {
  
  # ============================================================
  # 1. 基本检查
  # ============================================================
  
  stopifnot(
    inherits(
      obj,
      "Seurat"
    )
  )
  
  
  method <- match.arg(
    method
  )
  
  
  # ============================================================
  # 2. dependencies
  # ============================================================
  
  pkgs_need <- c(
    "Seurat",
    "Matrix",
    "dplyr",
    "tidyr",
    "ggplot2",
    "readr"
  )
  
  
  miss_pkg <- pkgs_need[
    !vapply(
      pkgs_need,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]
  
  
  if (length(miss_pkg) > 0) {
    
    stop(
      "Missing packages: ",
      paste(
        miss_pkg,
        collapse = ", "
      )
    )
  }
  
  
  # ============================================================
  # 3. metadata
  # ============================================================
  
  md <- obj@meta.data
  
  
  if (!viral_col %in% colnames(md)) {
    
    stop(
      "Column '",
      viral_col,
      "' not found in meta.data."
    )
  }
  
  
  if (!level_col %in% colnames(md)) {
    
    stop(
      "Column '",
      level_col,
      "' not found in meta.data."
    )
  }
  
  
  if (
    missing(genes) ||
    length(genes) == 0
  ) {
    
    stop(
      "'genes' cannot be empty."
    )
  }
  
  
  # ============================================================
  # 4. 获取表达矩阵
  #
  # Seurat v5:
  # 使用 layer，而不是 slot
  # ============================================================
  
  Seurat::DefaultAssay(
    obj
  ) <- assay
  
  
  expr <- Seurat::GetAssayData(
    object = obj,
    assay = assay,
    layer = layer
  )
  
  
  # ============================================================
  # 5. 检查基因
  # ============================================================
  
  genes_use <- intersect(
    genes,
    rownames(expr)
  )
  
  
  if (length(genes_use) == 0) {
    
    stop(
      "None of the input genes matched assay features."
    )
  }
  
  
  genes_missing <- setdiff(
    genes,
    genes_use
  )
  
  
  if (length(genes_missing) > 0) {
    
    warning(
      "These genes were not found and will be skipped: ",
      paste(
        genes_missing,
        collapse = ", "
      )
    )
  }
  
  
  # ============================================================
  # 6. cell filtering
  # ============================================================
  
  cells <- colnames(
    obj
  )
  
  
  # ============================================================
  # 7. figure title
  # ============================================================
  
  group_title <- "All selected cells"
  
  
  if (
    !is.null(group_include) &&
    length(group_include) > 0
  ) {
    
    group_title <- paste(
      group_include,
      collapse = ", "
    )
    
  } else if (!is.null(group_col)) {
    
    group_title <- paste0(
      "Cells grouped by ",
      group_col
    )
  }
  
  
  # ============================================================
  # 8. 根据 group 过滤
  # ============================================================
  
  if (!is.null(group_col)) {
    
    if (!group_col %in% colnames(md)) {
      
      stop(
        "Column '",
        group_col,
        "' not found in meta.data."
      )
    }
    
    
    if (!is.null(group_include)) {
      
      keep_g <- as.character(
        md[[group_col]]
      ) %in%
        as.character(
          group_include
        )
      
      
      cells <- rownames(md)[
        keep_g
      ]
    }
  }
  
  
  # ============================================================
  # 9. viral load & infection level
  # ============================================================
  
  vload <- suppressWarnings(
    as.numeric(
      as.character(
        md[
          cells,
          viral_col
        ]
      )
    )
  )
  
  
  lvl <- suppressWarnings(
    as.numeric(
      as.character(
        md[
          cells,
          level_col
        ]
      )
    )
  )
  
  
  # ============================================================
  # 10. complete cases
  # ============================================================
  
  keep <- is.finite(
    vload
  ) &
    is.finite(
      lvl
    )
  
  
  cells <- cells[
    keep
  ]
  
  
  vload <- vload[
    keep
  ]
  
  
  lvl <- lvl[
    keep
  ]
  
  
  # ============================================================
  # 11. 是否只保留感染细胞
  # ============================================================
  
  if (use_infected_only) {
    
    keep2 <- lvl > 0
    
    
    cells <- cells[
      keep2
    ]
    
    
    vload <- vload[
      keep2
    ]
    
    
    lvl <- lvl[
      keep2
    ]
  }
  
  
  # ============================================================
  # 12. 最低细胞数
  # ============================================================
  
  if (length(cells) < min_cells) {
    
    stop(
      "Too few cells after filtering: ",
      length(cells),
      " (< min_cells)."
    )
  }
  
  
  # ============================================================
  # 13. subset expression
  # ============================================================
  
  expr_sub <- expr[
    genes_use,
    cells,
    drop = FALSE
  ]
  
  
  # ============================================================
  # 14. long table
  # ============================================================
  
  df_long <- as.data.frame(
    as.matrix(
      t(
        expr_sub
      )
    )
  )
  
  
  df_long$cell <- rownames(
    df_long
  )
  
  
  df_long[[viral_col]] <-
    vload
  
  
  df_long[[level_col]] <-
    lvl
  
  
  df_long <- tidyr::pivot_longer(
    df_long,
    cols = dplyr::all_of(
      genes_use
    ),
    names_to = "gene",
    values_to = "expr"
  )
  
  
  # ============================================================
  # 15. infection level numeric
  # ============================================================
  
  df_long[[level_col]] <-
    suppressWarnings(
      as.numeric(
        as.character(
          df_long[[level_col]]
        )
      )
    )
  
  
  if (
    any(
      !is.finite(
        df_long[[level_col]]
      )
    )
  ) {
    
    stop(
      "Non-numeric infection_level detected after coercion."
    )
  }
  
  
  # ============================================================
  # 16. correlation per gene
  # ============================================================
  
  cor_table <- df_long |>
    
    dplyr::group_by(
      gene
    ) |>
    
    dplyr::summarise(
      
      n = dplyr::n(),
      
      r = suppressWarnings(
        stats::cor(
          expr,
          .data[[viral_col]],
          method = method,
          use = "pairwise.complete.obs"
        )
      ),
      
      p = suppressWarnings(
        stats::cor.test(
          expr,
          .data[[viral_col]],
          method = method
        )$p.value
      ),
      
      .groups = "drop"
    ) |>
    
    dplyr::mutate(
      
      label = sprintf(
        "r(%s)=%.2f, p=%.2g",
        method,
        r,
        p
      )
    )
  
  
  # ============================================================
  # 17. summarize across infection level
  # ============================================================
  
  df_sum <- df_long |>
    
    dplyr::group_by(
      gene,
      .data[[level_col]]
    ) |>
    
    dplyr::summarise(
      
      mean_expr = summarize_fun(
        expr,
        na.rm = TRUE
      ),
      
      se_expr = stats::sd(
        expr,
        na.rm = TRUE
      ) /
        sqrt(
          dplyr::n()
        ),
      
      n_cells = dplyr::n(),
      
      .groups = "drop"
    ) |>
    
    dplyr::rename(
      level = dplyr::all_of(
        level_col
      )
    ) |>
    
    dplyr::mutate(
      level = as.numeric(
        level
      )
    )
  
  
  # ============================================================
  # 18. correlation label position
  # ============================================================
  
  df_pos <- df_long |>
    
    dplyr::group_by(
      gene
    ) |>
    
    dplyr::summarise(
      
      x_min = min(
        .data[[level_col]],
        na.rm = TRUE
      ),
      
      x_max = max(
        .data[[level_col]],
        na.rm = TRUE
      ),
      
      y_min = min(
        expr,
        na.rm = TRUE
      ),
      
      y_max = max(
        expr,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) |>
    
    dplyr::left_join(
      
      cor_table |>
        dplyr::select(
          gene,
          label
        ),
      
      by = "gene"
    ) |>
    
    dplyr::mutate(
      
      x = x_min +
        0.05 *
        (
          x_max -
            x_min +
            1e-9
        ),
      
      y = y_max +
        0.06 *
        (
          y_max -
            y_min +
            1e-9
        )
    )
  
  
  # ============================================================
  # 19. X breaks
  # ============================================================
  
  x_breaks <- sort(
    unique(
      df_sum$level
    )
  )
  
  
  # ============================================================
  # 20. violin data
  #
  # 只使用 expression > 0
  # ============================================================
  
  df_violin <- df_long |>
    
    dplyr::filter(
      is.finite(
        expr
      ) &
        expr > 0
    )
  
  
  # ============================================================
  # 21. plot
  # ============================================================
  
  p <- ggplot2::ggplot()
  
  
  # ============================================================
  # 22. Violin
  # ============================================================
  
  if (isTRUE(show_cell_scatter)) {
    
    p <- p +
      
      ggplot2::geom_violin(
        
        data = df_violin,
        
        ggplot2::aes(
          x = .data[[level_col]],
          y = expr,
          group = .data[[level_col]]
        ),
        
        inherit.aes = FALSE,
        
        fill = "#D9902F",
        
        alpha = violin_alpha,
        
        color = NA,
        
        width = violin_width,
        
        trim = TRUE,
        
        bw = "nrd0",
        
        scale = "width"
      )
  }
  
  
  # ============================================================
  # 23. 单细胞点
  # ============================================================
  
  if (isTRUE(show_cell_points)) {
    
    p <- p +
      
      ggplot2::geom_point(
        
        data = df_long,
        
        ggplot2::aes(
          x = .data[[level_col]],
          y = expr
        ),
        
        inherit.aes = FALSE,
        
        position =
          ggplot2::position_jitter(
            width = jitter_width,
            height = 0
          ),
        
        size = cell_point_size,
        
        alpha = cell_point_alpha,
        
        color = "#D9902F"
      )
  }
  
  
  # ============================================================
  # 24. SE ribbon
  # ============================================================
  
  if (isTRUE(show_se)) {
    
    p <- p +
      
      ggplot2::geom_ribbon(
        
        data = df_sum,
        
        ggplot2::aes(
          x = level,
          ymin = mean_expr - se_expr,
          ymax = mean_expr + se_expr
        ),
        
        inherit.aes = FALSE,
        
        fill = "#9E2F2F",
        
        alpha = 0.18
      )
  }
  
  
  # ============================================================
  # 25. Mean line + points
  # ============================================================
  
  p <- p +
    
    ggplot2::geom_line(
      
      data = df_sum,
      
      ggplot2::aes(
        x = level,
        y = mean_expr
      ),
      
      inherit.aes = FALSE,
      
      linewidth = line_width,
      
      color = "#9E2F2F"
    ) +
    
    
    ggplot2::geom_point(
      
      data = df_sum,
      
      ggplot2::aes(
        x = level,
        y = mean_expr
      ),
      
      inherit.aes = FALSE,
      
      size = point_size,
      
      color = "#9E2F2F"
    ) +
    
    
    # ==========================================================
  # gene facet
  # ==========================================================
  
  ggplot2::facet_wrap(
    ~gene,
    scales = "free_y"
  ) +
    
    
    # ==========================================================
  # theme
  # ==========================================================
  
  ggplot2::theme_classic(
    base_size = base_size
  ) +
    
    
    # ==========================================================
  # labels
  # ==========================================================
  
  ggplot2::labs(
    
    title = group_title,
    
    x = "Infection level",
    
    y = paste0(
      "Expression (",
      layer,
      ")"
    )
  ) +
    
    
    # ==========================================================
  # 字体设置
  # ==========================================================
  
  ggplot2::theme(
    
    # --------------------------------------------------------
    # 最上方 group title
    # --------------------------------------------------------
    
    plot.title =
      ggplot2::element_text(
        
        hjust = 0.5,
        
        face = "bold",
        
        size = plot_title_size
      ),
    
    
    # --------------------------------------------------------
    # 坐标轴标题
    # --------------------------------------------------------
    
    axis.title =
      ggplot2::element_text(
        
        face = "bold",
        
        size = axis_title_size
      ),
    
    
    # --------------------------------------------------------
    # 坐标轴刻度
    # --------------------------------------------------------
    
    axis.text =
      ggplot2::element_text(
        
        size = axis_text_size
      ),
    
    
    # --------------------------------------------------------
    # 关键：
    # 基因名字体
    #
    # 例如：
    # NCL
    # SLC7A5
    # IGF1R
    # --------------------------------------------------------
    
    strip.text =
      ggplot2::element_text(
        
        face = "bold",
        
        size = gene_title_size
      )
  ) +
    
    
    # ==========================================================
  # Correlation label
  # ==========================================================
  
  ggplot2::geom_text(
    
    data = df_pos,
    
    ggplot2::aes(
      x = x,
      y = y,
      label = label
    ),
    
    inherit.aes = FALSE,
    
    hjust = 0,
    
    vjust = cor_label_vjust,
    
    size = cor_label_size,
    
    fontface = "bold"
  ) +
    
    
    # ==========================================================
  # X axis
  # ==========================================================
  
  ggplot2::scale_x_continuous(
    
    breaks = x_breaks,
    
    limits = range(
      x_breaks
    )
  ) +
    
    
    ggplot2::coord_cartesian(
      clip = "off"
    )
  
  
  # ============================================================
  # 26. 保存
  # ============================================================
  
  if (!is.null(save_dir)) {
    
    dir.create(
      save_dir,
      showWarnings = FALSE,
      recursive = TRUE
    )
    
    
    readr::write_csv(
      
      cor_table,
      
      file.path(
        save_dir,
        paste0(
          save_prefix,
          "_cor_table.csv"
        )
      )
    )
    
    
    ggplot2::ggsave(
      
      filename = file.path(
        save_dir,
        paste0(
          save_prefix,
          "_trend.pdf"
        )
      ),
      
      plot = p,
      
      width = 11,
      
      height = 6.5,
      
      useDingbats = FALSE
    )
  }
  
  
  # ============================================================
  # 27. return
  # ============================================================
  
  invisible(
    list(
      
      cor_table = cor_table,
      
      plot = p,
      
      data_long = df_long,
      
      summary_data = df_sum,
      
      genes_used = genes_use,
      
      gene_title_size = gene_title_size
    )
  )
}
