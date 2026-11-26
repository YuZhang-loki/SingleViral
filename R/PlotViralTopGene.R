#' Visualize dominant viral-gene cells across groups
#'
#' @description
#' Visualize cell-level viral load across user-defined groups and highlight
#' cells in which a specified viral gene is the dominant viral transcript.
#'
#' Viral-gene dominance is determined using raw UMI counts from the
#' \code{counts} layer of the specified Seurat assay. For each cell, the
#' expression of the target viral gene is compared with all viral genes
#' supplied through \code{virus_gene}. Cells satisfying the dominance
#' criterion can be highlighted using a distinct point shape.
#'
#' The function additionally supports:
#' \itemize{
#'   \item filtering selected infection levels;
#'   \item requiring a minimum target-gene UMI count;
#'   \item requiring a minimum total viral UMI count;
#'   \item allowing or excluding ties for the highest expressed viral gene;
#'   \item displaying a violin-density background;
#'   \item adding borders to ordinary or dominant-gene cells;
#'   \item independently controlling point size, transparency, colors, and legends;
#'   \item summarizing the number and fraction of dominant-gene cells in each group.
#' }
#'
#' @param obj A Seurat object.
#'
#' @param group_col Character. Metadata column used to define groups on the
#' x-axis.
#'
#' @param score_col Character. Metadata column containing the quantitative
#' score displayed on the y-axis, such as viral load.
#'
#' @param level_col Character. Metadata column defining infection levels used
#' for cell coloring and optional filtering.
#'
#' @param lvl_colors Named character vector specifying colors for infection
#' levels. Names must correspond to values in \code{level_col}.
#' For example:
#' \code{c("1" = "#4DBBD5", "2" = "#44C37F",
#' "3" = "#F3A13B", "4" = "#E64B35")}.
#'
#' @param gene Character. Target viral gene to identify as the dominant viral
#' gene in each cell.
#'
#' @param virus_gene Character vector containing all viral genes used for
#' determining viral-gene dominance.
#'
#' @param assay Character. Seurat assay containing the viral-gene expression
#' matrix. Default is \code{"RNA"}.
#'
#' @param infection_levels Optional character or numeric vector specifying
#' infection levels to retain. If \code{NULL}, all available infection levels
#' are included.
#'
#' @param min_gene_count Numeric. Minimum raw UMI count required for the target
#' gene before a cell can be classified as a dominant-gene cell.
#' Default is 2.
#'
#' @param min_total_viral_count Numeric. Minimum total viral UMI count across
#' all genes in \code{virus_gene} required before a cell can be classified as
#' a dominant-gene cell. Default is 0.
#'
#' @param gene_color Character. Color used for cells in which \code{gene} is
#' the dominant viral gene. Default is \code{"#B2182B"}.
#'
#' @param gene_shape Numeric. Point shape used for dominant-gene cells when
#' \code{gene_border = FALSE}. Default is 17.
#'
#' @param gene_size Numeric. Point size of dominant-gene cells in the main
#' plot. Default is 2.5.
#'
#' @param gene_alpha Numeric between 0 and 1. Transparency of dominant-gene
#' cells. Default is 0.7.
#'
#' @param gene_legend_label Optional character. Legend label for dominant-gene
#' cells. If \code{NULL}, the label is automatically generated as
#' \code{"<gene> dominant"}.
#'
#' @param gene_border Logical. If \code{TRUE}, dominant-gene cells are drawn
#' with a separate border using point shape 24. Default is \code{TRUE}.
#'
#' @param gene_border_color Character. Border color for dominant-gene cells
#' when \code{gene_border = TRUE}. Default is \code{"grey30"}.
#'
#' @param gene_border_width Numeric. Border width for dominant-gene cells.
#' Default is 0.2.
#'
#' @param point_size Numeric. Point size for ordinary cells in the main plot.
#' Default is 1.8.
#'
#' @param point_alpha Numeric between 0 and 1. Transparency of ordinary cells.
#' Default is 0.75.
#'
#' @param point_border Logical. If \code{TRUE}, ordinary cells are plotted
#' using shape 21 with infection-level color as the fill and
#' \code{point_border_color} as the border. Default is \code{TRUE}.
#'
#' @param point_border_color Character. Border color for ordinary cells when
#' \code{point_border = TRUE}. Default is \code{"grey30"}.
#'
#' @param point_border_width Numeric. Border width for ordinary cells.
#' Default is 0.1.
#'
#' @param show_violin Logical. If \code{TRUE}, overlays a violin-density
#' background for each group. Default is \code{TRUE}.
#'
#' @param violin_fill Character. Fill color of the violin background.
#' Default is \code{"grey70"}.
#'
#' @param violin_alpha Numeric between 0 and 1. Transparency of the violin
#' background. Default is 0.25.
#'
#' @param violin_width Numeric. Width of violin distributions.
#' Default is 0.75.
#'
#' @param violin_color Character or \code{NA}. Border color of violin
#' distributions. Use \code{NA} to remove violin borders.
#'
#' @param violin_linewidth Numeric. Width of violin borders.
#' Default is 0.4.
#'
#' @param violin_trim Logical. Passed to \code{ggplot2::geom_violin()}.
#' If \code{TRUE}, violin distributions are trimmed to the observed data range.
#'
#' @param jitter_width Numeric. Horizontal jitter width applied to cell points.
#' Default is 0.20.
#'
#' @param jitter_seed Integer. Random seed used for reproducible jittering.
#' Default is 123.
#'
#' @param include_ties Logical. Whether cells in which the target gene is tied
#' with another viral gene for the highest UMI count should be classified as
#' dominant-gene cells. If \code{FALSE}, the target gene must be strictly
#' higher than all other viral genes.
#'
#' @param y_n_intervals Integer. Number of equal intervals used to divide the
#' y-axis between 0 and the selected maximum. Default is 5.
#'
#' @param y_max Optional numeric. Maximum value of the y-axis. If \code{NULL},
#' the maximum value observed in \code{score_col} is used.
#'
#' @param group_levels Optional character vector defining the plotting order
#' of groups on the x-axis. If \code{NULL}, the order is inferred from the
#' metadata.
#'
#' @param show_n Logical. If \code{TRUE}, displays the total number of cells and
#' the number of dominant-gene cells beneath each x-axis group.
#'
#' @param title Optional character. Plot title. If \code{NULL}, an automatic
#' title of the form \code{"<gene> dominant infected cells"} is generated.
#'
#' @param x_title Optional character. X-axis title. If \code{NULL}, no x-axis
#' title is shown.
#'
#' @param y_title Optional character. Y-axis title. If \code{NULL},
#' \code{score_col} is used.
#'
#' @param base_size Numeric. Base font size used by
#' \code{ggplot2::theme_classic()}. Default is 14.
#'
#' @param legend_title_size Numeric. Font size of legend titles.
#' Default is 14.
#'
#' @param legend_text_size Numeric. Font size of legend labels.
#' Default is 14.
#'
#' @param legend_point_size Numeric. Size of ordinary-cell symbols in the
#' infection-level legend. This is independent of \code{point_size}.
#' Default is 5.
#'
#' @param legend_gene_size Numeric. Size of the dominant-gene symbol in the
#' legend. This is independent of \code{gene_size}. Default is 6.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{plot}: ggplot object.
#'   \item \code{data}: data frame used for plotting.
#'   \item \code{gene_top_cells}: subset of cells classified as dominant for
#'   the target viral gene.
#'   \item \code{group_stats}: per-group cell counts and dominant-gene fractions.
#'   \item \code{y_breaks}: y-axis break positions.
#'   \item \code{y_max}: final y-axis maximum.
#'   \item \code{infection_levels}: infection levels used in the analysis.
#'   \item \code{min_gene_count}: target-gene UMI threshold.
#'   \item \code{min_total_viral_count}: total viral UMI threshold.
#'   \item \code{colors_use}: infection-level colors used in the plot.
#'   \item \code{level_order}: final infection-level plotting order.
#'   \item \code{virus_gene_present}: viral genes found in the counts layer.
#'   \item \code{virus_gene_missing}: requested viral genes not found.
#' }
#'
#' @details
#' Dominant viral-gene classification is based exclusively on raw UMI counts
#' from the \code{counts} layer. For each cell, the UMI count of the target
#' gene is compared with the maximum UMI count among all other genes listed in
#' \code{virus_gene}.
#'
#' When \code{include_ties = FALSE}, a cell is classified as dominant only if:
#'
#' \deqn{
#' count_{target} > max(count_{other\ viral\ genes})
#' }
#'
#' and both:
#'
#' \deqn{
#' count_{target} >= min\_gene\_count
#' }
#'
#' and
#'
#' \deqn{
#' total\ viral\ UMI >= min\_total\_viral\_count
#' }
#'
#' are satisfied.
#'
#' @examples
#' \dontrun{
#'
#' custom_level_cols <- c(
#'   "1" = "#4DBBD5",
#'   "2" = "#44C37F",
#'   "3" = "#F3A13B",
#'   "4" = "#E64B35"
#' )
#'
#' res <- PlotViralTopGene(
#'   obj = data,
#'   group_col = "infection_level",
#'   score_col = "viral_load",
#'   level_col = "infection_level",
#'   lvl_colors = custom_level_cols,
#'   gene = "U83R",
#'   virus_gene = virus_gene,
#'   infection_levels = c("1", "2", "3", "4"),
#'   min_gene_count = 5,
#'   min_total_viral_count = 50,
#'   include_ties = FALSE,
#'   point_border = TRUE,
#'   gene_border = TRUE,
#'   show_violin = TRUE,
#'   y_title = "Viral load"
#' )
#'
#' res$plot
#' res$group_stats
#' res$gene_top_cells
#' }
#'
#' @export
PlotViralTopGene <- function(
    obj,
    group_col,
    score_col,
    level_col,
    lvl_colors,
    gene,
    virus_gene,
    assay = "RNA",
    
    # ==========================================================
    # 选择感染等级
    # NULL = 所有感染等级
    # ==========================================================
    infection_levels = NULL,
    
    # ==========================================================
    # 目标基因最低原始 UMI
    # ==========================================================
    min_gene_count = 2,
    
    # ==========================================================
    # 单细胞总体病毒 UMI 最低阈值
    # ==========================================================
    min_total_viral_count = 0,
    
    # ==========================================================
    # gene-highest 特殊细胞
    # ==========================================================
    gene_color = "#B2182B",
    gene_shape = 17,
    gene_size = 2.5,
    gene_alpha = 0.7,
    
    # 特殊细胞图例名称
    gene_legend_label = NULL,
    
    # ==========================================================
    # 特殊三角形是否加边界
    # TRUE 时自动使用 shape = 24
    # ==========================================================
    gene_border = TRUE,
    gene_border_color = "grey30",
    gene_border_width = 0.2,
    
    # ==========================================================
    # 普通细胞
    # ==========================================================
    point_size = 1.8,
    point_alpha = 0.75,
    
    # ==========================================================
    # 普通细胞是否加边界
    #
    # TRUE:
    # shape = 21
    # fill = infection level
    # colour = point_border_color
    #
    # FALSE:
    # 仍使用 shape = 21
    # 但没有边框
    # ==========================================================
    point_border = TRUE,
    point_border_color = "grey30",
    point_border_width = 0.1,
    
    # ==========================================================
    # 灰色 violin 虚影
    # ==========================================================
    show_violin = TRUE,
    violin_fill = "grey70",
    violin_alpha = 0.25,
    violin_width = 0.75,
    violin_color = NA,
    violin_linewidth = 0.4,
    violin_trim = TRUE,
    
    # ==========================================================
    # jitter
    # ==========================================================
    jitter_width = 0.20,
    jitter_seed = 123,
    
    # ==========================================================
    # 是否允许并列最高
    # ==========================================================
    include_ties = FALSE,
    
    # ==========================================================
    # Y轴
    # ==========================================================
    y_n_intervals = 5,
    y_max = NULL,
    
    # ==========================================================
    # 分组顺序
    # ==========================================================
    group_levels = NULL,
    
    # 是否显示每组数量
    show_n = TRUE,
    
    # ==========================================================
    # 标题
    # ==========================================================
    title = NULL,
    x_title = NULL,
    y_title = NULL,
    
    # ==========================================================
    # 主图字体
    # ==========================================================
    base_size = 14,
    
    # ==========================================================
    # 图例文字
    # ==========================================================
    legend_title_size = 14,
    legend_text_size = 14,
    
    # ==========================================================
    # 图例中的细胞点大小
    # 与主图大小独立
    # ==========================================================
    legend_point_size = 5,
    legend_gene_size = 6
) {
  
  # ============================================================
  # 1. 基本检查
  # ============================================================
  
  if (!inherits(obj, "Seurat")) {
    stop("'obj' must be a Seurat object.")
  }
  
  meta <- obj@meta.data
  
  required_cols <- c(
    group_col,
    score_col,
    level_col
  )
  
  missing_cols <- setdiff(
    required_cols,
    colnames(meta)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      "The following metadata columns were not found: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  
  # ============================================================
  # 2. 参数检查
  # ============================================================
  
  if (
    !is.numeric(min_gene_count) ||
    length(min_gene_count) != 1 ||
    is.na(min_gene_count) ||
    min_gene_count < 0
  ) {
    stop(
      "'min_gene_count' must be a single non-negative numeric value."
    )
  }
  
  if (
    !is.numeric(min_total_viral_count) ||
    length(min_total_viral_count) != 1 ||
    is.na(min_total_viral_count) ||
    min_total_viral_count < 0
  ) {
    stop(
      "'min_total_viral_count' must be a single non-negative numeric value."
    )
  }
  
  if (
    !is.numeric(gene_alpha) ||
    length(gene_alpha) != 1 ||
    is.na(gene_alpha) ||
    gene_alpha < 0 ||
    gene_alpha > 1
  ) {
    stop(
      "'gene_alpha' must be between 0 and 1."
    )
  }
  
  if (
    !is.numeric(point_alpha) ||
    length(point_alpha) != 1 ||
    is.na(point_alpha) ||
    point_alpha < 0 ||
    point_alpha > 1
  ) {
    stop(
      "'point_alpha' must be between 0 and 1."
    )
  }
  
  if (
    !is.numeric(violin_alpha) ||
    length(violin_alpha) != 1 ||
    is.na(violin_alpha) ||
    violin_alpha < 0 ||
    violin_alpha > 1
  ) {
    stop(
      "'violin_alpha' must be between 0 and 1."
    )
  }
  
  if (
    !is.logical(point_border) ||
    length(point_border) != 1
  ) {
    stop(
      "'point_border' must be TRUE or FALSE."
    )
  }
  
  if (
    !is.logical(gene_border) ||
    length(gene_border) != 1
  ) {
    stop(
      "'gene_border' must be TRUE or FALSE."
    )
  }
  
  
  # ============================================================
  # 3. 特殊细胞图例名称
  # ============================================================
  
  if (is.null(gene_legend_label)) {
    gene_legend_label <- paste0(
      gene,
      " dominant"
    )
  }
  
  
  # ============================================================
  # 4. 根据 infection_levels 筛选
  # ============================================================
  
  if (!is.null(infection_levels)) {
    
    infection_levels <- as.character(
      infection_levels
    )
    
    available_levels <- unique(
      as.character(
        meta[[level_col]]
      )
    )
    
    missing_levels <- setdiff(
      infection_levels,
      available_levels
    )
    
    if (length(missing_levels) > 0) {
      warning(
        "The following infection levels were not found: ",
        paste(
          missing_levels,
          collapse = ", "
        )
      )
    }
    
    keep_cells <- rownames(meta)[
      as.character(
        meta[[level_col]]
      ) %in% infection_levels
    ]
    
    meta <- meta[
      keep_cells,
      ,
      drop = FALSE
    ]
    
    if (nrow(meta) == 0) {
      stop(
        "No cells remain after filtering by infection_levels."
      )
    }
    
    message(
      "Selected infection levels: ",
      paste(
        infection_levels,
        collapse = ", "
      )
    )
    
    message(
      "Cells after infection-level filtering: ",
      nrow(meta)
    )
  }
  
  
  # ============================================================
  # 5. 获取原始 UMI counts
  # ============================================================
  
  counts <- SeuratObject::LayerData(
    object = obj,
    assay = assay,
    layer = "counts"
  )
  
  message(
    "Viral gene dominance is calculated using raw UMI counts."
  )
  
  
  # ============================================================
  # 6. 检查病毒基因
  # ============================================================
  
  virus_gene_present <- intersect(
    virus_gene,
    rownames(counts)
  )
  
  virus_gene_missing <- setdiff(
    virus_gene,
    rownames(counts)
  )
  
  if (length(virus_gene_present) == 0) {
    stop(
      "None of the genes in 'virus_gene' were found in the counts layer."
    )
  }
  
  if (!gene %in% rownames(counts)) {
    stop(
      paste0(
        "Target gene '",
        gene,
        "' was not found in the counts layer."
      )
    )
  }
  
  if (!gene %in% virus_gene_present) {
    stop(
      paste0(
        "Target gene '",
        gene,
        "' is not included in 'virus_gene'."
      )
    )
  }
  
  message(
    "Virus genes provided: ",
    length(virus_gene)
  )
  
  message(
    "Virus genes found in counts layer: ",
    length(virus_gene_present)
  )
  
  if (length(virus_gene_missing) > 0) {
    message(
      "Virus genes not found: ",
      paste(
        virus_gene_missing,
        collapse = ", "
      )
    )
  }
  
  
  # ============================================================
  # 7. 对齐细胞
  # ============================================================
  
  common_cells <- intersect(
    rownames(meta),
    colnames(counts)
  )
  
  if (length(common_cells) == 0) {
    stop(
      "No common cells were found between metadata and counts matrix."
    )
  }
  
  meta_plot <- meta[
    common_cells,
    ,
    drop = FALSE
  ]
  
  viral_counts <- counts[
    virus_gene_present,
    common_cells,
    drop = FALSE
  ]
  
  
  # ============================================================
  # 8. 获取目标基因 UMI
  # ============================================================
  
  gene_count <- as.numeric(
    viral_counts[
      gene,
      ,
      drop = TRUE
    ]
  )
  
  
  # ============================================================
  # 9. 病毒 UMI 信息
  # ============================================================
  
  viral_dense <- as.matrix(
    viral_counts
  )
  
  total_viral_count <- colSums(
    viral_dense,
    na.rm = TRUE
  )
  
  max_viral_count <- apply(
    viral_dense,
    2,
    max,
    na.rm = TRUE
  )
  
  n_detected_viral_genes <- colSums(
    viral_dense > 0,
    na.rm = TRUE
  )
  
  
  # ============================================================
  # 10. gene fraction
  # ============================================================
  
  gene_fraction <- rep(
    0,
    length(gene_count)
  )
  
  positive_total <- total_viral_count > 0
  
  gene_fraction[positive_total] <-
    gene_count[positive_total] /
    total_viral_count[positive_total]
  
  
  # ============================================================
  # 11. 判断目标 gene 是否为最高表达
  # ============================================================
  
  if (include_ties) {
    
    gene_is_top <-
      gene_count == max_viral_count
    
  } else {
    
    other_genes <- setdiff(
      virus_gene_present,
      gene
    )
    
    if (length(other_genes) == 0) {
      
      gene_is_top <- rep(
        TRUE,
        length(gene_count)
      )
      
    } else {
      
      other_max <- apply(
        viral_dense[
          other_genes,
          ,
          drop = FALSE
        ],
        2,
        max,
        na.rm = TRUE
      )
      
      gene_is_top <-
        gene_count > other_max
    }
  }
  
  
  # ============================================================
  # 12. 加入 UMI 阈值
  # ============================================================
  
  gene_is_top <-
    gene_is_top &
    gene_count >= min_gene_count &
    total_viral_count >= min_total_viral_count
  
  
  # ============================================================
  # 13. 构建绘图数据
  # ============================================================
  
  plot_data <- data.frame(
    
    cell = common_cells,
    
    group = meta_plot[[group_col]],
    
    score = as.numeric(
      meta_plot[[score_col]]
    ),
    
    level = meta_plot[[level_col]],
    
    gene_count = gene_count,
    
    max_viral_count = max_viral_count,
    
    total_viral_count = total_viral_count,
    
    gene_fraction = gene_fraction,
    
    n_detected_viral_genes =
      n_detected_viral_genes,
    
    gene_is_top = gene_is_top,
    
    stringsAsFactors = FALSE
  )
  
  
  # ============================================================
  # 14. 去除 NA
  # ============================================================
  
  plot_data <- plot_data[
    !is.na(plot_data$group) &
      !is.na(plot_data$score) &
      !is.na(plot_data$level),
    ,
    drop = FALSE
  ]
  
  if (nrow(plot_data) == 0) {
    stop(
      "No valid cells remain after removing NA values."
    )
  }
  
  
  # ============================================================
  # 15. group 顺序
  # ============================================================
  
  if (!is.null(group_levels)) {
    
    group_levels_use <- group_levels[
      group_levels %in%
        unique(
          as.character(
            plot_data$group
          )
        )
    ]
    
    plot_data$group <- factor(
      as.character(
        plot_data$group
      ),
      levels = group_levels_use
    )
    
  } else if (is.factor(meta_plot[[group_col]])) {
    
    original_levels <- levels(
      meta_plot[[group_col]]
    )
    
    group_levels_use <- original_levels[
      original_levels %in%
        unique(
          as.character(
            plot_data$group
          )
        )
    ]
    
    plot_data$group <- factor(
      as.character(
        plot_data$group
      ),
      levels = group_levels_use
    )
    
  } else {
    
    plot_data$group <- factor(
      as.character(
        plot_data$group
      ),
      levels = unique(
        as.character(
          plot_data$group
        )
      )
    )
  }
  
  
  # ============================================================
  # 16. infection level 顺序
  # ============================================================
  
  if (is.null(names(lvl_colors))) {
    stop(
      "'lvl_colors' must be a named vector."
    )
  }
  
  if (!is.null(infection_levels)) {
    
    level_order <- infection_levels[
      infection_levels %in%
        names(lvl_colors)
    ]
    
  } else {
    
    level_order <- names(
      lvl_colors
    )
  }
  
  
  present_levels <- unique(
    as.character(
      plot_data$level
    )
  )
  
  
  level_order <- level_order[
    level_order %in%
      present_levels
  ]
  
  
  missing_color_levels <- setdiff(
    present_levels,
    names(lvl_colors)
  )
  
  
  if (length(missing_color_levels) > 0) {
    
    stop(
      "No color supplied in 'lvl_colors' for levels: ",
      paste(
        missing_color_levels,
        collapse = ", "
      )
    )
  }
  
  
  plot_data$level <- factor(
    as.character(
      plot_data$level
    ),
    levels = level_order
  )
  
  
  # ============================================================
  # 17. 最终颜色顺序
  #
  # 强制按照 factor level 顺序排列
  # ============================================================
  
  colors_use <- lvl_colors[
    level_order
  ]
  
  
  # ============================================================
  # 18. 普通 / 特殊细胞
  # ============================================================
  
  normal_data <- plot_data[
    !plot_data$gene_is_top,
    ,
    drop = FALSE
  ]
  
  gene_top_data <- plot_data[
    plot_data$gene_is_top,
    ,
    drop = FALSE
  ]
  
  gene_top_data$special_type <-
    gene_legend_label
  
  
  # ============================================================
  # 19. 输出信息
  # ============================================================
  
  message(
    "Total cells used for plotting: ",
    nrow(plot_data)
  )
  
  message(
    "Minimum ",
    gene,
    " UMI required: ",
    min_gene_count
  )
  
  message(
    "Minimum total viral UMI required: ",
    min_total_viral_count
  )
  
  message(
    "Cells with ",
    gene,
    " as dominant viral gene: ",
    nrow(gene_top_data)
  )
  
  
  # ============================================================
  # 20. 每组统计
  # ============================================================
  
  group_names <- levels(
    droplevels(
      plot_data$group
    )
  )
  
  group_stats <- data.frame(
    group = group_names,
    n_total = 0,
    n_gene_top = 0,
    gene_top_fraction = 0,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(group_names)) {
    
    g <- group_names[i]
    
    temp <- plot_data[
      as.character(
        plot_data$group
      ) == g,
      ,
      drop = FALSE
    ]
    
    group_stats$n_total[i] <-
      nrow(temp)
    
    group_stats$n_gene_top[i] <-
      sum(
        temp$gene_is_top,
        na.rm = TRUE
      )
    
    if (group_stats$n_total[i] > 0) {
      
      group_stats$gene_top_fraction[i] <-
        group_stats$n_gene_top[i] /
        group_stats$n_total[i]
    }
  }
  
  
  # ============================================================
  # 21. X轴标签
  # ============================================================
  
  if (show_n) {
    
    x_labels <- setNames(
      
      paste0(
        group_stats$group,
        "\n",
        "n = ",
        group_stats$n_total,
        "\n",
        "\u25B2 n = ",
        group_stats$n_gene_top
      ),
      
      group_stats$group
    )
    
  } else {
    
    x_labels <- setNames(
      group_stats$group,
      group_stats$group
    )
  }
  
  
  # ============================================================
  # 22. Y轴最大值
  # ============================================================
  
  data_max <- max(
    plot_data$score,
    na.rm = TRUE
  )
  
  if (is.null(y_max)) {
    
    final_y_max <- data_max
    
  } else {
    
    if (
      !is.numeric(y_max) ||
      length(y_max) != 1 ||
      is.na(y_max) ||
      y_max <= 0
    ) {
      stop(
        "'y_max' must be a single positive numeric value."
      )
    }
    
    final_y_max <- y_max
    
    if (data_max > y_max) {
      warning(
        "Some score values are larger than y_max and ",
        "will not be displayed."
      )
    }
  }
  
  if (final_y_max <= 0) {
    stop(
      "Maximum value of score_col must be greater than 0."
    )
  }
  
  
  # ============================================================
  # 23. Y轴均分
  # ============================================================
  
  if (
    !is.numeric(y_n_intervals) ||
    length(y_n_intervals) != 1 ||
    y_n_intervals < 1
  ) {
    stop(
      "'y_n_intervals' must be a positive integer."
    )
  }
  
  y_n_intervals <-
    as.integer(
      y_n_intervals
    )
  
  y_breaks <- seq(
    from = 0,
    to = final_y_max,
    length.out = y_n_intervals + 1
  )
  
  
  # ============================================================
  # 24. jitter
  # ============================================================
  
  pos_normal <- ggplot2::position_jitter(
    width = jitter_width,
    height = 0,
    seed = jitter_seed
  )
  
  pos_gene <- ggplot2::position_jitter(
    width = jitter_width,
    height = 0,
    seed = jitter_seed + 1
  )
  
  
  # ============================================================
  # 25. 开始绘图
  # ============================================================
  
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = group,
      y = score
    )
  )
  
  
  # ============================================================
  # 26. violin
  # ============================================================
  
  if (show_violin) {
    
    p <- p +
      
      ggplot2::geom_violin(
        
        ggplot2::aes(
          group = group
        ),
        
        fill = violin_fill,
        
        color = violin_color,
        
        alpha = violin_alpha,
        
        width = violin_width,
        
        linewidth = violin_linewidth,
        
        trim = violin_trim,
        
        scale = "width",
        
        show.legend = FALSE
      )
  }
  
  
  # ============================================================
  # 27. 普通细胞
  #
  # 统一使用 shape 21 + fill
  #
  # 这样主图与右侧图例完全使用同一套颜色系统
  # ============================================================
  
  if (nrow(normal_data) > 0) {
    
    if (point_border) {
      
      # --------------------------------------------------------
      # 有边框
      # --------------------------------------------------------
      
      p <- p +
        
        ggplot2::geom_point(
          
          data = normal_data,
          
          ggplot2::aes(
            fill = level
          ),
          
          shape = 21,
          
          color = point_border_color,
          
          stroke = point_border_width,
          
          size = point_size,
          
          alpha = point_alpha,
          
          position = pos_normal,
          
          show.legend = TRUE
        )
      
    } else {
      
      # --------------------------------------------------------
      # 无边框
      # --------------------------------------------------------
      
      p <- p +
        
        ggplot2::geom_point(
          
          data = normal_data,
          
          ggplot2::aes(
            fill = level
          ),
          
          shape = 21,
          
          color = NA,
          
          stroke = 0,
          
          size = point_size,
          
          alpha = point_alpha,
          
          position = pos_normal,
          
          show.legend = TRUE
        )
    }
  }
  
  
  # ============================================================
  # 28. dominant gene 特殊细胞
  # ============================================================
  
  if (nrow(gene_top_data) > 0) {
    
    if (gene_border) {
      
      # --------------------------------------------------------
      # 有边框三角形
      # --------------------------------------------------------
      
      p <- p +
        
        ggplot2::geom_point(
          
          data = gene_top_data,
          
          ggplot2::aes(
            shape = special_type
          ),
          
          fill = gene_color,
          
          color = gene_border_color,
          
          stroke = gene_border_width,
          
          size = gene_size,
          
          alpha = gene_alpha,
          
          position = pos_gene,
          
          show.legend = TRUE
        )
      
    } else {
      
      # --------------------------------------------------------
      # 无边框三角形
      # --------------------------------------------------------
      
      p <- p +
        
        ggplot2::geom_point(
          
          data = gene_top_data,
          
          ggplot2::aes(
            shape = special_type
          ),
          
          color = gene_color,
          
          size = gene_size,
          
          alpha = gene_alpha,
          
          position = pos_gene,
          
          show.legend = TRUE
        )
    }
  }
  
  
  # ============================================================
  # 29. infection level fill 颜色
  #
  # 左侧普通点和右侧图例都由这里控制
  # ============================================================
  
  p <- p +
    
    ggplot2::scale_fill_manual(
      
      values = colors_use,
      
      breaks = level_order,
      
      limits = level_order,
      
      drop = FALSE,
      
      name = level_col
    )
  
  
  # ============================================================
  # 30. dominant shape
  # ============================================================
  
  gene_shape_use <- if (
    gene_border
  ) {
    24
  } else {
    gene_shape
  }
  
  
  p <- p +
    
    ggplot2::scale_shape_manual(
      
      values = setNames(
        gene_shape_use,
        gene_legend_label
      ),
      
      name = NULL
    )
  
  
  # ============================================================
  # 31. X轴
  # ============================================================
  
  p <- p +
    
    ggplot2::scale_x_discrete(
      labels = x_labels,
      drop = TRUE
    )
  
  
  # ============================================================
  # 32. Y轴
  # ============================================================
  
  p <- p +
    
    ggplot2::scale_y_continuous(
      
      limits = c(
        0,
        final_y_max
      ),
      
      breaks = y_breaks,
      
      labels = scales::label_number(
        accuracy = NULL,
        big.mark = ","
      ),
      
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.03
        )
      )
    )
  
  
  # ============================================================
  # 33. 标题
  # ============================================================
  
  if (is.null(y_title)) {
    y_title <- score_col
  }
  
  if (is.null(x_title)) {
    x_title <- NULL
  }
  
  if (is.null(title)) {
    
    title <- paste0(
      gene,
      " dominant infected cells"
    )
  }
  
  p <- p +
    
    ggplot2::labs(
      x = x_title,
      y = y_title,
      title = title
    )
  
  
  # ============================================================
  # 34. 图例
  #
  # 关键：
  # 右侧普通细胞颜色按照 colors_use 强制设置
  # ============================================================
  
  if (point_border) {
    
    p <- p +
      
      ggplot2::guides(
        
        fill = ggplot2::guide_legend(
          
          order = 1,
          
          override.aes = list(
            
            size = legend_point_size,
            
            alpha = 1,
            
            shape = 21,
            
            fill = unname(
              colors_use
            ),
            
            color = point_border_color,
            
            stroke = point_border_width
          )
        ),
        
        
        shape = ggplot2::guide_legend(
          
          order = 2,
          
          override.aes = list(
            
            size = legend_gene_size,
            
            alpha = 1,
            
            shape = gene_shape_use,
            
            fill = if (gene_border) {
              gene_color
            } else {
              NA
            },
            
            color = if (gene_border) {
              gene_border_color
            } else {
              gene_color
            },
            
            stroke = if (gene_border) {
              gene_border_width
            } else {
              0
            }
          )
        )
      )
    
  } else {
    
    p <- p +
      
      ggplot2::guides(
        
        fill = ggplot2::guide_legend(
          
          order = 1,
          
          override.aes = list(
            
            size = legend_point_size,
            
            alpha = 1,
            
            shape = 21,
            
            fill = unname(
              colors_use
            ),
            
            color = NA,
            
            stroke = 0
          )
        ),
        
        
        shape = ggplot2::guide_legend(
          
          order = 2,
          
          override.aes = list(
            
            size = legend_gene_size,
            
            alpha = 1,
            
            shape = gene_shape_use,
            
            fill = if (gene_border) {
              gene_color
            } else {
              NA
            },
            
            color = if (gene_border) {
              gene_border_color
            } else {
              gene_color
            },
            
            stroke = if (gene_border) {
              gene_border_width
            } else {
              0
            }
          )
        )
      )
  }
  
  
  # ============================================================
  # 35. Theme
  # ============================================================
  
  p <- p +
    
    ggplot2::theme_classic(
      base_size = base_size
    ) +
    
    ggplot2::theme(
      
      axis.text.x = ggplot2::element_text(
        color = "black",
        size = base_size - 1,
        lineheight = 1.15
      ),
      
      axis.text.y = ggplot2::element_text(
        color = "black",
        size = base_size - 1
      ),
      
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = base_size
      ),
      
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = base_size
      ),
      
      axis.line = ggplot2::element_line(
        linewidth = 0.6,
        color = "black"
      ),
      
      axis.ticks = ggplot2::element_line(
        color = "black"
      ),
      
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      
      legend.title = ggplot2::element_text(
        size = legend_title_size
      ),
      
      legend.text = ggplot2::element_text(
        size = legend_text_size,
        color = "black"
      ),
      
      legend.key.height = grid::unit(
        1.4,
        "lines"
      ),
      
      legend.key.width = grid::unit(
        1.4,
        "lines"
      ),
      
      legend.spacing.y = grid::unit(
        0.3,
        "cm"
      ),
      
      plot.margin = ggplot2::margin(
        t = 10,
        r = 15,
        b = 15,
        l = 10
      )
    )
  
  
  # ============================================================
  # 36. 返回
  # ============================================================
  
  return(
    list(
      
      plot = p,
      
      data = plot_data,
      
      gene_top_cells = gene_top_data,
      
      group_stats = group_stats,
      
      y_breaks = y_breaks,
      
      y_max = final_y_max,
      
      infection_levels = infection_levels,
      
      min_gene_count = min_gene_count,
      
      min_total_viral_count =
        min_total_viral_count,
      
      point_border = point_border,
      
      gene_border = gene_border,
      
      legend_point_size =
        legend_point_size,
      
      legend_gene_size =
        legend_gene_size,
      
      colors_use =
        colors_use,
      
      level_order =
        level_order,
      
      virus_gene_present =
        virus_gene_present,
      
      virus_gene_missing =
        virus_gene_missing
    )
  )
}


