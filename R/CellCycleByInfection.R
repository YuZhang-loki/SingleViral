#' @title Cell-cycle Scoring and Visualization by Infection Level
#'
#' @description
#' Performs Seurat-based cell-cycle scoring (S and G2M phases) and visualizes
#' both *cell-cycle scores* and *phase proportions* across different infection
#' levels (or any grouping variable) in a Seurat object.
#'
#' The function generates:
#' 1. A violin plot of S.Score and G2M.Score by group.
#' 2. A stacked bar plot showing proportions of G1, S, and G2M phases by group.
#'
#' @param obj A Seurat object.
#' @param level_col Metadata column used for grouping (e.g., "infection_level").
#' @param species One of "human", "mouse", or "other".
#'        When "human" or "mouse", cell-cycle gene sets are read from the
#'        package inst/ directory. When "other", user must supply s.genes and g2m.genes.
#' @param s.genes S phase gene vector (required only when species = "other").
#' @param g2m.genes G2M phase gene vector (required only when species = "other").
#' @param cell_type Optional metadata column name used for subsetting specific cell types.
#' @param celltype_include Optional character vector of cell types to retain.
#'        If either argument is NULL, subsetting is skipped.
#' @param assay Assay to use (default = NULL, meaning DefaultAssay(obj)).
#' @param layer Data layer to use: "data" (default) or "counts".
#' @param title_violin Title for violin plot.
#' @param title_phase Title for stacked bar plot.
#' @param palette Colors for G1, S, and G2M.
#' @param lvl_colors Custom colors for violin plot groups.
#' @param save_dir Folder to save plots. NULL = do not save.
#' @param file_violin Filename of violin plot.
#' @param file_phase Filename of stacked phase plot.
#' @param width,height Plot dimensions.
#' @param show_points Whether to overlay points on violin plot.
#'
#' @return A list containing violin_plot, phase_plot, and phase_table.
#' @export
CellCycleByInfection <- function(
    obj,
    level_col   = "infection_level",
    species     = c("human", "mouse", "other"),
    s.genes     = NULL,
    g2m.genes   = NULL,
    cell_type   = NULL,
    celltype_include = NULL,
    assay       = NULL,
    layer        = "data",
    title_violin = "Cell-cycle scores by infection level",
    title_phase  = "Cell-cycle phase proportions by infection level",
    palette      = c("G1" = "#6C757D", "S" = "#4EA8DE", "G2M" = "#F45B69"),
    lvl_colors   = NULL,
    save_dir     = NULL,
    file_violin  = "cellcycle_scores_violin.pdf",
    file_phase   = "cellcycle_phase_proportion.pdf",
    width        = 8, height = 6,
    show_points  = TRUE
){
  stopifnot(inherits(obj, "Seurat"))
  md <- obj@meta.data
  if (!(level_col %in% colnames(md))) {
    stop("Column '", level_col, "' not found in obj@meta.data.")
  }

  ## -------------------------------
  ## 1) Optional subsetting by cell type
  ## -------------------------------
  if (!is.null(cell_type) && !is.null(celltype_include)) {
    if (!(cell_type %in% colnames(md))) {
      stop("Column '", cell_type, "' not found in obj@meta.data.")
    }
    keep_cells <- rownames(md)[md[[cell_type]] %in% celltype_include]
    if (length(keep_cells) == 0) {
      stop("No cells found for specified 'celltype_include' values.")
    }
    obj <- subset(obj, cells = keep_cells)
    md  <- obj@meta.data
  }

  ## -------------------------------
  ## 2) Load species-specific gene sets using readLines()
  ## -------------------------------
  species <- match.arg(species)

  if (species %in% c("human", "mouse")) {
    pkg <- utils::packageName()
    if (is.null(pkg) || is.na(pkg)) pkg <- "SingleViral"  # your package name

    if (species == "human") {
      s_file <- system.file("s_genes_human.txt",  package = pkg)
      g_file <- system.file("g2m_genes_human.txt", package = pkg)
    } else {  # mouse
      s_file <- system.file("s_genes_mouse.txt",  package = pkg)
      g_file <- system.file("g2m_genes_mouse.txt", package = pkg)
    }

    if (s_file == "" || g_file == "") {
      stop("Cannot find built-in ", species, " gene files in package '", pkg, "'.")
    }

    s.genes   <- unique(readLines(s_file, warn = FALSE))
    g2m.genes <- unique(readLines(g_file, warn = FALSE))

    s.genes   <- s.genes[nzchar(s.genes)]
    g2m.genes <- g2m.genes[nzchar(g2m.genes)]

  } else {
    if (is.null(s.genes) || is.null(g2m.genes)) {
      stop("For species='other', both s.genes and g2m.genes must be supplied.")
    }
  }

  ## -------------------------------
  ## 3) Select assay and layer
  ## -------------------------------
  if (!is.null(assay)) DefaultAssay(obj) <- assay
  dat <- GetAssayData(obj, layer = layer)
  if (nrow(dat) == 0) stop("Expression matrix in layer '", layer, "' is empty.")

  ## -------------------------------
  ## 4) Run Seurat cell-cycle scoring
  ## -------------------------------
  obj <- Seurat::CellCycleScoring(
    object       = obj,
    s.features   = s.genes,
    g2m.features = g2m.genes,
    set.ident    = FALSE
  )

  cols_use <- if (!is.null(lvl_colors)) lvl_colors else NULL

  ## Violin plot
  p_violin <- Seurat::VlnPlot(
    obj,
    features = c("S.Score", "G2M.Score"),
    group.by = level_col,
    pt.size  = if (show_points) 0.2 else 0,
    cols     = cols_use
  ) +
    patchwork::plot_annotation(
      title = title_violin,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold")
      )
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10))

  ## -------------------------------
  ## 5) Calculate phase proportions
  ## -------------------------------
  df_phase <- obj@meta.data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(
      Phase = as.character(Phase),
      !!rlang::sym(level_col) := as.character(!!rlang::sym(level_col))
    ) %>%
    dplyr::count(!!rlang::sym(level_col), Phase, name = "n") %>%
    dplyr::group_by(!!rlang::sym(level_col)) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()

  df_phase$Phase <- factor(df_phase$Phase, levels = c("G1", "S", "G2M"))

  ## Stacked bar plot
  p_phase <- ggplot(df_phase, aes(x = !!rlang::sym(level_col), y = prop, fill = Phase)) +
    geom_col(width = 0.7, color = "white") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = palette) +
    labs(x = level_col, y = "Proportion", title = title_phase) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 20))

  ## -------------------------------
  ## 6) Save plots if requested
  ## -------------------------------
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggsave(file.path(save_dir, file_violin), p_violin,
           width = width, height = height, useDingbats = FALSE)
    ggsave(file.path(save_dir, file_phase),  p_phase,
           width = width, height = height, useDingbats = FALSE)
  }

  ## -------------------------------
  ## Return results
  ## -------------------------------
  list(
    violin_plot = p_violin,
    phase_plot  = p_phase,
    phase_table = df_phase
  )
}
