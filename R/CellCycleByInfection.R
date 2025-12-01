#' @title Cell-cycle Scoring and Visualization by Infection Level
#'
#' @description
#' Performs Seurat-based cell-cycle scoring (S, G2M phases) and visualizes
#' both the *cell-cycle scores* and *phase proportions* across different
#' infection levels (or other grouping variables) within a Seurat object.
#'
#' The function produces:
#' 1. A violin plot showing S.Score and G2M.Score distributions by group.
#' 2. A stacked bar plot showing the proportion of cells in G1, S, and G2M phases by group.
#'
#' @param obj A Seurat object containing single-cell expression data.
#' @param level_col Character string specifying the metadata column used for grouping
#'        (e.g., `"infection_level"` or `"manual_label"`).
#' @param s.genes Character vector of S phase gene symbols.
#'        Defaults to `Seurat::cc.genes.updated.2019$s.genes`.
#' @param g2m.genes Character vector of G2/M phase gene symbols.
#'        Defaults to `Seurat::cc.genes.updated.2019$g2m.genes`.
#' @param assay Character string specifying which assay to use
#'        (e.g., `"RNA"`, `"SCT"`). Default: `NULL`, uses `DefaultAssay(obj)`.
#' @param slot Character string specifying which assay slot to use.
#'        Options: `"data"` (default, normalized values) or `"counts"` (raw counts).
#' @param title_violin Character string, title for the violin plot.
#'        Default: `"Cell-cycle scores by infection level"`.
#' @param title_phase Character string, title for the stacked phase proportion plot.
#'        Default: `"Cell-cycle phase proportions by infection level"`.
#' @param palette Named character vector specifying fill colors for cell-cycle phases.
#'        Default: `c("G1"="#6C757D", "S"="#4EA8DE", "G2M"="#F45B69")`.
#' @param lvl_colors Optional character vector specifying custom colors for infection-level groups
#'        (applied to violin plots only). Default: `NULL`.
#' @param save_dir Optional character string for output directory to save plots (e.g., `"results/cellcycle"`).
#'        If `NULL`, plots are not saved.
#' @param file_violin Character string specifying the filename for the violin plot PDF.
#'        Default: `"cellcycle_scores_violin.pdf"`.
#' @param file_phase Character string specifying the filename for the phase proportion plot PDF.
#'        Default: `"cellcycle_phase_proportion.pdf"`.
#' @param width Numeric, plot width in inches for saved figures. Default: `7`.
#' @param height Numeric, plot height in inches for saved figures. Default: `5`.
#' @param show_points Logical, whether to overlay points on the violin plot. Default: `TRUE`.
#'
#' @details
#' This function wraps Seurat's \code{\link[Seurat]{CellCycleScoring}} and
#' automatically generates cell-cycle scores (`S.Score`, `G2M.Score`)
#' and categorical phases (`Phase`), grouped by a user-defined infection level
#' or experimental condition.
#'
#' The function then visualizes these results in:
#' \itemize{
#'   \item A violin plot (`S.Score`, `G2M.Score`) by infection level.
#'   \item A stacked bar chart of G1/S/G2M phase proportions by infection level.
#' }
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{`object`}{The input Seurat object with updated metadata columns:
#'         `S.Score`, `G2M.Score`, and `Phase`.}
#'   \item{`violin_plot`}{A ggplot object showing the S and G2M scores per group.}
#'   \item{`phase_plot`}{A ggplot object showing phase proportions (G1/S/G2M).}
#'   \item{`phase_table`}{A data frame containing counts and proportions of each phase by group.}
#' }
#'
#' @seealso
#' \code{\link[Seurat]{CellCycleScoring}}, \code{\link[ggplot2]{ggplot}},
#' \code{\link[patchwork]{plot_annotation}}
#'
#' @import Seurat ggplot2 dplyr patchwork tibble scales rlang
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' res <- CellCycleByInfection(
#'   obj = seurat_obj,
#'   level_col = "infection_level",
#'   save_dir = "results/cellcycle",
#'   lvl_colors = c("#7FB3D5","#5DADE2","#48C9B0","#F5B041","#D35400")
#' )
#'
#' # Access components
#' res$violin_plot
#' res$phase_plot
#' res$phase_table
#' }
#'
#' @export
CellCycleByInfection <- function(
    obj,
    level_col   = "infection_level",
    s.genes     = Seurat::cc.genes.updated.2019$s.genes,
    g2m.genes   = Seurat::cc.genes.updated.2019$g2m.genes,
    assay       = NULL,
    slot        = "data",
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

  if (!is.null(assay)) DefaultAssay(obj) <- assay
  if (slot == "counts") dat <- GetAssayData(obj, slot = "counts") else dat <- GetAssayData(obj, slot = "data")
  if (nrow(dat) == 0) stop("No expression matrix available in slot '", slot, "'.")

  obj <- Seurat::CellCycleScoring(
    object    = obj,
    s.features   = s.genes,
    g2m.features = g2m.genes,
    set.ident    = FALSE
  )

  cols_use <- if (!is.null(lvl_colors)) lvl_colors else NULL

  p_violin <- Seurat::VlnPlot(
    obj,
    features = c("S.Score", "G2M.Score"),
    group.by = level_col,
    pt.size  = if (show_points) 0.2 else 0,
    cols     = cols_use
  )

  p_violin <- p_violin +
    patchwork::plot_annotation(
      title = title_violin,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold")
      )
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, size = 10)
    )

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

  p_phase <- ggplot(df_phase, aes(x = !!rlang::sym(level_col), y = prop, fill = Phase)) +
    geom_col(width = 0.7, color = "white") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = palette) +
    xlab(level_col) + ylab("Proportion") +
    ggtitle(title_phase) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 20))

  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, file_violin), p_violin, width = width, height = height, useDingbats = FALSE)
    ggplot2::ggsave(file.path(save_dir, file_phase),  p_phase,  width = width, height = height, useDingbats = FALSE)
  }

  list(
    object = obj,
    violin_plot = p_violin,
    phase_plot  = p_phase,
    phase_table = df_phase
  )
}
