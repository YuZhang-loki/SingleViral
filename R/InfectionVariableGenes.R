#' Differential expression between infected vs uninfected within selected groups,
#' with publication-ready volcano plot
#'
#' @description
#' Filters cells by timepoints and manual labels, converts an infection-level
#' column into a binary infection status (case vs control), runs
#' \code{Seurat::FindMarkers()} between the two groups, removes optional
#' mitochondrial/viral genes, classifies genes as Up/Down/No change, and
#' produces a volcano plot. Labeled genes can be specified manually
#' (\code{label_genes}) or chosen as the top-N by smallest adjusted p-value
#' (\code{top_n_pval}). Results and plots can be saved to disk.
#'
#' @param obj A \strong{Seurat} object.
#' @param timepoint_col Character. Meta-data column name for timepoints.
#' @param manual_col Character. Meta-data column name for manual cell-type labels.
#' @param level_col Character. Meta-data column name storing integer infection levels.
#' @param time_values Character vector of timepoints to include (subset).
#' @param manual_include Character vector of manual labels (cell types) to include.
#' @param case_levels Integer vector. Levels in \code{level_col} treated as \emph{infected}.
#' @param control_levels Integer vector. Levels in \code{level_col} treated as \emph{uninfected}.
#' @param time_levels Optional character vector giving the factor order for timepoints
#'   (used only for consistent metadata ordering; set \code{NULL} to skip).
#' @param virus_genes Character vector of gene symbols to exclude from DE (e.g., viral genes).
#' @param mt_pattern Character regex for mitochondrial genes to exclude (e.g., \code{"^MT-"}).
#'
#' @param test.use, min.pct, logfc.threshold, only.pos, pseudocount.use, assay
#'  Passed to \code{Seurat::FindMarkers()}.
#'
#' @param p_cut Numeric FDR threshold for volcano classification (default 0.05).
#' @param log2fc_cut Numeric absolute log2FC threshold for volcano classification (default 0.5).
#' @param label_genes Optional character vector of genes to always label on the plot.
#' @param top_n_pval Integer. Label the top-N genes with smallest adjusted p-values
#'   (ignored if \code{label_genes} is provided).
#' @param point_size Numeric point size in volcano scatter.
#' @param alpha_nochange Alpha for points classified as "No change".
#' @param color_up,color_down,color_nochange,line_color Colors for volcano plot.
#' @param repel_max_overlaps Integer passed to \code{ggrepel::geom_text_repel(max.overlaps=)}.
#'
#' @param save_csv File name for writing the DE table (CSV). Set \code{NULL} to skip.
#' @param save_plot File name for saving the volcano PDF. Set \code{NULL} to skip.
#' @param save_dir Output directory. Created if missing.
#' @param width,height Numeric plot size in inches for saved figures.
#' @param plot_title Character. Title for the volcano plot.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{n_case}, \code{n_control} — numbers of infected/uninfected cells used
#'   \item \code{cells_used} — character vector of cell barcodes
#'   \item \code{deg_table} — tibble of DE results with annotated columns
#'   \item \code{volcano} — a \code{ggplot} object
#'   \item \code{object} — the filtered Seurat object used for DE
#' }
#'
#' @examples
#' \dontrun{
#' res <- InfectionVariableGenes(
#'   obj            = data,
#'   time_values    = c("24hpi","48hpi"),
#'   manual_include = c("CD4+Tcells"),
#'   case_levels    = c(1,2,3,4),
#'   control_levels = 0,
#'   label_genes    = c("HSPA8","SDF2L1"),
#'   top_n_pval     = 30,
#'   save_dir       = "results/DE_Tcells"
#' )
#' res$volcano
#' }
#'
#' @seealso \code{\link[Seurat]{FindMarkers}}, \code{\link[ggrepel]{geom_text_repel}}
#'
#' @importFrom Seurat FindMarkers Idents
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate case_when arrange slice_head pull filter %>%
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline scale_color_manual
#'   theme_bw theme element_blank element_text xlab ylab ggtitle ggsave
#' @importFrom ggrepel geom_text_repel
#' @importFrom readr write_csv
#' @export
InfectionVariableGenes <- function(
    obj,
    # —— Metadata columns —— #
    timepoint_col   = "timepoint",
    manual_col      = "manual_label",
    level_col       = "infection_level",
    # —— Selection filters —— #
    time_values     = c("24hpi"),
    manual_include  = c("CD4+Tcells"),
    case_levels     = c(1,2,3,4),
    control_levels  = c(0),
    time_levels     = c("control","24hpi","48hpi","72hpi","96hpi"),
    # —— Optional filters —— #
    virus_genes     = character(0),
    mt_pattern      = "^MT-",
    # —— FindMarkers parameters —— #
    test.use        = "wilcox",
    min.pct         = 0.1,
    logfc.threshold = 0,
    only.pos        = FALSE,
    pseudocount.use = 1,
    assay           = NULL,
    # —— Volcano thresholds and labeling —— #
    p_cut           = 0.05,
    log2fc_cut      = 0.5,
    label_genes     = NULL,
    top_n_pval      = 30,
    point_size      = 1.5,
    alpha_nochange  = 0.6,
    color_up        = "#F45B69",
    color_down      = "#4EA8DE",
    color_nochange  = "gray70",
    line_color      = "#999999",
    repel_max_overlaps = 100,
    # —— Saving options —— #
    save_csv   = "DEG_result.csv",
    save_plot  = "Volcano_plot.pdf",
    save_dir   = "./results",
    width      = 7,
    height     = 6,
    # —— Volcano plot title —— #
    plot_title = "Volcano Plot: Infected vs Uninfected"
){
  stopifnot(inherits(obj, "Seurat"))
  md <- obj@meta.data
  need <- c(timepoint_col, manual_col, level_col)
  miss <- setdiff(need, colnames(md))
  if (length(miss) > 0)
    stop("Missing metadata columns: ", paste(miss, collapse = ", "))

  ## Step 1: Filter by timepoints
  cells_tp <- rownames(md[ md[[timepoint_col]] %in% time_values, , drop = FALSE ])
  if (length(cells_tp) == 0)
    stop("No cells found for the given 'time_values'.")
  obj_tp  <- subset(obj, cells = cells_tp)

  ## Step 2: Filter by cell types (manual labels)
  md2 <- obj_tp@meta.data
  cells_group <- rownames(md2[ md2[[manual_col]] %in% manual_include, , drop = FALSE ])
  if (length(cells_group) == 0)
    stop("No cells found for the given 'manual_include'.")
  obj_sel <- subset(obj_tp, cells = cells_group)

  ## Step 3: Convert infection_level to infection_status
  lev <- obj_sel@meta.data[[level_col]]
  if (all(grepl("^\\d+$", as.character(lev)))) lev <- as.integer(as.character(lev))
  inf_stat <- ifelse(lev %in% case_levels, "infected",
                     ifelse(lev %in% control_levels, "uninfected", NA))
  obj_sel$infection_status <- inf_stat

  keep_cells <- colnames(obj_sel)[!is.na(obj_sel$infection_status)]
  if (length(keep_cells) == 0)
    stop("No valid cells remain after assigning infection_status.")
  obj_sel <- subset(obj_sel, cells = keep_cells)

  ## Step 4: Optional ordering of timepoints
  if (!is.null(time_levels)) {
    tp <- obj_sel@meta.data[[timepoint_col]]
    obj_sel@meta.data[[timepoint_col]] <- factor(tp, levels = intersect(time_levels, unique(as.character(tp))))
  }

  ## Step 5: Differential expression (infected vs uninfected)
  Seurat::Idents(obj_sel) <- "infection_status"
  deg <- Seurat::FindMarkers(
    object          = obj_sel,
    ident.1         = "infected",
    ident.2         = "uninfected",
    test.use        = test.use,
    min.pct         = min.pct,
    logfc.threshold = logfc.threshold,
    only.pos        = only.pos,
    pseudocount.use = pseudocount.use,
    assay           = assay
  ) %>% tibble::rownames_to_column("gene")

  if (!"p_val_adj" %in% colnames(deg)) {
    if ("p_val" %in% colnames(deg)) {
      deg$p_val_adj <- p.adjust(deg$p_val, method = "BH")
    } else {
      stop("Differential expression results are missing 'p_val' or 'p_val_adj'.")
    }
  }

  ## Step 6: Remove viral or mitochondrial genes (optional)
  if (length(virus_genes) > 0)
    deg <- deg[ !(deg$gene %in% virus_genes), , drop = FALSE ]
  if (!is.na(mt_pattern))
    deg <- deg[ !grepl(mt_pattern, deg$gene), , drop = FALSE ]

  ## Step 7: Classify Up / Down / No change
  deg <- dplyr::mutate(
    .data = deg,
    neg_log10_p = -log10(p_val_adj),
    change = dplyr::case_when(
      p_val_adj < p_cut & avg_log2FC >  log2fc_cut ~ "Up",
      p_val_adj < p_cut & avg_log2FC < -log2fc_cut ~ "Down",
      TRUE                                         ~ "No change"
    )
  )

  ## Step 8: Determine labeled genes
  if (!is.null(label_genes)) {
    lab_set <- intersect(label_genes, deg$gene)
  } else if (!is.null(top_n_pval) && top_n_pval > 0) {
    top_n_pval <- as.integer(top_n_pval)
    if (is.na(top_n_pval) || top_n_pval < 0) top_n_pval <- 0
    k <- min(top_n_pval, nrow(deg))
    lab_set <- deg %>% dplyr::arrange(p_val_adj) %>% dplyr::slice_head(n = k) %>% dplyr::pull(gene)
  } else {
    lab_set <- character(0)
  }
  deg$label_flag <- deg$gene %in% lab_set

  ## Step 9: Volcano plot
  p_vol <- ggplot2::ggplot(deg, ggplot2::aes(x = avg_log2FC, y = neg_log10_p)) +
    ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = line_color) +
    ggplot2::geom_vline(xintercept = c(-log2fc_cut, log2fc_cut), linetype = "dashed", color = line_color) +
    ggplot2::geom_point(ggplot2::aes(color = change),
                        size = point_size,
                        alpha = ifelse(deg$change == "No change", alpha_nochange, 1)) +
    ggplot2::scale_color_manual(values = c("Up" = color_up, "Down" = color_down, "No change" = color_nochange)) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = "right",
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 16)) +
    ggplot2::xlab("Log2FC") + ggplot2::ylab("-Log10(FDR)") +
    ggplot2::ggtitle(plot_title)

  # Add labels only for Up/Down (skip "No change")
  if (any(deg$label_flag)) {
    label_df <- dplyr::filter(deg, label_flag & change != "No change")
    if (nrow(label_df) > 0) {
      p_vol <- p_vol +
        ggrepel::geom_text_repel(
          data = label_df,
          ggplot2::aes(label = gene, color = change),
          size = 3,
          max.overlaps = repel_max_overlaps,
          box.padding = 0.5,
          show.legend = FALSE
        )
    }
  }

  ## Step 10: Save outputs
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  if (!is.null(save_csv))
    readr::write_csv(deg, file.path(save_dir, save_csv))
  if (!is.null(save_plot))
    ggplot2::ggsave(file.path(save_dir, save_plot), p_vol, width = width, height = height, useDingbats = FALSE)

  list(
    n_case     = sum(obj_sel$infection_status == "infected"),
    n_control  = sum(obj_sel$infection_status == "uninfected"),
    cells_used = colnames(obj_sel),
    deg_table  = deg,
    volcano    = p_vol,
    object     = obj_sel
  )
}
