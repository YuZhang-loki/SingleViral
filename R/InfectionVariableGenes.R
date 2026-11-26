#' Identify Infection-Associated Host Genes Using Differential Expression
#' and Viral-Load Correlation
#'
#' @description
#' `InfectionVariableGenes()` identifies and prioritizes host genes
#' associated with viral transcriptional burden by integrating
#' differential expression analysis, gene expression–viral-load
#' correlation analysis, and a composite rank-based scoring strategy.
#'
#' The function compares infected and uninfected cells within selected
#' experimental conditions and cell populations, evaluates correlations
#' between differentially expressed genes and viral transcriptional
#' burden, and ranks candidate genes according to their differential
#' expression significance, correlation strength, and expression
#' prevalence.
#'
#' The function also generates differential expression volcano plots,
#' correlation–significance plots, and infection-level gene expression
#' trend plots to facilitate the identification and interpretation of
#' infection-associated host factors.
#'
#' @details
#' The analysis consists of six major steps.
#'
#' **1. Cell selection and infection status assignment**
#'
#' Cells are selected according to the experimental conditions specified
#' by `group_include` and the cell populations specified by
#' `celltype_include`.
#'
#' Selected cells are classified into two groups:
#'
#' - Infected: cells with infection levels specified by `case_levels`.
#' - Uninfected: cells with infection levels specified by
#'   `control_levels`.
#'
#' Cells outside these levels are excluded from subsequent analyses.
#'
#' **2. Differential expression analysis**
#'
#' Differentially expressed genes (DEGs) between infected and
#' uninfected cells are identified using `Seurat::FindMarkers()`.
#'
#' Genes are classified as:
#'
#' - Up: adjusted P value < `deg_p_cut` and
#'   log2 fold change > `deg_log2fc_cut`.
#' - Down: adjusted P value < `deg_p_cut` and
#'   log2 fold change < -`deg_log2fc_cut`.
#' - No change: genes not meeting either criterion.
#'
#' Viral genes and mitochondrial genes can be excluded from the
#' analysis using `virus_genes`, `mt_genes`, and `mt_pattern`.
#'
#' **3. Viral-load correlation analysis**
#'
#' For genes classified as significantly upregulated or downregulated,
#' correlations between gene expression and viral transcriptional burden
#' are calculated using either Spearman or Pearson correlation.
#'
#' Correlation analysis is performed using cells belonging to the
#' infection levels specified by `correlation_levels`.
#'
#' Correlation P values are adjusted using the Benjamini–Hochberg
#' procedure.
#'
#' Genes are classified as:
#'
#' - Positive: adjusted correlation P value < `cor_p_cut` and
#'   correlation coefficient > `cor_cut`.
#' - Negative: adjusted correlation P value < `cor_p_cut` and
#'   correlation coefficient < -`cor_cut`.
#' - No correlation: genes not meeting either criterion.
#'
#' **4. Integrated gene classification**
#'
#' Differential expression and viral-load correlation results are
#' integrated to classify genes into four categories:
#'
#' - Up-Positive: upregulated genes positively correlated with viral load.
#' - Up-Negative: upregulated genes negatively correlated with viral load.
#' - Down-Positive: downregulated genes positively correlated with viral load.
#' - Down-Negative: downregulated genes negatively correlated with viral load.
#'
#' Genes that do not satisfy both differential expression and
#' correlation criteria are classified as "Not significant".
#'
#' **5. Composite gene ranking**
#'
#' Genes satisfying both differential expression and correlation
#' thresholds are prioritized using an equal-weight percentile-ranking
#' strategy based on three indicators:
#'
#' - Correlation strength: absolute correlation coefficient, `abs(r)`.
#' - Differential expression significance: negative adjusted DEG
#'   P value, `-p_val_adj`.
#' - Expression prevalence: `pct.1` for upregulated genes and
#'   `pct.2` for downregulated genes.
#'
#' Each indicator is transformed into a percentile rank using
#' `dplyr::percent_rank()`.
#'
#' The composite ranking score is calculated as:
#'
#' `rank_score = (rho_rank + deg_rank + pct_rank) / 3`
#'
#' where:
#'
#' - `rho_rank` is the percentile rank of `abs(r)`.
#' - `deg_rank` is the percentile rank of `-p_val_adj`.
#' - `pct_rank` is the percentile rank of expression prevalence.
#'
#' Higher scores indicate stronger viral-load associations,
#' greater differential expression significance, and higher
#' expression prevalence.
#'
#' All significantly correlated DEGs are ranked together,
#' regardless of correlation direction.
#'
#' Genes are subsequently ranked within positive and negative
#' correlation groups to identify top candidates.
#'
#' **6. Visualization**
#'
#' Three types of visualizations are generated:
#'
#' - DEG volcano plot: displays log2 fold changes against
#'   -log10(adjusted DEG P values).
#' - Correlation–significance plot: displays viral-load correlation
#'   coefficients against -log10(adjusted DEG P values).
#' - Infection-level expression trends: displays gene expression
#'   across infection levels for top positively and negatively
#'   correlated genes.
#'
#' In the correlation–significance plot, significant positive
#' correlations are represented by a red gradient based on `pct.1`,
#' whereas significant negative correlations are represented by
#' a blue gradient based on `pct.2`.
#'
#' Gene expression trend plots display mean expression across
#' infection levels, with optional standard error ribbons,
#' violin distributions, and individual cell expression values.
#'
#' @param obj
#' A Seurat object containing single-cell gene expression data,
#' infection-level annotations, and viral transcriptional burden scores.
#'
#' @param group_col
#' Character string specifying the metadata column containing
#' experimental conditions or infection time points.
#' Default is `"timepoint"`.
#'
#' @param cell_type
#' Character string specifying the metadata column containing
#' cell-type annotations.
#' Default is `"manual_label"`.
#'
#' @param level_col
#' Character string specifying the metadata column containing
#' infection-level annotations.
#' Default is `"infection_level"`.
#'
#' @param viral_col
#' Character string specifying the metadata column containing
#' viral transcriptional burden scores.
#' Default is `"viral_load"`.
#'
#' @param group_include
#' Character vector specifying the experimental conditions
#' included in the analysis.
#' If `NULL`, all conditions are included.
#' Default is `"24hpi"`.
#'
#' @param celltype_include
#' Character vector specifying the cell types included in the analysis.
#' If `NULL`, all cell types are included.
#' Default is `"CD4+Tcells"`.
#'
#' @param case_levels
#' Numeric or character vector specifying infection levels
#' assigned to the infected group.
#' Default is `1:4`.
#'
#' @param control_levels
#' Numeric or character vector specifying infection levels
#' assigned to the uninfected group.
#' Default is `0`.
#'
#' @param group_levels
#' Optional character vector specifying the preferred ordering
#' of experimental conditions.
#' Default is
#' `c("control", "24hpi", "48hpi", "72hpi", "96hpi")`.
#' If `NULL`, the existing metadata values are retained.
#'
#' @param correlation_levels
#' Numeric or character vector specifying the infection levels
#' included in viral-load correlation analysis.
#' Default is `0:4`.
#'
#' @param min_cells
#' Minimum number of cells required for correlation analysis.
#' Default is `50`.
#' Must be at least 3.
#'
#' @param virus_genes
#' Character vector containing viral genes to exclude from
#' differential expression and correlation analyses.
#' Default is `character(0)`.
#'
#' @param mt_genes
#' Optional character vector containing mitochondrial genes
#' to exclude from the analysis.
#' Default is `NULL`.
#'
#' @param mt_pattern
#' Regular expression used to identify and exclude
#' mitochondrial genes.
#' Default is `"^MT-"`.
#' Set to `NULL` to disable pattern-based exclusion.
#'
#' @param test.use
#' Statistical test used by `Seurat::FindMarkers()`
#' for differential expression analysis.
#' Default is `"wilcox"`.
#'
#' @param min.pct
#' Minimum fraction of cells expressing a gene in at least
#' one comparison group for differential expression testing.
#' Default is `0.1`.
#'
#' @param logfc.threshold
#' Minimum absolute log fold-change threshold used for
#' gene filtering in `Seurat::FindMarkers()`.
#' Default is `0`.
#'
#' @param only.pos
#' Logical value indicating whether only positively
#' differentially expressed genes should be returned.
#' Default is `FALSE`.
#'
#' @param pseudocount.use
#' Pseudocount used in fold-change calculations performed
#' by `Seurat::FindMarkers()`.
#' Default is `1`.
#'
#' @param assay
#' Character string specifying the Seurat assay used
#' for gene expression analysis.
#' Default is `"RNA"`.
#'
#' @param layer
#' Character string specifying the Seurat assay layer
#' used for correlation analysis and expression trend plots.
#' Default is `"data"`.
#'
#' @param deg_p_cut
#' Adjusted P-value threshold for identifying significant DEGs.
#' Default is `0.1`.
#'
#' @param deg_log2fc_cut
#' Absolute log2 fold-change threshold for classifying
#' significantly upregulated and downregulated genes.
#' Default is `0.5`.
#'
#' @param cor_method
#' Correlation method used to evaluate associations between
#' gene expression and viral transcriptional burden.
#' Available options are `"spearman"` and `"pearson"`.
#' Default is `"spearman"`.
#'
#' @param cor_p_cut
#' Adjusted P-value threshold for identifying significant
#' viral-load correlations.
#' Default is `0.05`.
#'
#' @param cor_cut
#' Minimum absolute correlation coefficient required
#' to classify a gene as positively or negatively correlated.
#' Default is `0.3`.
#'
#' @param label_genes
#' Optional character vector specifying genes to label
#' in the DEG volcano and correlation–significance plots.
#' If `NULL`, genes are selected automatically.
#' Default is `NULL`.
#'
#' @param top_n_pval
#' Maximum number of significant DEGs selected for
#' automatic labeling in the DEG volcano plot.
#' Genes are selected according to adjusted P values.
#' Default is `30`.
#'
#' @param point_size
#' Point size used in the DEG volcano plot.
#' Default is `1.5`.
#'
#' @param alpha_nochange
#' Transparency of nonsignificant genes in the
#' DEG volcano and correlation–significance plots.
#' Default is `0.6`.
#'
#' @param color_up
#' Color used for upregulated genes in the DEG volcano plot.
#' Default is `"#F45B69"`.
#'
#' @param color_down
#' Color used for downregulated genes in the DEG volcano plot.
#' Default is `"#4EA8DE"`.
#'
#' @param color_nochange
#' Color used for nonsignificant genes.
#' Default is `"gray70"`.
#'
#' @param line_color
#' Color of significance and effect-size threshold lines.
#' Default is `"#999999"`.
#'
#' @param repel_max_overlaps
#' Maximum number of overlapping labels allowed by
#' `ggrepel::geom_text_repel()`.
#' Default is `100`.
#'
#' @param cor_point_size
#' Point size used in the correlation–significance plot.
#' Default is `2`.
#'
#' @param cor_label_top_n
#' Requested number of top-ranked genes to label
#' in the correlation–significance plot.
#'
#' When `label_genes = NULL`, labels are selected from
#' positive and negative correlation groups separately.
#'
#' Default is `20`.
#'
#' @param trend_top_n
#' Maximum number of top-ranked genes selected separately
#' from positive and negative correlation groups
#' for expression trend visualization.
#' Default is `10`.
#'
#' @param trend_genes
#' Optional character vector used to restrict the
#' automatically selected top-ranked genes displayed
#' in expression trend plots.
#'
#' This parameter filters the selected top-ranked genes
#' rather than replacing the ranking procedure.
#' Default is `NULL`.
#'
#' @param summarize_fun
#' Function used to summarize gene expression values
#' within each infection level.
#' Default is `mean`.
#'
#' The function must accept an `na.rm` argument.
#'
#' @param show_se
#' Logical value indicating whether standard error
#' ribbons should be displayed in expression trend plots.
#' Default is `TRUE`.
#'
#' @param trend_point_size
#' Point size used for summarized expression values
#' in expression trend plots.
#' Default is `2.2`.
#'
#' @param trend_line_width
#' Line width used for summarized expression trends.
#' Default is `1.2`.
#'
#' @param violin_alpha
#' Transparency of violin distributions in expression
#' trend plots.
#' Default is `0.25`.
#'
#' @param show_violin
#' Logical value indicating whether violin distributions
#' should be displayed in expression trend plots.
#' Default is `TRUE`.
#'
#' @param violin_width
#' Width of violin distributions.
#' Default is `0.85`.
#'
#' @param show_cell_points
#' Logical value indicating whether individual cell
#' expression values should be displayed.
#' Default is `TRUE`.
#'
#' @param cell_point_alpha
#' Transparency of individual cell expression points.
#' Default is `0.15`.
#'
#' @param cell_point_size
#' Point size of individual cell expression values.
#' Default is `0.1`.
#'
#' @param jitter_width
#' Horizontal jitter width applied to individual
#' cell expression points.
#' Default is `0.12`.
#'
#' @param cor_label_size
#' Text size used for correlation coefficients
#' and composite ranking scores in trend plots.
#' Default is `3.5`.
#'
#' @param gene_title_size
#' Font size used for gene names in expression
#' trend plot facet labels.
#' Default is `14`.
#'
#' @param base_size
#' Base font size used in ggplot2 visualizations.
#' Default is `12`.
#'
#' @param save_dir
#' Directory used to save analysis tables and plots.
#' Default is `"./results"`.
#' Set to `NULL` to disable file saving.
#'
#' @param save_deg_csv
#' Filename for the complete differential expression table.
#' Default is `"DEG_result.csv"`.
#'
#' @param save_sig_deg_csv
#' Filename for significant differentially expressed genes.
#' Default is `"Significant_DEG.csv"`.
#'
#' @param save_cor_csv
#' Filename for viral-load correlation results.
#' Default is `"Correlation_result.csv"`.
#'
#' @param save_combined_csv
#' Filename for the integrated DEG and correlation results.
#' Default is `"DEG_Correlation_result.csv"`.
#'
#' @param save_sig_combined_csv
#' Filename for significant infection-associated genes
#' with composite ranking scores.
#' Default is `"Significant_DEG_Correlation.csv"`.
#'
#' @param save_volcano
#' Filename for the DEG volcano plot.
#' Default is `"Volcano_plot.pdf"`.
#'
#' @param save_correlation
#' Filename for the correlation–significance plot.
#' Default is `"DEG_Correlation_plot.pdf"`.
#'
#' @param save_trend
#' Base filename for gene expression trend plots.
#'
#' Two separate PDF files are generated by appending
#' `"_Positive.pdf"` and `"_Negative.pdf"` to the base name.
#'
#' Default is `"Gene_Expression_Trend.pdf"`.
#'
#' @param width
#' Width of DEG volcano and correlation–significance
#' plots in inches.
#' Default is `7`.
#'
#' @param height
#' Height of DEG volcano and correlation–significance
#' plots in inches.
#' Default is `6`.
#'
#' @param trend_width
#' Width of gene expression trend plots in inches.
#' Default is `12`.
#'
#' @param trend_height
#' Height of gene expression trend plots in inches.
#' Default is `8`.
#'
#' @param plot_title
#' Title displayed in the DEG volcano plot.
#' Default is `"Volcano Plot: Infected vs Uninfected"`.
#'
#' @return
#' An invisible named list containing the following components:
#'
#' \describe{
#'   \item{`n_case`}{
#'     Number of infected cells used for differential expression analysis.
#'   }
#'   \item{`n_control`}{
#'     Number of uninfected cells used for differential expression analysis.
#'   }
#'   \item{`n_correlation_cells`}{
#'     Number of cells included in correlation analysis.
#'   }
#'   \item{`cells_used`}{
#'     Cell barcodes included in the filtered analysis object.
#'   }
#'   \item{`correlation_cells`}{
#'     Cell barcodes included in viral-load correlation analysis.
#'   }
#'   \item{`deg_table`}{
#'     Complete differential expression results after gene exclusions.
#'   }
#'   \item{`deg_significant`}{
#'     Significantly upregulated and downregulated genes.
#'   }
#'   \item{`cor_table`}{
#'     Correlation coefficients, P values, adjusted P values,
#'     and correlation classifications.
#'   }
#'   \item{`combined_table`}{
#'     Integrated differential expression and viral-load
#'     correlation results, including composite ranking metrics
#'     for significant candidates.
#'   }
#'   \item{`significant_combined`}{
#'     Significant infection-associated genes ranked
#'     by their composite scores.
#'   }
#'   \item{`significant_genes`}{
#'     Character vector containing the significant
#'     infection-associated gene names.
#'   }
#'   \item{`top_positive`}{
#'     Top-ranked genes positively correlated with viral load.
#'   }
#'   \item{`top_negative`}{
#'     Top-ranked genes negatively correlated with viral load.
#'   }
#'   \item{`volcano`}{
#'     ggplot object containing the DEG volcano plot.
#'   }
#'   \item{`correlation_plot`}{
#'     ggplot object containing the correlation–significance plot.
#'   }
#'   \item{`trend_plot_positive`}{
#'     ggplot object showing expression trends of positively
#'     correlated genes, or `NULL` if no genes are selected.
#'   }
#'   \item{`trend_plot_negative`}{
#'     ggplot object showing expression trends of negatively
#'     correlated genes, or `NULL` if no genes are selected.
#'   }
#'   \item{`trend_data_positive`}{
#'     Cell-level expression data for positive correlation trends.
#'   }
#'   \item{`trend_data_negative`}{
#'     Cell-level expression data for negative correlation trends.
#'   }
#'   \item{`trend_summary_positive`}{
#'     Infection-level expression summaries for positive genes.
#'   }
#'   \item{`trend_summary_negative`}{
#'     Infection-level expression summaries for negative genes.
#'   }
#'   \item{`genes_used`}{
#'     Significant DEGs included in correlation testing.
#'   }
#'   \item{`object`}{
#'     Filtered Seurat object used for the analysis,
#'     including the `infection_status` metadata column.
#'   }
#' }
#'
#' @note
#' The integrated ranking score is intended to prioritize
#' candidate infection-associated host genes rather than
#' establish causal relationships with viral infection.
#'
#' Differential expression and correlation analyses are
#' performed at the individual-cell level.
#' For datasets containing multiple biological replicates,
#' potential donor-level or sample-level dependencies
#' should be considered when interpreting statistical significance.
#'
#' The function performs correlation analysis only for genes
#' passing the differential expression thresholds.
#' Consequently, genes associated with viral load but not
#' classified as significant DEGs are not included in the
#' integrated ranking.
#'
#' The `layer` parameter controls expression extraction for
#' correlation and trend analyses. Differential expression
#' is performed through `Seurat::FindMarkers()` using the
#' selected assay and its applicable expression data.
#'
#' @examples
#' \dontrun{
#'
#' # Basic infection-associated gene analysis
#' result <- InfectionVariableGenes(
#'   obj = seurat_obj,
#'   group_col = "timepoint",
#'   cell_type = "manual_label",
#'   level_col = "infection_level",
#'   viral_col = "viral_load",
#'   group_include = "24hpi",
#'   celltype_include = "CD4+Tcells",
#'   case_levels = 1:4,
#'   control_levels = 0
#' )
#'
#' # Examine significant infection-associated genes
#' head(result$significant_combined)
#'
#' # Extract top positively and negatively correlated genes
#' positive_genes <- result$top_positive
#' negative_genes <- result$top_negative
#'
#' # Display the DEG volcano plot
#' print(result$volcano)
#'
#' # Display the correlation-significance plot
#' print(result$correlation_plot)
#'
#' # Display positive correlation expression trends
#' print(result$trend_plot_positive)
#'
#' # Analyze multiple infection time points
#' result2 <- InfectionVariableGenes(
#'   obj = seurat_obj,
#'   group_include = c("24hpi", "48hpi", "72hpi"),
#'   celltype_include = "CD4+Tcells",
#'   cor_method = "spearman",
#'   deg_p_cut = 0.05,
#'   deg_log2fc_cut = 0.5,
#'   cor_p_cut = 0.05,
#'   cor_cut = 0.3,
#'   trend_top_n = 10,
#'   save_dir = "./SingleViral_results"
#' )
#'
#' # Analyze selected cell populations
#' result3 <- InfectionVariableGenes(
#'   obj = seurat_obj,
#'   group_include = "48hpi",
#'   celltype_include = c("CD4+Tcells", "CD8+Tcells"),
#'   virus_genes = viral_genes,
#'   mt_pattern = "^MT-",
#'   save_dir = "./results/T_cells"
#' )
#'
#' # Disable automatic file saving
#' result4 <- InfectionVariableGenes(
#'   obj = seurat_obj,
#'   save_dir = NULL
#' )
#'
#' }
#'
#' @seealso
#' [ViralLoadScoring()],
#' [TargetInfectionLevel()],
#' [GeneViralLoadCorrelation()],
#' [InfectionDynamic()]
#'
#' @importFrom Seurat Assays DefaultAssay Idents FindMarkers GetAssayData
#' @importFrom dplyr case_when bind_rows left_join percent_rank
#' @importFrom dplyr group_by summarise n
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon
#' @importFrom ggplot2 geom_violin geom_text geom_hline geom_vline
#' @importFrom ggplot2 geom_vline scale_color_manual scale_alpha_manual
#' @importFrom ggplot2 scale_color_gradient scale_x_continuous
#' @importFrom ggplot2 theme_bw theme_classic theme labs
#' @importFrom ggplot2 element_blank element_text
#' @importFrom ggplot2 position_jitter facet_wrap ggsave
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggnewscale new_scale_color
#' @importFrom readr write_csv
#' @importFrom tibble rownames_to_column
#' @importFrom stats cor.test p.adjust sd
#' @importFrom grid unit
#'
#' @export
InfectionVariableGenes <- function(
    obj,
    group_col = "timepoint", cell_type = "manual_label",
    level_col = "infection_level", viral_col = "viral_load",
    group_include = "24hpi", celltype_include = "CD4+Tcells",
    case_levels = 1:4, control_levels = 0,
    group_levels = c("control","24hpi","48hpi","72hpi","96hpi"),
    correlation_levels = 0:4, min_cells = 50,
    virus_genes = character(0), mt_genes = NULL, mt_pattern = "^MT-",
    test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0,
    only.pos = FALSE, pseudocount.use = 1,
    assay = "RNA", layer = "data",
    deg_p_cut = 0.1, deg_log2fc_cut = 0.5,
    cor_method = c("spearman","pearson"), cor_p_cut = 0.05, cor_cut = 0.3,
    label_genes = NULL, top_n_pval = 30,
    point_size = 1.5, alpha_nochange = 0.6,
    color_up = "#F45B69", color_down = "#4EA8DE",
    color_nochange = "gray70", line_color = "#999999",
    repel_max_overlaps = 100,
    cor_point_size = 2, cor_label_top_n = 20,
    trend_top_n = 10, trend_genes = NULL,
    summarize_fun = mean, show_se = TRUE,
    trend_point_size = 2.2, trend_line_width = 1.2,
    violin_alpha = 0.25, show_violin = TRUE, violin_width = 0.85,
    show_cell_points = TRUE, cell_point_alpha = 0.15,
    cell_point_size = 0.1, jitter_width = 0.12,
    cor_label_size = 3.5, gene_title_size = 14, base_size = 12,
    save_dir = "./results",
    save_deg_csv = "DEG_result.csv",
    save_sig_deg_csv = "Significant_DEG.csv",
    save_cor_csv = "Correlation_result.csv",
    save_combined_csv = "DEG_Correlation_result.csv",
    save_sig_combined_csv = "Significant_DEG_Correlation.csv",
    save_volcano = "Volcano_plot.pdf",
    save_correlation = "DEG_Correlation_plot.pdf",
    save_trend = "Gene_Expression_Trend.pdf",
    width = 7, height = 6, trend_width = 12, trend_height = 8,
    plot_title = "Volcano Plot: Infected vs Uninfected"
) {
  stopifnot(inherits(obj, "Seurat"))
  cor_method <- match.arg(cor_method)
  packages <- c("Seurat","dplyr","ggplot2","ggrepel","ggnewscale","readr","tibble","Matrix")
  missing_pkgs <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing_pkgs)) stop("Missing packages: ", paste(missing_pkgs, collapse=", "))
  if (!assay %in% Seurat::Assays(obj)) stop("Assay not found: ", assay)
  if (length(min_cells) != 1L || !is.finite(min_cells) || min_cells < 3)
    stop("min_cells must be >= 3.")
  if (length(trend_top_n) != 1L || !is.finite(trend_top_n) || trend_top_n < 0)
    stop("trend_top_n must be >= 0.")
  if (length(intersect(as.character(case_levels), as.character(control_levels))))
    stop("case_levels and control_levels overlap.")
  md <- obj@meta.data
  need <- c(group_col, level_col, viral_col, cell_type)
  missing_cols <- setdiff(need, colnames(md))
  if (length(missing_cols)) stop("Missing metadata columns: ", paste(missing_cols, collapse=", "))
  
  keep <- rep(TRUE, nrow(md))
  if (!is.null(group_include))
    keep <- keep & as.character(md[[group_col]]) %in% as.character(group_include)
  if (!is.null(cell_type) && !is.null(celltype_include))
    keep <- keep & as.character(md[[cell_type]]) %in% as.character(celltype_include)
  cells <- rownames(md)[keep]
  if (!length(cells)) stop("No cells after group/cell type filtering.")
  obj_sel <- subset(obj, cells=cells)
  lev <- as.character(obj_sel@meta.data[[level_col]])
  obj_sel$infection_status <- ifelse(
    lev %in% as.character(case_levels), "infected",
    ifelse(lev %in% as.character(control_levels), "uninfected", NA_character_)
  )
  cells <- colnames(obj_sel)[!is.na(obj_sel$infection_status)]
  if (!length(cells)) stop("No valid infection levels.")
  obj_sel <- subset(obj_sel, cells=cells)
  n_case <- sum(obj_sel$infection_status == "infected")
  n_control <- sum(obj_sel$infection_status == "uninfected")
  if (min(n_case,n_control) < 3) stop("At least 3 cells per DEG group are needed.")
  if (!is.null(group_levels)) {
    tp <- as.character(obj_sel@meta.data[[group_col]])
    obj_sel@meta.data[[group_col]] <- factor(tp, levels=unique(c(group_levels,tp)))
  }
  Seurat::DefaultAssay(obj_sel) <- assay
  Seurat::Idents(obj_sel) <- "infection_status"
  deg <- Seurat::FindMarkers(
    obj_sel, ident.1="infected", ident.2="uninfected",
    test.use=test.use, min.pct=min.pct,
    logfc.threshold=logfc.threshold, only.pos=only.pos,
    pseudocount.use=pseudocount.use, assay=assay
  )
  if (!nrow(deg)) stop("FindMarkers returned no genes.")
  deg <- tibble::rownames_to_column(deg, "gene")
  if (!"avg_log2FC" %in% names(deg)) {
    if ("avg_logFC" %in% names(deg)) deg$avg_log2FC <- deg$avg_logFC
    else stop("No avg_log2FC/avg_logFC column.")
  }
  if (!"p_val_adj" %in% names(deg)) {
    if (!"p_val" %in% names(deg)) stop("DEG results have no P values.")
    deg$p_val_adj <- stats::p.adjust(deg$p_val, "BH")
  }
  if (!all(c("pct.1","pct.2") %in% names(deg)))
    stop("DEG results must contain pct.1 and pct.2.")
  excluded <- unique(c(virus_genes,mt_genes))
  deg <- deg[!deg$gene %in% excluded, , drop=FALSE]
  if (!is.null(mt_pattern) && nzchar(mt_pattern))
    deg <- deg[!grepl(mt_pattern,deg$gene),,drop=FALSE]
  if (!nrow(deg)) stop("No genes after exclusions.")
  deg$neg_log10_p <- -log10(pmax(deg$p_val_adj,.Machine$double.xmin))
  deg$change <- dplyr::case_when(
    !is.na(deg$p_val_adj) & deg$p_val_adj < deg_p_cut &
      deg$avg_log2FC > deg_log2fc_cut ~ "Up",
    !is.na(deg$p_val_adj) & deg$p_val_adj < deg_p_cut &
      deg$avg_log2FC < -deg_log2fc_cut ~ "Down",
    TRUE ~ "No change"
  )
  deg_sig <- deg[deg$change != "No change",,drop=FALSE]
  message("DEG: ", nrow(deg_sig), "; infected cells: ",n_case,
          "; uninfected cells: ",n_control)
  
  lab_set <- if (!is.null(label_genes)) intersect(label_genes,deg$gene) else {
    head(deg_sig$gene[order(deg_sig$p_val_adj)],max(0,top_n_pval))
  }
  deg$label_flag <- deg$gene %in% lab_set
  p_vol <- ggplot2::ggplot(deg,ggplot2::aes(x=avg_log2FC,y=neg_log10_p)) +
    ggplot2::geom_hline(yintercept=-log10(deg_p_cut),linetype="dashed",color=line_color) +
    ggplot2::geom_vline(xintercept=c(-deg_log2fc_cut,deg_log2fc_cut),
                        linetype="dashed",color=line_color) +
    ggplot2::geom_point(ggplot2::aes(color=change,alpha=change),size=point_size) +
    ggplot2::scale_color_manual(values=c(Up=color_up,Down=color_down,
                                         "No change"=color_nochange)) +
    ggplot2::scale_alpha_manual(values=c(Up=1,Down=1,"No change"=alpha_nochange),
                                guide="none") +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(panel.grid=ggplot2::element_blank(),
                   plot.title=ggplot2::element_text(hjust=.5)) +
    ggplot2::labs(x="Log2FC",y="-Log10(FDR)",title=plot_title,color="DEG")
  label_df <- deg[deg$label_flag & deg$change != "No change",,drop=FALSE]
  if (nrow(label_df))
    p_vol <- p_vol + ggrepel::geom_text_repel(
      data=label_df,ggplot2::aes(label=gene,color=change),
      size=3,max.overlaps=repel_max_overlaps,show.legend=FALSE
    )
  
  md_cor <- obj_sel@meta.data
  v_all <- suppressWarnings(as.numeric(as.character(md_cor[[viral_col]])))
  keep_cor <- as.character(md_cor[[level_col]]) %in%
    as.character(correlation_levels) & is.finite(v_all)
  cells_cor <- rownames(md_cor)[keep_cor]
  viral_load <- v_all[keep_cor]
  if (length(cells_cor) < min_cells) stop("Too few cells for correlation.")
  if (!is.finite(stats::sd(viral_load)) || stats::sd(viral_load)==0)
    stop("Viral load has no variation.")
  expr <- Seurat::GetAssayData(obj_sel,assay=assay,layer=layer)
  genes_use <- intersect(deg_sig$gene,rownames(expr))
  cor_one <- function(gene) {
    x <- as.numeric(expr[gene,cells_cor,drop=TRUE])
    ok <- is.finite(x) & is.finite(viral_load)
    x <- x[ok]; y <- viral_load[ok]; n <- length(x)
    if (n < min_cells || length(unique(x))<2 || length(unique(y))<2)
      return(data.frame(gene=gene,n=n,r=NA_real_,p=NA_real_))
    ct <- tryCatch(suppressWarnings(stats::cor.test(
      x,y,method=cor_method,exact=FALSE
    )),error=function(e) NULL)
    if (is.null(ct)) return(data.frame(gene=gene,n=n,r=NA_real_,p=NA_real_))
    data.frame(gene=gene,n=n,r=unname(ct$estimate),p=ct$p.value)
  }
  cor_table <- if (length(genes_use))
    dplyr::bind_rows(lapply(genes_use,cor_one)) else
      data.frame(gene=character(),n=integer(),r=numeric(),p=numeric())
  cor_table$p_adj <- stats::p.adjust(cor_table$p,method="BH")
  cor_table$cor_change <- dplyr::case_when(
    !is.na(cor_table$p_adj) & cor_table$p_adj < cor_p_cut &
      cor_table$r > cor_cut ~ "Positive",
    !is.na(cor_table$p_adj) & cor_table$p_adj < cor_p_cut &
      cor_table$r < -cor_cut ~ "Negative",
    TRUE ~ "No correlation"
  )
  cor_table$label <- sprintf("%s r=%.2f, FDR=%.2g",
                             cor_method,cor_table$r,cor_table$p_adj)
  combined_table <- dplyr::left_join(deg_sig,cor_table,by="gene")
  combined_table$classification <- dplyr::case_when(
    combined_table$change=="Up" & combined_table$cor_change=="Positive" ~ "Up-Positive",
    combined_table$change=="Up" & combined_table$cor_change=="Negative" ~ "Up-Negative",
    combined_table$change=="Down" & combined_table$cor_change=="Positive" ~ "Down-Positive",
    combined_table$change=="Down" & combined_table$cor_change=="Negative" ~ "Down-Negative",
    TRUE ~ "Not significant"
  )
  # Equal-weight percentile ranking, computed across ALL significant correlated DEGs.
  # Up genes use pct.1, down genes use pct.2.
  sig <- combined_table[combined_table$classification!="Not significant" &
                          is.finite(combined_table$r) &
                          is.finite(combined_table$p_val_adj),,drop=FALSE]
  sig$pct_weight <- ifelse(sig$change=="Up",sig$pct.1,sig$pct.2)
  rank01 <- function(x) {
    if (!length(x)) return(numeric(0))
    if (length(x)==1L) return(1)
    dplyr::percent_rank(x)
  }
  sig$rho_rank <- rank01(abs(sig$r))
  sig$deg_rank <- rank01(-sig$p_val_adj)
  sig$pct_rank <- rank01(sig$pct_weight)
  sig$rank_score <- rowMeans(
    sig[,c("rho_rank","deg_rank","pct_rank"),drop=FALSE],na.rm=FALSE
  )
  sig <- sig[order(-sig$rank_score,sig$p_val_adj,-abs(sig$r)),,drop=FALSE]
  sig$rank_overall <- seq_len(nrow(sig))
  sig$rank_direction <- ave(seq_len(nrow(sig)),
                            ifelse(sig$r>0,"Positive","Negative"),
                            FUN=seq_along)
  combined_table <- dplyr::left_join(
    combined_table,
    sig[,c("gene","pct_weight","rho_rank","deg_rank","pct_rank",
           "rank_score","rank_overall","rank_direction"),drop=FALSE],
    by="gene"
  )
  significant_combined <- sig
  positive_top <- head(sig[sig$r>0,,drop=FALSE],trend_top_n)
  negative_top <- head(sig[sig$r<0,,drop=FALSE],trend_top_n)
  message("Correlated DEG: ",nrow(sig),"; positive top: ",nrow(positive_top),
          "; negative top: ",nrow(negative_top))
  
  # ============================================================
  # Correlation-FDR volcano: two independent colorbar legends
  # Positive r: red gradient based on pct.1
  # Negative r: blue gradient based on pct.2
  # ============================================================
  plot_df <- combined_table
  plot_df$neg_log10_p <- -log10(pmax(plot_df$p_val_adj, .Machine$double.xmin))
  plot_df$plot_direction <- dplyr::case_when(
    plot_df$classification != "Not significant" &
      is.finite(plot_df$r) & plot_df$r > 0 ~ "Positive",
    plot_df$classification != "Not significant" &
      is.finite(plot_df$r) & plot_df$r < 0 ~ "Negative",
    TRUE ~ "Not significant"
  )
  bg <- plot_df[plot_df$plot_direction == "Not significant", , drop=FALSE]
  pos <- plot_df[plot_df$plot_direction == "Positive", , drop=FALSE]
  neg <- plot_df[plot_df$plot_direction == "Negative", , drop=FALSE]
  
  p_cor <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept=-log10(deg_p_cut),
                        linetype="dashed", color=line_color) +
    ggplot2::geom_vline(xintercept=c(-cor_cut,cor_cut),
                        linetype="dashed", color=line_color) +
    ggplot2::geom_vline(xintercept=0, linetype="dotted", color="grey70") +
    ggplot2::geom_point(data=bg,
                        ggplot2::aes(x=r,y=neg_log10_p),
                        color=color_nochange, alpha=alpha_nochange,
                        size=cor_point_size, na.rm=TRUE) +
    ggplot2::geom_point(data=pos,
                        ggplot2::aes(x=r,y=neg_log10_p,color=pct.1),
                        size=cor_point_size, na.rm=TRUE) +
    ggplot2::scale_color_gradient(
      name="pct1", low="#FADBD8", high="#B71C1C",
      limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1),
      guide=ggplot2::guide_colorbar(
        title.position="top", barheight=grid::unit(4,"cm"),
        barwidth=grid::unit(0.6,"cm"),
        frame.colour=NA, ticks.colour="black"
      )
    ) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(data=neg,
                        ggplot2::aes(x=r,y=neg_log10_p,color=pct.2),
                        size=cor_point_size, na.rm=TRUE) +
    ggplot2::scale_color_gradient(
      name="pct2", low="#D6EAF8", high="#154360",
      limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1),
      guide=ggplot2::guide_colorbar(
        title.position="top", barheight=grid::unit(4,"cm"),
        barwidth=grid::unit(0.6,"cm"),
        frame.colour=NA, ticks.colour="black"
      )
    ) +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(
      panel.grid=ggplot2::element_blank(),
      plot.title=ggplot2::element_text(hjust=0.5,size=15),
      legend.position="right", legend.box="vertical",
      legend.title=ggplot2::element_text(size=12,face="bold"),
      legend.text=ggplot2::element_text(size=10),
      legend.spacing.y=grid::unit(0.5,"cm")
    ) +
    ggplot2::labs(
      x=paste0(cor_method," correlation with viral load (r)"),
      y="-Log10(DEG FDR)",
      title="DEG significance vs viral-load correlation"
    )
  
  # Label top-ranked genes, balanced between correlation directions.
  if (is.null(label_genes)) {
    k_pos <- ceiling(cor_label_top_n/2)
    k_neg <- floor(cor_label_top_n/2)
    cor_labels <- c(head(positive_top$gene,k_pos),
                    head(negative_top$gene,k_neg))
  } else {
    cor_labels <- intersect(label_genes,significant_combined$gene)
  }
  cor_label_df <- plot_df[plot_df$gene %in% cor_labels,,drop=FALSE]
  if (nrow(cor_label_df)) {
    p_cor <- p_cor + ggrepel::geom_text_repel(
      data=cor_label_df,
      ggplot2::aes(x=r,y=neg_log10_p,label=gene),
      size=3, color="black", max.overlaps=repel_max_overlaps,
      show.legend=FALSE, na.rm=TRUE
    )
  }
  
  # Trend plots are generated separately for positive and negative genes.
  make_trend <- function(selected, title) {
    if (!nrow(selected)) return(list(plot=NULL,long=NULL,summary=NULL))
    selected <- selected[selected$gene %in% rownames(expr),,drop=FALSE]
    if (!nrow(selected)) return(list(plot=NULL,long=NULL,summary=NULL))
    gs <- selected$gene
    mat <- as.matrix(expr[gs,cells_cor,drop=FALSE])
    dat <- data.frame(
      cell=rep(cells_cor,each=length(gs)),
      gene=rep(gs,times=length(cells_cor)),
      expr=as.vector(mat),
      level=rep(suppressWarnings(as.numeric(as.character(
        md_cor[cells_cor,level_col]))),each=length(gs)),
      viral_load=rep(viral_load,each=length(gs))
    )
    dat <- dat[is.finite(dat$level) & is.finite(dat$expr),,drop=FALSE]
    dat$gene <- factor(dat$gene,levels=gs)
    summ <- dat |>
      dplyr::group_by(gene,level) |>
      dplyr::summarise(
        mean_expr=summarize_fun(expr,na.rm=TRUE),
        se_expr=stats::sd(expr,na.rm=TRUE)/sqrt(sum(is.finite(expr))),
        n_cells=dplyr::n(),.groups="drop"
      )
    summ$se_expr[!is.finite(summ$se_expr)] <- 0
    ann <- selected[,c("gene","r","p_adj","rank_score"),drop=FALSE]
    ann$gene <- factor(ann$gene,levels=gs)
    ann$label <- sprintf("r=%.2f | Score=%.3f",ann$r,ann$rank_score)
    # Use facet strip for gene names; annotation at top of each facet.
    ann$x <- min(dat$level,na.rm=TRUE)
    ann$y <- Inf
    p <- ggplot2::ggplot()
    if (show_violin) {
      vv <- dat[dat$expr>0,,drop=FALSE]
      if (nrow(vv))
        p <- p + ggplot2::geom_violin(
          data=vv,ggplot2::aes(x=level,y=expr,group=interaction(gene,level)),
          fill="#D9902F",alpha=violin_alpha,color=NA,
          width=violin_width,trim=TRUE,scale="width"
        )
    }
    if (show_cell_points)
      p <- p + ggplot2::geom_point(
        data=dat,ggplot2::aes(x=level,y=expr),
        position=ggplot2::position_jitter(width=jitter_width,height=0),
        size=cell_point_size,alpha=cell_point_alpha,color="#D9902F"
      )
    if (show_se)
      p <- p + ggplot2::geom_ribbon(
        data=summ,ggplot2::aes(x=level,ymin=mean_expr-se_expr,
                               ymax=mean_expr+se_expr,group=gene),
        fill="#9E2F2F",alpha=.18
      )
    p <- p +
      ggplot2::geom_line(data=summ,ggplot2::aes(x=level,y=mean_expr,group=gene),
                         linewidth=trend_line_width,color="#9E2F2F") +
      ggplot2::geom_point(data=summ,ggplot2::aes(x=level,y=mean_expr),
                          size=trend_point_size,color="#9E2F2F") +
      ggplot2::geom_text(data=ann,ggplot2::aes(x=x,y=y,label=label),
                         hjust=0,vjust=1.15,size=cor_label_size,
                         fontface="bold",inherit.aes=FALSE) +
      ggplot2::facet_wrap(~gene,scales="free_y") +
      ggplot2::scale_x_continuous(breaks=sort(unique(dat$level))) +
      ggplot2::theme_classic(base_size=base_size) +
      ggplot2::theme(plot.title=ggplot2::element_text(hjust=.5,face="bold"),
                     axis.title=ggplot2::element_text(face="bold"),
                     strip.text=ggplot2::element_text(face="bold",size=gene_title_size)) +
      ggplot2::labs(title=title,x="Infection level",
                    y=paste0("Expression (",layer,")"))
    list(plot=p,long=dat,summary=summ)
  }
  if (!is.null(trend_genes)) {
    positive_top <- positive_top[positive_top$gene %in% trend_genes,,drop=FALSE]
    negative_top <- negative_top[negative_top$gene %in% trend_genes,,drop=FALSE]
  }
  pos_out <- make_trend(positive_top,"Positive viral-load correlation")
  neg_out <- make_trend(negative_top,"Negative viral-load correlation")
  
  if (!is.null(save_dir)) {
    dir.create(save_dir,recursive=TRUE,showWarnings=FALSE)
    save_csv <- function(x,name) {
      if (!is.null(name)) readr::write_csv(x,file.path(save_dir,name))
    }
    save_pdf <- function(p,name,w,h) {
      if (!is.null(p) && !is.null(name))
        ggplot2::ggsave(file.path(save_dir,name),plot=p,
                        width=w,height=h,useDingbats=FALSE)
    }
    save_csv(deg,save_deg_csv)
    save_csv(deg_sig,save_sig_deg_csv)
    save_csv(cor_table,save_cor_csv)
    save_csv(combined_table,save_combined_csv)
    save_csv(significant_combined,save_sig_combined_csv)
    save_csv(positive_top,"Top_positive_genes.csv")
    save_csv(negative_top,"Top_negative_genes.csv")
    save_pdf(p_vol,save_volcano,width,height)
    save_pdf(p_cor,save_correlation,width,height)
    # save_trend is treated as a base name for two separate PDF files.
    if (!is.null(save_trend)) {
      base <- sub("\\.[Pp][Dd][Ff]$","",save_trend)
      save_pdf(pos_out$plot,paste0(base,"_Positive.pdf"),trend_width,trend_height)
      save_pdf(neg_out$plot,paste0(base,"_Negative.pdf"),trend_width,trend_height)
    }
  }
  invisible(list(
    n_case=n_case,n_control=n_control,
    n_correlation_cells=length(cells_cor),
    cells_used=colnames(obj_sel),correlation_cells=cells_cor,
    deg_table=deg,deg_significant=deg_sig,cor_table=cor_table,
    combined_table=combined_table,significant_combined=significant_combined,
    significant_genes=significant_combined$gene,
    top_positive=positive_top,top_negative=negative_top,
    volcano=p_vol,correlation_plot=p_cor,
    trend_plot_positive=pos_out$plot,trend_plot_negative=neg_out$plot,
    trend_data_positive=pos_out$long,trend_data_negative=neg_out$long,
    trend_summary_positive=pos_out$summary,trend_summary_negative=neg_out$summary,
    genes_used=genes_use,object=obj_sel
  ))
}
