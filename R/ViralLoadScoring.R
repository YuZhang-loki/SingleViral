#' Quantify Viral Transcriptional Burden at Single-Cell Resolution
#'
#' @description
#' `ViralLoadScoring()` quantifies viral transcriptional burden in individual
#' cells by integrating viral RNA abundance and viral transcriptional breadth
#' into a unified score.
#'
#' The function operates on Seurat objects and supports viral transcript
#' quantification, optional background correction, reference-based intensity
#' normalization, integrated viral-load scoring, and UMAP visualization.
#'
#' Unlike approaches relying solely on viral transcript counts, this function
#' considers both the relative abundance of viral transcripts and the
#' proportion of detected viral genes to characterize viral transcriptional
#' heterogeneity at single-cell resolution.
#'
#' @details
#' The calculation consists of the following steps:
#'
#' **1. Viral transcript quantification**
#'
#' For each cell, the total viral UMI count is calculated by summing the raw
#' counts of all detected viral genes. Viral transcript abundance is then
#' normalized by the total UMI count, including both viral and host transcripts.
#'
#' **2. Background correction**
#'
#' When `background_method = "negative_reference"`, the background viral
#' fraction is estimated as the median viral fraction across the cells
#' specified by `neg_cells`. This value is subtracted from the viral fraction
#' of each cell, with negative results truncated to zero.
#'
#' When `background_method = "none"`, no background correction is applied.
#'
#' **3. Viral expression intensity**
#'
#' Viral expression intensity is calculated as:
#'
#' `I = log1p(1e6 * F_adj)`
#'
#' where `F_adj` represents the background-adjusted viral fraction.
#'
#' **4. Viral transcriptional breadth**
#'
#' Viral transcriptional breadth is defined as the proportion of viral genes
#' detected in each cell:
#'
#' `B = D / G`
#'
#' where `D` is the number of detected viral genes and `G` is the total
#' number of viral genes included in the analysis.
#'
#' **5. Reference intensity normalization**
#'
#' Viral expression intensity is normalized using a reference threshold
#' determined by the specified quantile of positive viral intensities,
#' excluding negative reference cells when provided.
#'
#' The normalized intensity is calculated as:
#'
#' `I_scaled = min(I / I_ref, 1)`
#'
#' where `I_ref` is the reference intensity threshold.
#'
#' **6. Integrated viral-load scoring**
#'
#' The integrated viral transcriptional burden score is calculated as:
#'
#' `V = I_scaled^alpha * B^(1 - alpha)`
#'
#' where `alpha` controls the relative contributions of viral expression
#' intensity and transcriptional breadth.
#'
#' With the default `alpha = 0.5`, the score corresponds to the geometric
#' mean of the two components.
#'
#' The resulting score ranges from 0 to 1, with higher values indicating
#' greater relative viral RNA abundance and broader viral gene detection.
#'
#' Cells without detectable viral UMIs and cells explicitly specified
#' as negative references are assigned a final viral-load score of zero.
#'
#' **7. Reference model**
#'
#' The function can generate a reference model containing the viral gene set,
#' background correction parameters, intensity normalization threshold,
#' missing-gene handling strategy, and integration coefficient.
#'
#' When `reference_model` is supplied, the stored scoring parameters are
#' reused rather than recalculated, facilitating consistent scoring
#' across datasets.
#'
#' The reference model does not perform batch-effect correction.
#'
#' **8. Missing viral genes**
#'
#' Three strategies are supported:
#'
#' - `"zero"`: Missing viral genes are treated as having zero expression
#'   and retained in the transcriptional breadth denominator.
#' - `"filter"`: Missing viral genes are excluded from the analysis,
#'   and the breadth denominator is adjusted accordingly.
#' - `"error"`: The function stops if any predefined viral genes
#'   are absent from the expression matrix.
#'
#' @param obj
#' A Seurat object containing single-cell RNA sequencing data.
#' The selected assay must contain a `"counts"` layer with raw UMI counts.
#'
#' @param viral_genes
#' A character vector specifying the viral genes used for viral-load
#' quantification. Gene identifiers must match the row names of the
#' expression matrix. When `reference_model` is provided, the gene set
#' stored in the model is used instead.
#'
#' @param neg_cells
#' Optional character vector specifying negative reference cell barcodes.
#' These cells are used to estimate background viral signals when
#' `background_method = "negative_reference"`.
#' Cells specified in `neg_cells` are assigned a final viral-load score
#' of zero, regardless of their observed viral transcript counts.
#' Default is `NULL`.
#'
#' @param background_method
#' Character string specifying the background correction strategy.
#' Available options are:
#' \describe{
#'   \item{`"none"`}{No background correction (default).}
#'   \item{`"negative_reference"`}{
#'     Estimate the background viral fraction from negative reference cells.
#'   }
#' }
#'
#' @param assay
#' Character string specifying the Seurat assay used for quantification.
#' If `NULL`, the default assay of `obj` is used.
#' Default is `NULL`.
#'
#' @param alpha
#' Numeric value between 0 and 1 controlling the relative contributions
#' of viral expression intensity and transcriptional breadth to the
#' integrated score.
#' Larger values place greater emphasis on viral expression intensity.
#' Default is `0.5`.
#'
#' @param intensity_quantile
#' Numeric value between 0 and 1 specifying the quantile used to
#' determine the reference viral expression intensity.
#' Default is `0.95`.
#'
#' @param missing_gene_method
#' Character string specifying how viral genes absent from the expression
#' matrix are handled.
#' Available options are:
#' \describe{
#'   \item{`"zero"`}{
#'     Retain missing viral genes in the breadth denominator
#'     and assign them zero expression (default).
#'   }
#'   \item{`"filter"`}{
#'     Exclude missing viral genes and adjust the breadth denominator.
#'   }
#'   \item{`"error"`}{
#'     Stop the calculation if any predefined viral genes are missing.
#'   }
#' }
#'
#' @param reference_model
#' Optional reference model generated by a previous call to
#' `ViralLoadScoring()`.
#'
#' When provided, the stored viral gene set, integration coefficient,
#' intensity quantile, intensity reference threshold, background
#' correction settings, and missing-gene handling strategy are reused.
#'
#' This allows viral-load scores to be calculated using consistent
#' reference parameters across different datasets.
#' Default is `NULL`.
#'
#' @param return_model
#' Logical value indicating whether the reference model should be
#' included in the returned results.
#' Default is `TRUE`.
#'
#' @param plot
#' Logical value indicating whether UMAP-based visualizations should
#' be generated.
#'
#' When `TRUE`, the function generates feature plots showing the
#' integrated viral-load score, normalized viral expression intensity,
#' and viral transcriptional breadth.
#'
#' The Seurat object must contain a UMAP dimensional reduction.
#' Default is `TRUE`.
#'
#' @param group.col
#' Optional character string specifying a metadata column used to
#' split the viral-load UMAP visualization across experimental conditions,
#' sample groups, or cell populations.
#' Default is `NULL`.
#'
#' @return
#' A named list containing the following components:
#'
#' \describe{
#'   \item{`obj`}{
#'     Updated Seurat object containing viral-load metrics in its metadata.
#'   }
#'   \item{`assay`}{
#'     Assay used for viral transcript quantification.
#'   }
#'   \item{`genes_used`}{
#'     Viral genes included in the analysis.
#'   }
#'   \item{`missing_genes`}{
#'     Viral genes absent from the expression matrix.
#'   }
#'   \item{`neg_cells_used`}{
#'     Negative reference cells matched to the Seurat object.
#'   }
#'   \item{`group.col`}{
#'     Metadata column used for grouped visualization.
#'   }
#'   \item{`background_method`}{
#'     Background correction strategy used.
#'   }
#'   \item{`background_fraction`}{
#'     Estimated background viral fraction.
#'   }
#'   \item{`viral_load`}{
#'     Integrated viral transcriptional burden score for each cell.
#'   }
#'   \item{`viral_intensity`}{
#'     Log-transformed, background-adjusted viral expression intensity.
#'   }
#'   \item{`viral_intensity_scaled`}{
#'     Reference-normalized viral expression intensity.
#'   }
#'   \item{`viral_breadth`}{
#'     Proportion of detected viral genes in each cell.
#'   }
#'   \item{`viral_UMI`}{
#'     Total viral UMI count per cell.
#'   }
#'   \item{`viral_fraction`}{
#'     Fraction of total UMIs assigned to viral genes.
#'   }
#'   \item{`viral_fraction_adjusted`}{
#'     Background-adjusted viral fraction.
#'   }
#'   \item{`detected_genes`}{
#'     Number of detected viral genes per cell.
#'   }
#'   \item{`plots`}{
#'     List of ggplot objects generated when `plot = TRUE`,
#'     or `NULL` otherwise. When generated, the list contains
#'     `global`, `split`, `viral_load`, `intensity`, and `breadth`.
#'     The `split` component is `NULL` when `group.col = NULL`.
#'   }
#'   \item{`reference_model`}{
#'     Reference model containing scoring parameters.
#'     Included only when `return_model = TRUE`.
#'   }
#' }
#'
#' The updated Seurat object contains the following metadata columns:
#'
#' - `viral_load`: Integrated viral transcriptional burden score.
#' - `viral_load_intensity`: Log-transformed viral expression intensity.
#' - `viral_load_intensity_scaled`: Normalized viral expression intensity.
#' - `viral_load_breadth`: Viral transcriptional breadth.
#' - `viral_load_UMI`: Total viral UMI count.
#' - `viral_load_fraction`: Viral UMI fraction.
#' - `viral_load_fraction_adjusted`: Background-adjusted viral fraction.
#' - `viral_load_detected_genes`: Number of detected viral genes.
#'
#' @note
#' The viral-load score represents relative viral transcriptional burden
#' rather than a direct measurement of infectious virus production
#' or productive viral replication.
#'
#' The function requires raw UMI counts and is intended for single-cell
#' transcriptomic datasets with identifiable viral gene annotations.
#'
#' When a reference model is used, the viral gene set and scoring
#' parameters should be compatible with the target dataset.
#'
#' @examples
#' \dontrun{
#'
#' # Define viral genes
#' viral_genes <- c("IE1", "U79", "gB")
#'
#' # Basic viral-load quantification
#' result <- ViralLoadScoring(
#'   obj = seurat_obj,
#'   viral_genes = viral_genes,
#'   assay = "RNA"
#' )
#'
#' # Extract the updated Seurat object
#' seurat_obj <- result$obj
#'
#' # Examine viral-load scores
#' head(seurat_obj$viral_load)
#'
#' # Viral-load scoring with negative reference cells
#' result <- ViralLoadScoring(
#'   obj = seurat_obj,
#'   viral_genes = viral_genes,
#'   neg_cells = control_cells,
#'   background_method = "negative_reference",
#'   alpha = 0.5,
#'   intensity_quantile = 0.95
#' )
#'
#' # Save the reference model
#' model <- result$reference_model
#'
#' # Apply the reference model to another dataset
#' result2 <- ViralLoadScoring(
#'   obj = seurat_obj2,
#'   viral_genes = viral_genes,
#'   reference_model = model
#' )
#'
#' # Compare viral transcriptional burden across conditions
#' result3 <- ViralLoadScoring(
#'   obj = seurat_obj,
#'   viral_genes = viral_genes,
#'   group.col = "timepoint",
#'   plot = TRUE
#' )
#'
#' # Display the global viral-load UMAP
#' print(result3$plots$global)
#'
#' # Display viral-load scores by time point
#' print(result3$plots$split)
#'
#' }
#'
#' @seealso
#' [TargetInfectionLevel()],
#' [InfectionDynamic()],
#' [InfectionVariableGenes()],
#' [CellVariableInfection()]
#'
#' @importFrom Seurat DefaultAssay GetAssayData FeaturePlot
#' @importFrom Matrix colSums
#' @importFrom ggplot2 ggtitle theme element_text
#' @importFrom stats median quantile
#' @importFrom methods as
#'
#' @export
ViralLoadScoring <- function(
    obj,
    viral_genes,
    
    # ==========================================================
    # Negative reference cells
    # ==========================================================
    neg_cells = NULL,
    background_method = c(
      "none",
      "negative_reference"
    ),
    
    # ==========================================================
    # Assay
    # ==========================================================
    assay = NULL,
    
    # ==========================================================
    # Integrated viral score
    # ==========================================================
    alpha = 0.5,
    intensity_quantile = 0.95,
    
    # ==========================================================
    # Missing viral genes
    # ==========================================================
    missing_gene_method = c(
      "zero",
      "filter",
      "error"
    ),
    
    # ==========================================================
    # Reference model
    # ==========================================================
    reference_model = NULL,
    return_model = TRUE,
    
    # ==========================================================
    # Visualization
    # ==========================================================
    plot = TRUE,
    group.col = NULL
) {
  
  ############################################################
  ## 1. Package checks
  ############################################################
  
  if (!inherits(obj, "Seurat")) {
    stop("obj must be a Seurat object.")
  }
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.")
  }
  
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  
  ############################################################
  ## 2. Parameter validation
  ############################################################
  
  background_method <- match.arg(
    background_method
  )
  
  missing_gene_method <- match.arg(
    missing_gene_method
  )
  
  if (!is.numeric(alpha) ||
      length(alpha) != 1 ||
      !is.finite(alpha) ||
      alpha < 0 ||
      alpha > 1) {
    
    stop("alpha must be between 0 and 1.")
    
  }
  
  if (!is.numeric(intensity_quantile) ||
      length(intensity_quantile) != 1 ||
      !is.finite(intensity_quantile) ||
      intensity_quantile <= 0 ||
      intensity_quantile >= 1) {
    
    stop("intensity_quantile must be between 0 and 1.")
    
  }
  
  ############################################################
  ## 3. Select assay
  ############################################################
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(obj)
  }
  
  ############################################################
  ## 4. Extract counts
  ############################################################
  
  counts <- Seurat::GetAssayData(
    obj,
    assay = assay,
    layer = "counts"
  )
  
  if (!inherits(counts, "dgCMatrix")) {
    
    counts <- methods::as(
      counts,
      "dgCMatrix"
    )
    
  }
  
  if (anyNA(counts) ||
      any(!is.finite(counts@x)) ||
      any(counts@x < 0)) {
    
    stop("Invalid count matrix.")
    
  }
  
  ############################################################
  ## 5. Viral gene selection
  ############################################################
  
  if (is.null(reference_model)) {
    
    genes <- unique(
      as.character(viral_genes)
    )
    
    genes <- genes[
      !is.na(genes) & nzchar(genes)
    ]
    
  } else {
    
    required_fields <- c(
      "genes",
      "alpha",
      "intensity_reference"
    )
    
    if (!is.list(reference_model) ||
        !all(required_fields %in% names(reference_model))) {
      
      stop("Invalid reference_model.")
      
    }
    
    genes <- reference_model$genes
    
    alpha <- reference_model$alpha
    
    intensity_quantile <-
      reference_model$intensity_quantile
    
    if (!is.null(
      reference_model$missing_gene_method
    )) {
      
      missing_gene_method <-
        reference_model$missing_gene_method
      
    }
    
    if (!is.null(
      reference_model$background_method
    )) {
      
      background_method <-
        reference_model$background_method
      
    }
    
  }
  
  if (length(genes) == 0) {
    stop("No viral genes provided.")
  }
  
  ############################################################
  ## 6. Check missing viral genes
  ############################################################
  
  present_genes <- intersect(
    genes,
    rownames(counts)
  )
  
  missing_genes <- setdiff(
    genes,
    rownames(counts)
  )
  
  if (length(missing_genes) > 0) {
    
    message(
      "Total viral genes: ",
      length(genes)
    )
    
    message(
      "Present viral genes: ",
      length(present_genes)
    )
    
    message(
      "Missing viral genes: ",
      length(missing_genes)
    )
    
    if (missing_gene_method == "error") {
      
      stop(
        "Missing viral genes: ",
        paste(
          missing_genes,
          collapse = ", "
        )
      )
      
    }
    
    if (missing_gene_method == "filter") {
      
      message(
        "Missing genes will be excluded."
      )
      
      genes <- present_genes
      
    }
    
    if (missing_gene_method == "zero") {
      
      warning(
        length(missing_genes),
        " viral genes are absent from the expression matrix. ",
        "They will be treated as zero expression ",
        "and retained in the breadth denominator."
      )
      
    }
    
  }
  
  if (length(genes) == 0) {
    stop("No viral genes remain after filtering.")
  }
  
  ############################################################
  ## 7. Extract viral counts
  ############################################################
  
  present_genes <- intersect(
    genes,
    rownames(counts)
  )
  
  viral_counts <- counts[
    present_genes,
    ,
    drop = FALSE
  ]
  
  G <- length(genes)
  
  ############################################################
  ## 8. Identify negative cells
  ############################################################
  
  neg_cells_used <- character(0)
  
  if (!is.null(neg_cells)) {
    
    neg_cells_used <- intersect(
      colnames(counts),
      as.character(neg_cells)
    )
    
    if (length(neg_cells_used) == 0) {
      
      stop(
        "None of 'neg_cells' found among object cells."
      )
      
    }
    
    message(
      "Negative reference cells: ",
      length(neg_cells_used)
    )
    
  }
  
  ############################################################
  ## 9. Basic viral metrics
  ############################################################
  
  viral_UMI <- Matrix::colSums(
    viral_counts
  )
  
  total_UMI <- Matrix::colSums(
    counts
  )
  
  viral_fraction <- viral_UMI /
    pmax(total_UMI, 1)
  
  detected_genes <- Matrix::colSums(
    viral_counts > 0
  )
  
  ############################################################
  ## 10. Background correction
  ############################################################
  
  background_fraction <- 0
  
  if (is.null(reference_model)) {
    
    if (background_method == "negative_reference") {
      
      if (length(neg_cells_used) == 0) {
        
        stop(
          "neg_cells must be provided when ",
          "background_method = 'negative_reference'."
        )
        
      }
      
      background_fraction <- stats::median(
        viral_fraction[neg_cells_used],
        na.rm = TRUE
      )
      
    }
    
  } else {
    
    if (!is.null(
      reference_model$background_fraction
    )) {
      
      background_fraction <-
        reference_model$background_fraction
      
    } else if (
      background_method == "negative_reference"
    ) {
      
      stop(
        "reference_model lacks background_fraction."
      )
      
    }
    
  }
  
  viral_fraction_adjusted <- pmax(
    viral_fraction - background_fraction,
    0
  )
  
  ############################################################
  ## 11. Viral expression intensity
  ############################################################
  
  viral_intensity <- log1p(
    viral_fraction_adjusted * 1e6
  )
  
  ############################################################
  ## 12. Viral transcriptional breadth
  ############################################################
  
  viral_breadth <- detected_genes / G
  
  ############################################################
  ## 13. Reference scaling
  ############################################################
  
  if (is.null(reference_model)) {
    
    positive_intensity <- viral_intensity[
      viral_intensity > 0
    ]
    
    if (length(neg_cells_used) > 0) {
      
      positive_intensity <- viral_intensity[
        viral_intensity > 0 &
          !names(viral_intensity) %in% neg_cells_used
      ]
      
    }
    
    if (length(positive_intensity) == 0) {
      
      stop(
        "No positive viral intensity available ",
        "for reference scaling."
      )
      
    }
    
    intensity_reference <- as.numeric(
      stats::quantile(
        positive_intensity,
        probs = intensity_quantile,
        names = FALSE
      )
    )
    
  } else {
    
    intensity_reference <-
      reference_model$intensity_reference
    
  }
  
  if (!is.numeric(intensity_reference) ||
      length(intensity_reference) != 1 ||
      !is.finite(intensity_reference) ||
      intensity_reference <= 0) {
    
    stop("Invalid intensity reference.")
    
  }
  
  viral_intensity_scaled <- pmin(
    viral_intensity / intensity_reference,
    1
  )
  
  ############################################################
  ## 14. Integrated viral score
  ############################################################
  
  viral_load <- (
    viral_intensity_scaled ^ alpha
  ) * (
    viral_breadth ^ (1 - alpha)
  )
  
  ############################################################
  ## 15. Reset negative cells
  ############################################################
  
  if (length(neg_cells_used) > 0) {
    
    viral_load[neg_cells_used] <- 0
    
    viral_intensity[neg_cells_used] <- 0
    
    viral_intensity_scaled[neg_cells_used] <- 0
    
  }
  
  viral_load[viral_UMI == 0] <- 0
  
  ############################################################
  ## 16. Store metadata
  ############################################################
  
  cell_names <- colnames(obj)
  
  obj$viral_load <-
    viral_load[cell_names]
  
  obj$viral_load_intensity <-
    viral_intensity[cell_names]
  
  obj$viral_load_intensity_scaled <-
    viral_intensity_scaled[cell_names]
  
  obj$viral_load_breadth <-
    viral_breadth[cell_names]
  
  obj$viral_load_UMI <-
    viral_UMI[cell_names]
  
  obj$viral_load_fraction <-
    viral_fraction[cell_names]
  
  obj$viral_load_fraction_adjusted <-
    viral_fraction_adjusted[cell_names]
  
  obj$viral_load_detected_genes <-
    detected_genes[cell_names]
  
  ############################################################
  ## 17. Reference model
  ############################################################
  
  model <- list(
    
    version = "VEBIS_1.2",
    
    genes = genes,
    
    alpha = alpha,
    
    intensity_quantile =
      intensity_quantile,
    
    intensity_reference =
      intensity_reference,
    
    missing_gene_method =
      missing_gene_method,
    
    background_method =
      background_method,
    
    background_fraction =
      background_fraction,
    
    normalization = "total_CPM"
    
  )
  
  ############################################################
  ## 18. Visualization
  ############################################################
  
  plots <- NULL
  
  if (isTRUE(plot)) {
    
    p1 <- Seurat::FeaturePlot(
      obj,
      features = "viral_load",
      cols = c("grey90", "#E64B35"),
      order = TRUE
    ) +
      ggplot2::ggtitle(
        "Viral load scoring"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 18,
          face = "bold",
          hjust = 0.5
        )
      )
    
    p2 <- NULL
    
    if (!is.null(group.col)) {
      
      if (!group.col %in%
          colnames(obj@meta.data)) {
        
        stop(
          "group.col '",
          group.col,
          "' not found in obj@meta.data."
        )
        
      }
      
      p2 <- Seurat::FeaturePlot(
        obj,
        features = "viral_load",
        cols = c("grey90", "#E64B35"),
        split.by = group.col,
        order = TRUE
      )
      
    }
    
    plots <- list(
      global = p1,
      split = p2,
      viral_load = p1,
      intensity = Seurat::FeaturePlot(
        obj,
        features = "viral_load_intensity_scaled",
        cols = c("grey90", "#E64B35"),
        order = TRUE
      ),
      breadth = Seurat::FeaturePlot(
        obj,
        features = "viral_load_breadth",
        cols = c("grey90", "#E64B35"),
        order = TRUE
      )
    )
    
  }
  
  ############################################################
  ## 19. Return results
  ############################################################
  
  result <- list(
    
    obj = obj,
    
    assay = assay,
    
    genes_used = genes,
    
    missing_genes = missing_genes,
    
    neg_cells_used = neg_cells_used,
    
    group.col = group.col,
    
    background_method = background_method,
    
    background_fraction = background_fraction,
    
    viral_load = viral_load,
    
    viral_intensity = viral_intensity,
    
    viral_intensity_scaled =
      viral_intensity_scaled,
    
    viral_breadth = viral_breadth,
    
    viral_UMI = viral_UMI,
    
    viral_fraction = viral_fraction,
    
    viral_fraction_adjusted =
      viral_fraction_adjusted,
    
    detected_genes = detected_genes,
    
    plots = plots
    
  )
  
  if (isTRUE(return_model)) {
    
    result$reference_model <- model
    
  }
  
  return(result)
  
}
