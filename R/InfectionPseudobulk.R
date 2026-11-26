#' Patient-paired pseudobulk analysis across viral infection levels
#'
#' @description
#' Aggregates raw single-cell UMI counts by patient and infection level and
#' performs patient-paired host differential expression analysis with DESeq2.
#' By default, cells are selected using existing infection-level annotations.
#' Optionally, cells can be restricted to a common viral-load interval across
#' eligible patients and infected cells can be reclassified into global
#' quantile-based infection levels. The function supports Hallmark GSEA,
#' visualizations, and optional result export.
#'
#' @details
#' Each pseudobulk sample represents one patient-by-infection-level combination.
#' For each target level, only patients with QC-passing pseudobulk samples in
#' both the reference and target levels are included. Differential expression
#' is fitted using the paired design `~ patient_id + infection_level`.
#' This design accounts for patient-specific baseline differences and avoids
#' treating cells from the same patient as independent biological replicates.
#' It does not, by itself, control for differences in cell-type composition;
#' use `celltype_include` for cell-type-specific analyses when appropriate.
#'
#' Genes listed in `viral_genes` are excluded from host expression analysis.
#' If no viral genes are supplied, no automatic viral-gene removal occurs.
#' In overlap mode, only infected cells within the common viral-load interval
#' are retained, together with reference cells from eligible patients.
#' Infection levels are then recalculated using pooled viral-load quantiles.
#'
#' @param obj A Seurat object containing raw UMI counts and cell metadata.
#' @param patient_col Metadata column identifying patients or donors.
#' @param level_col Metadata column with existing infection-level labels.
#'   Used when `use_viral_overlap = FALSE`.
#' @param viral_col Metadata column containing numeric viral-load scores.
#'   Required for common-range overlap analysis.
#' @param cell_type Metadata column containing cell-type annotations.
#' @param celltype_include Optional character vector of cell types to retain;
#'   `NULL` includes all cell types.
#' @param ref_level Reference infection level, typically `0`.
#' @param case_levels Infection levels compared individually against `ref_level`.
#' @param n_levels Number of positive infection levels generated in overlap mode.
#'   In that mode, `case_levels` must equal `1:n_levels`.
#' @param zero_tol Nonnegative tolerance used to identify positive viral load
#'   in overlap mode.
#' @param use_viral_overlap Logical; if `FALSE` (default), use the supplied
#'   infection-level labels without overlap filtering or reclassification.
#'   If `TRUE`, apply common viral-load interval filtering and reclassification.
#' @param overlap_method Method for defining each patient's viral-load range:
#'   `"quantile"` or `"minmax"`. Used only in overlap mode.
#' @param overlap_quantiles Two increasing quantile probabilities defining
#'   patient-specific bounds when `overlap_method = "quantile"`.
#' @param min_infected_cells_per_patient Minimum infected cells required for a
#'   patient to be eligible for overlap analysis.
#' @param assay Seurat assay containing expression counts.
#' @param layer Layer containing nonnegative integer raw UMI counts, normally
#'   `"counts"`. Split layers must be joined before analysis.
#' @param min_cells_per_sample Minimum cells per patient-level pseudobulk sample.
#' @param min_patients Minimum number of paired patients required per contrast.
#' @param min_count Minimum count per gene in a sample for gene filtering.
#' @param min_samples Minimum number of samples meeting `min_count` for a gene
#'   to be retained.
#' @param viral_genes Character vector of viral gene names to exclude from
#'   host expression analyses. Names must match count-matrix row names.
#' @param fit_type Dispersion fitting method passed to `DESeq2::DESeq()`.
#' @param sf_type Size-factor estimation method for DESeq2. When `"ratio"`
#'   cannot be used because no gene is positive in every sample, the function
#'   falls back to `"poscounts"`.
#' @param alpha Significance level passed to `DESeq2::results()`.
#' @param padj_cut Adjusted p-value threshold for classifying DE genes.
#' @param logfc_cut Absolute log2 fold-change threshold for classifying DE genes.
#' @param run_joint_model Logical; additionally fit a joint model using patients
#'   with QC-passing samples at all requested levels.
#' @param run_gsea Logical; perform Hallmark gene set enrichment analysis.
#' @param species Species used for bundled Hallmark gene sets: `"human"`,
#'   `"mouse"`, or `"other"`. For `"other"`, provide `gmt_path`.
#' @param gmt_path Optional path to a GMT gene-set file. If `NULL`, bundled
#'   Hallmark sets or an available `msigdbr` fallback are used.
#' @param ranks_from GSEA ranking statistic: `"avg_log2FC"` (DESeq2 log2 fold
#'   change) or `"pval_adj_signed"` (signed -log10 adjusted p-value).
#' @param top_up Maximum positively enriched pathways displayed in GSEA plots.
#' @param top_down Maximum negatively enriched pathways displayed in GSEA plots.
#' @param drop_prefix Logical; remove `HALLMARK_` from displayed pathway names.
#' @param title_prefix Optional prefix for GSEA plot titles.
#' @param plot_volcano Logical; generate volcano plots for completed contrasts.
#' @param label_volcano Logical; label selected significant genes on volcano plots.
#' @param volcano_top_up Maximum upregulated genes labeled per volcano plot.
#' @param volcano_top_down Maximum downregulated genes labeled per volcano plot.
#' @param volcano_label_size Text size of volcano-plot gene labels.
#' @param plot_pca Logical; generate a pseudobulk PCA plot when feasible.
#' @param plot_gsea Logical; generate GSEA bar plots when GSEA is enabled.
#' @param plot_overlap Logical; generate a common viral-load interval plot in
#'   overlap mode. Has no effect when `use_viral_overlap = FALSE`.
#' @param plot_level_distribution Logical; plot infection-level composition
#'   across patients.
#' @param col_up Color used for upregulated genes and positively enriched sets.
#' @param col_down Color used for downregulated genes and negatively enriched sets.
#' @param save_results Logical; save results and generated plots to disk.
#' @param outdir Directory used when `save_results = TRUE`.
#'
#' @return A named list containing:
#' \describe{
#'   \item{analysis_mode, use_viral_overlap}{Selected analysis mode and its flag.}
#'   \item{patient_ranges, overlap_interval, retention_summary}{Overlap-mode
#'     diagnostics; may be `NULL` when overlap analysis is disabled.}
#'   \item{level_boundaries, level_distribution, cell_metadata}{Level
#'     boundaries (when reclassified), patient-level distribution, and retained
#'     cell metadata.}
#'   \item{pseudobulk_counts}{Aggregated host-gene count matrix.}
#'   \item{sample_metadata, sample_metadata_all}{Pseudobulk sample metadata
#'     after and before sample-level QC filtering.}
#'   \item{excluded_viral_genes}{Viral genes matched and excluded from analysis.}
#'   \item{DE, DDS}{Named lists of DESeq2 result tables and fitted objects
#'     for each infection-level comparison.}
#'   \item{GSEA}{Named list of GSEA results.}
#'   \item{paired_samples, comparison_summary}{Matched sample metadata and
#'     comparison status summaries.}
#'   \item{volcano_plots, gsea_plots, pca_plot, overlap_plot,
#'     level_distribution_plot}{Generated ggplot objects, or `NULL` when
#'     disabled or unavailable.}
#'   \item{joint_dds, joint_results}{Optional joint-model fit and contrasts.}
#' }
#'
#' @examples
#' \dontrun{
#' res <- InfectionPseudobulk(
#'   obj = seu,
#'   patient_col = "patient_id",
#'   level_col = "infection_level",
#'   ref_level = 0,
#'   case_levels = 1:2,
#'   viral_genes = c("gag", "pol", "env"),
#'   run_gsea = FALSE
#' )
#' res$comparison_summary
#' res$volcano_plots[["1"]]
#' }
#'
#' @export

InfectionPseudobulk <- function(
    obj,
    
    # ============================================================
    # Metadata
    # ============================================================
    patient_col = "patient_id",
    level_col = "infection_level",
    viral_col = "viral_load",
    cell_type = "manual_label",
    celltype_include = NULL,
    
    # ============================================================
    # Infection levels
    # ============================================================
    ref_level = 0,
    case_levels = 1:4,
    n_levels = 4,
    zero_tol = 1e-12,
    
    # ============================================================
    # Optional viral-load overlap
    # ============================================================
    use_viral_overlap = FALSE,
    overlap_method = c("quantile", "minmax"),
    overlap_quantiles = c(0.05, 0.95),
    min_infected_cells_per_patient = 20,
    
    # ============================================================
    # Expression
    # ============================================================
    assay = "RNA",
    layer = "counts",
    
    # ============================================================
    # Sample filtering
    # ============================================================
    min_cells_per_sample = 10,
    min_patients = 3,
    
    # ============================================================
    # Gene filtering
    # ============================================================
    min_count = 10,
    min_samples = 3,
    viral_genes = character(0),  # Explicit viral gene symbols to exclude
    
    # ============================================================
    # DESeq2
    # ============================================================
    fit_type = "parametric",
    sf_type = "ratio",
    alpha = 0.05,
    
    # ============================================================
    # DEG thresholds
    # ============================================================
    padj_cut = 0.05,
    logfc_cut = 0.5,
    
    # ============================================================
    # Analysis
    # ============================================================
    run_joint_model = FALSE,
    
    # ============================================================
    # Hallmark GSEA
    # ============================================================
    run_gsea = TRUE,
    species = c("human", "mouse", "other"),
    gmt_path = NULL,
    ranks_from = c("avg_log2FC", "pval_adj_signed"),
    top_up = 10,
    top_down = 10,
    drop_prefix = TRUE,
    title_prefix = "",
    
    # ============================================================
    # Visualization
    # ============================================================
    plot_volcano = TRUE,
    label_volcano = TRUE,
    volcano_top_up = 10,
    volcano_top_down = 10,
    volcano_label_size = 3.5,
    plot_pca = TRUE,
    plot_gsea = TRUE,
    plot_overlap = FALSE,
    plot_level_distribution = TRUE,
    
    col_up = "#D74C4C",
    col_down = "#3B6FB3",
    
    # ============================================================
    # Output
    # ============================================================
    save_results = FALSE,
    outdir = "Pseudobulk_results"
) {
  
  # ============================================================
  # 1. Basic checks
  # ============================================================
  
  overlap_method <- match.arg(overlap_method)
  ranks_from <- match.arg(ranks_from)
  species <- match.arg(species)
  
  if (!inherits(obj, "Seurat")) {
    stop("'obj' must be a Seurat object.")
  }
  
  required_packages <- c(
    "SeuratObject",
    "DESeq2",
    "Matrix",
    "ggplot2"
  )
  
  if (run_gsea) {
    required_packages <- c(
      required_packages,
      "clusterProfiler"
    )
  }
  
  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]
  
  if (length(missing_packages) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing_packages, collapse = ", ")
    )
  }
  
  if (!assay %in% SeuratObject::Assays(obj)) {
    stop("Assay not found: ", assay)
  }
  
  md <- obj[[]]
  
  required_cols <- if (use_viral_overlap) {
    c(patient_col, viral_col)
  } else {
    c(patient_col, level_col)
  }
  
  missing_cols <- setdiff(
    required_cols,
    colnames(md)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing metadata columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!is.null(celltype_include) &&
      !cell_type %in% colnames(md)) {
    stop("Cell-type column not found: ", cell_type)
  }
  
  if (min_patients < 2) {
    stop("min_patients must be >= 2.")
  }
  
  if (min_cells_per_sample < 1) {
    stop("min_cells_per_sample must be >= 1.")
  }
  
  if (min_count < 0 || min_samples < 1) {
    stop("Invalid gene-filtering thresholds.")
  }
  
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1.")
  }
  
  if (padj_cut <= 0 || padj_cut > 1) {
    stop("padj_cut must be in (0, 1].")
  }
  
  if (logfc_cut < 0) {
    stop("logfc_cut must be non-negative.")
  }
  
  if (zero_tol < 0) {
    stop("zero_tol must be non-negative.")
  }
  
  if (!sf_type %in% c(
    "ratio",
    "poscounts",
    "iterate"
  )) {
    stop("Invalid sf_type.")
  }
  
  if (!fit_type %in% c(
    "parametric",
    "local",
    "mean",
    "glmGamPoi"
  )) {
    stop("Invalid fit_type.")
  }
  
  if (top_up < 0 || top_down < 0) {
    stop("top_up and top_down must be non-negative.")
  }

  if (label_volcano && plot_volcano &&
      !requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Please install ggrepel to label volcano plots.")
  }

  if (any(!is.finite(c(volcano_top_up, volcano_top_down))) ||
      any(c(volcano_top_up, volcano_top_down) < 0) ||
      any(c(volcano_top_up, volcano_top_down) !=
          floor(c(volcano_top_up, volcano_top_down)))) {
    stop("volcano_top_up and volcano_top_down must be non-negative integers.")
  }

  if (!is.finite(volcano_label_size) || volcano_label_size <= 0) {
    stop("volcano_label_size must be positive.")
  }
  
  ref_chr <- as.character(ref_level)
  
  case_chr <- unique(
    as.character(case_levels)
  )
  
  case_chr <- setdiff(
    case_chr,
    ref_chr
  )
  
  if (length(case_chr) == 0) {
    stop("No valid case levels.")
  }
  
  if (use_viral_overlap) {
    
    if (n_levels < 1 ||
        length(n_levels) != 1 ||
        n_levels != as.integer(n_levels)) {
      stop("n_levels must be a positive integer.")
    }
    
    if (length(case_chr) != n_levels ||
        !setequal(
          case_chr,
          as.character(seq_len(n_levels))
        )) {
      
      stop(
        "When use_viral_overlap = TRUE, ",
        "case_levels must equal 1:n_levels."
      )
    }
    
    if (ref_chr %in% as.character(seq_len(n_levels))) {
      stop(
        "ref_level must not overlap with 1:n_levels."
      )
    }
    
    if (min_infected_cells_per_patient < 1) {
      stop(
        "min_infected_cells_per_patient must be >= 1."
      )
    }
    
    if (length(overlap_quantiles) != 2 ||
        anyNA(overlap_quantiles) ||
        overlap_quantiles[1] < 0 ||
        overlap_quantiles[2] > 1 ||
        overlap_quantiles[1] >= overlap_quantiles[2]) {
      
      stop(
        "overlap_quantiles must contain two ",
        "increasing values between 0 and 1."
      )
    }
    
    case_chr <- as.character(seq_len(n_levels))
  }
  
  analysis_mode <- if (use_viral_overlap) {
    "Common_viral_load_overlap"
  } else {
    "Original_infection_levels"
  }
  
  message("Analysis mode: ", analysis_mode)
  
  # ============================================================
  # 2. Load Hallmark gene sets
  # ============================================================
  
  hallmark <- NULL
  
  if (run_gsea) {
    
    message("Loading Hallmark gene sets...")
    
    user_gmt <- !is.null(gmt_path) &&
      length(gmt_path) == 1 &&
      !is.na(gmt_path) &&
      nzchar(gmt_path)
    
    if (!user_gmt) {
      
      if (species == "human") {
        
        gmt_path <- system.file(
          "h.all.v2025.1.Hs.symbols.gmt",
          package = "SingleViral"
        )
        
      } else if (species == "mouse") {
        
        gmt_path <- system.file(
          "mh.all.v2025.1.Mm.symbols.gmt",
          package = "SingleViral"
        )
        
      } else {
        
        stop(
          "For species = 'other', ",
          "you must provide a valid gmt_path."
        )
      }
    }
    
    if (nzchar(gmt_path) && file.exists(gmt_path)) {
      
      message("Using GMT file from: ", gmt_path)
      
      hallmark <- clusterProfiler::read.gmt(
        gmt_path
      )
      
      colnames(hallmark)[1:2] <- c(
        "gs_name",
        "gene_symbol"
      )
      
      hallmark <- hallmark[
        ,
        c("gs_name", "gene_symbol"),
        drop = FALSE
      ]
      
    } else {
      
      if (user_gmt) {
        stop("Specified GMT file not found: ", gmt_path)
      }
      
      if (!requireNamespace(
        "msigdbr",
        quietly = TRUE
      )) {
        stop(
          "GMT file not found and msigdbr ",
          "is not installed."
        )
      }
      
      sp_name <- if (species == "human") {
        "Homo sapiens"
      } else {
        "Mus musculus"
      }
      
      message(
        "GMT not found. Falling back to ",
        "msigdbr Hallmark for ",
        sp_name
      )
      
      hallmark <- msigdbr::msigdbr(
        species = sp_name,
        category = "H"
      )
      
      hallmark <- unique(
        hallmark[
          ,
          c("gs_name", "gene_symbol"),
          drop = FALSE
        ]
      )
    }
    
    hallmark <- unique(hallmark)
    
    hallmark <- hallmark[
      !is.na(hallmark$gene_symbol) &
        hallmark$gene_symbol != "",
      ,
      drop = FALSE
    ]
    
    if (nrow(hallmark) == 0) {
      stop("Hallmark gene sets are empty.")
    }
  }
  
  # ============================================================
  # 3. Extract raw UMI counts
  # ============================================================
  
  message("Extracting raw UMI counts...")
  
  available_layers <- SeuratObject::Layers(
    obj[[assay]],
    search = NA
  )
  
  split_layers <- grep(
    paste0("^", layer, "\\."),
    available_layers,
    value = TRUE
  )
  
  if (length(split_layers) > 0) {
    stop(
      "Split count layers detected. ",
      "Please run SeuratObject::JoinLayers() first."
    )
  }
  
  if (!layer %in% available_layers) {
    stop("Layer not found: ", layer)
  }
  
  counts <- SeuratObject::LayerData(
    object = obj,
    assay = assay,
    layer = layer
  )
  
  if (!inherits(counts, "sparseMatrix")) {
    counts <- Matrix::Matrix(
      counts,
      sparse = TRUE
    )
  }
  
  if (is.null(rownames(counts)) ||
      is.null(colnames(counts))) {
    stop("Count matrix must have gene and cell names.")
  }
  
  if (anyNA(counts@x) ||
      any(!is.finite(counts@x)) ||
      any(counts@x < 0) ||
      any(abs(counts@x - round(counts@x)) > 1e-8)) {
    
    stop(
      "Pseudobulk requires non-negative ",
      "integer raw counts."
    )
  }
  
  # Exclude explicitly specified viral genes from host-expression analyses.
  # Keep the original counts matrix intact for cell/sample aggregation.
  if (!is.character(viral_genes) || anyNA(viral_genes)) {
    stop("viral_genes must be a character vector without NA values.")
  }
  viral_genes <- unique(viral_genes[nzchar(viral_genes)])
  excluded_viral_genes <- intersect(rownames(counts), viral_genes)
  if (length(viral_genes) > 0) {
    message("Viral genes matched and excluded from host analyses: ",
            length(excluded_viral_genes))
    if (length(excluded_viral_genes) == 0) {
      warning("None of the supplied viral_genes matched count-matrix row names.")
    }
  } else {
    message("No viral_genes supplied: no viral genes will be excluded automatically.")
  }

  # ============================================================
  # 4. Prepare metadata
  # ============================================================
  
  md <- md[
    colnames(counts),
    ,
    drop = FALSE
  ]
  
  md$pb_patient <- as.character(
    md[[patient_col]]
  )
  
  if (viral_col %in% colnames(md)) {
    
    md$pb_viral_load <- suppressWarnings(
      as.numeric(as.character(md[[viral_col]]))
    )
    
  } else {
    
    md$pb_viral_load <- NA_real_
  }
  
  if (level_col %in% colnames(md)) {
    
    md$pb_original_level <- as.character(
      md[[level_col]]
    )
    
  } else {
    
    md$pb_original_level <- NA_character_
  }
  
  keep <- !is.na(md$pb_patient) &
    md$pb_patient != ""
  
  if (!is.null(celltype_include)) {
    
    keep <- keep &
      as.character(md[[cell_type]]) %in%
      as.character(celltype_include)
  }
  
  # ============================================================
  # 5A. Original infection-level mode
  # ============================================================
  
  patient_ranges <- NULL
  overlap_interval <- NULL
  retention_summary <- NULL
  level_boundaries <- NULL
  md_before_overlap <- NULL
  
  if (!use_viral_overlap) {
    
    message(
      "Using original infection levels..."
    )
    
    keep <- keep &
      !is.na(md$pb_original_level) &
      md$pb_original_level %in%
      c(ref_chr, case_chr)
    
    keep[is.na(keep)] <- FALSE
    
    md <- md[
      keep,
      ,
      drop = FALSE
    ]
    
    if (nrow(md) == 0) {
      stop(
        "No cells remain after filtering ",
        "original infection levels."
      )
    }
    
    md$pb_level <- md$pb_original_level
    
    md$pb_infected <- md$pb_level != ref_chr
    
    md$pb_in_common_range <- NA
    
    md$pb_keep <- TRUE
    
    message(
      "Cells retained: ",
      nrow(md)
    )
    
  } else {
    
    # ==========================================================
    # 5B. Common viral-load overlap mode
    # ==========================================================
    
    message(
      "Calculating common viral-load interval..."
    )
    
    keep <- keep &
      is.finite(md$pb_viral_load) &
      md$pb_viral_load >= 0
    
    keep[is.na(keep)] <- FALSE
    
    md <- md[
      keep,
      ,
      drop = FALSE
    ]
    
    if (nrow(md) == 0) {
      stop("No valid cells remain.")
    }
    
    md$pb_infected <- md$pb_viral_load > zero_tol
    
    patients_all <- sort(
      unique(md$pb_patient)
    )
    
    # --------------------------------------------------------
    # Patient-specific viral-load ranges
    # --------------------------------------------------------
    
    patient_ranges <- lapply(
      patients_all,
      function(p) {
        
        x <- md$pb_viral_load[
          md$pb_patient == p &
            md$pb_infected
        ]
        
        x <- x[is.finite(x)]
        
        n_inf <- length(x)
        
        if (n_inf == 0) {
          
          return(data.frame(
            patient_id = p,
            n_infected = 0L,
            n_uninfected = sum(
              md$pb_patient == p &
                !md$pb_infected
            ),
            viral_min = NA_real_,
            viral_max = NA_real_,
            viral_q_low = NA_real_,
            viral_q_high = NA_real_,
            range_low = NA_real_,
            range_high = NA_real_,
            eligible = FALSE,
            stringsAsFactors = FALSE
          ))
        }
        
        q <- stats::quantile(
          x,
          probs = overlap_quantiles,
          na.rm = TRUE,
          names = FALSE
        )
        
        low <- if (overlap_method == "quantile") {
          q[1]
        } else {
          min(x)
        }
        
        high <- if (overlap_method == "quantile") {
          q[2]
        } else {
          max(x)
        }
        
        data.frame(
          patient_id = p,
          n_infected = n_inf,
          n_uninfected = sum(
            md$pb_patient == p &
              !md$pb_infected
          ),
          viral_min = min(x),
          viral_max = max(x),
          viral_q_low = q[1],
          viral_q_high = q[2],
          range_low = low,
          range_high = high,
          eligible =
            n_inf >= min_infected_cells_per_patient &&
            low < high,
          stringsAsFactors = FALSE
        )
      }
    )
    
    patient_ranges <- do.call(
      rbind,
      patient_ranges
    )
    
    rownames(patient_ranges) <- NULL
    
    eligible_patients <- patient_ranges$patient_id[
      patient_ranges$eligible
    ]
    
    if (length(eligible_patients) < min_patients) {
      
      stop(
        "Insufficient eligible patients. ",
        "Found: ",
        length(eligible_patients),
        "; required: ",
        min_patients
      )
    }
    
    # --------------------------------------------------------
    # Common viral-load interval
    # --------------------------------------------------------
    
    ranges_use <- patient_ranges[
      patient_ranges$eligible,
      ,
      drop = FALSE
    ]
    
    common_low <- max(
      ranges_use$range_low
    )
    
    common_high <- min(
      ranges_use$range_high
    )
    
    if (!is.finite(common_low) ||
        !is.finite(common_high) ||
        common_low >= common_high) {
      
      stop(
        "No common viral-load interval exists ",
        "across eligible patients."
      )
    }
    
    overlap_interval <- data.frame(
      method = overlap_method,
      lower = common_low,
      upper = common_high,
      n_patients = length(eligible_patients),
      stringsAsFactors = FALSE
    )
    
    message(
      "Common viral-load interval: [",
      signif(common_low, 5),
      ", ",
      signif(common_high, 5),
      "]"
    )
    
    # --------------------------------------------------------
    # Retain common-range infected cells
    # --------------------------------------------------------
    
    md$pb_in_common_range <- md$pb_infected &
      md$pb_viral_load >= common_low &
      md$pb_viral_load <= common_high
    
    md$pb_keep <- md$pb_patient %in%
      eligible_patients &
      (
        !md$pb_infected |
          md$pb_in_common_range
      )
    
    retention_summary <- lapply(
      eligible_patients,
      function(p) {
        
        idx <- md$pb_patient == p
        
        n_inf <- sum(
          md$pb_infected[idx]
        )
        
        n_common <- sum(
          md$pb_in_common_range[idx]
        )
        
        data.frame(
          patient_id = p,
          infected_before = n_inf,
          infected_after = n_common,
          retention_rate = if (n_inf > 0) {
            n_common / n_inf
          } else {
            NA_real_
          },
          uninfected_cells = sum(
            !md$pb_infected[idx]
          ),
          stringsAsFactors = FALSE
        )
      }
    )
    
    retention_summary <- do.call(
      rbind,
      retention_summary
    )
    
    rownames(retention_summary) <- NULL
    
    md_before_overlap <- md
    
    md <- md[
      md$pb_keep,
      ,
      drop = FALSE
    ]
    
    if (nrow(md) == 0) {
      stop(
        "No cells remain after overlap filtering."
      )
    }
    
    # --------------------------------------------------------
    # Global quantile reclassification
    # --------------------------------------------------------
    
    message(
      "Reclassifying infected cells..."
    )
    
    infected_idx <- md$pb_infected
    
    infected_load <- md$pb_viral_load[
      infected_idx
    ]
    
    if (length(infected_load) < n_levels) {
      stop(
        "Too few infected cells for reclassification."
      )
    }
    
    probs <- seq(
      0,
      1,
      length.out = n_levels + 1
    )
    
    breaks <- as.numeric(
      stats::quantile(
        infected_load,
        probs = probs,
        na.rm = TRUE,
        names = FALSE
      )
    )
    
    if (length(unique(breaks)) <
        n_levels + 1) {
      
      stop(
        "Duplicated viral-load quantile boundaries. ",
        "Reduce n_levels or inspect viral_load."
      )
    }
    
    md$pb_level <- ref_chr
    
    md$pb_level[infected_idx] <- as.character(
      cut(
        infected_load,
        breaks = breaks,
        labels = seq_len(n_levels),
        include.lowest = TRUE,
        right = TRUE
      )
    )
    
    if (anyNA(md$pb_level)) {
      stop(
        "NA infection levels detected."
      )
    }
    
    level_boundaries <- data.frame(
      infection_level = case_chr,
      lower = breaks[-length(breaks)],
      upper = breaks[-1],
      stringsAsFactors = FALSE
    )
    
    message(
      "New infection-level boundaries:"
    )
    
    print(level_boundaries)
    
    message(
      "Patient-specific infected-cell retention:"
    )
    
    print(retention_summary)
  }
  
  # ============================================================
  # 6. Update count matrix
  # ============================================================
  
  counts <- counts[
    ,
    rownames(md),
    drop = FALSE
  ]
  
  if (ncol(counts) == 0) {
    stop("No cells remain after filtering.")
  }
  
  # ============================================================
  # 7. Infection-level distribution
  # ============================================================
  
  patients_use <- sort(
    unique(md$pb_patient)
  )
  
  level_distribution <- as.data.frame(
    table(
      patient_id = factor(
        md$pb_patient,
        levels = patients_use
      ),
      infection_level = factor(
        md$pb_level,
        levels = c(ref_chr, case_chr)
      )
    )
  )
  
  colnames(level_distribution)[3] <- "n_cells"
  
  level_distribution$percentage <- ave(
    level_distribution$n_cells,
    level_distribution$patient_id,
    FUN = function(x) {
      
      if (sum(x) > 0) {
        100 * x / sum(x)
      } else {
        rep(0, length(x))
      }
    }
  )
  
  message("Infection-level distribution:")
  
  print(level_distribution)
  
  # ============================================================
  # 8. Construct pseudobulk samples
  # ============================================================
  
  message("Constructing pseudobulk samples...")
  
  sample_key <- unique(
    md[
      ,
      c("pb_patient", "pb_level"),
      drop = FALSE
    ]
  )
  
  sample_key$sample_id <- paste0(
    "PB_",
    seq_len(nrow(sample_key))
  )
  
  key <- paste(
    sample_key$pb_patient,
    sample_key$pb_level,
    sep = "\r"
  )
  
  cell_key <- paste(
    md$pb_patient,
    md$pb_level,
    sep = "\r"
  )
  
  md$sample_id <- sample_key$sample_id[
    match(cell_key, key)
  ]
  
  # ============================================================
  # 9. Aggregate raw UMI counts
  # ============================================================
  
  sample_factor <- factor(
    md$sample_id,
    levels = sample_key$sample_id
  )
  
  aggregation_matrix <- Matrix::sparseMatrix(
    i = seq_len(ncol(counts)),
    j = as.integer(sample_factor),
    x = 1,
    dims = c(
      ncol(counts),
      nrow(sample_key)
    )
  )
  
  pb_counts <- counts %*% aggregation_matrix
  
  colnames(pb_counts) <- sample_key$sample_id
  
  # ============================================================
  # 10. Pseudobulk metadata and QC
  # ============================================================
  
  cell_numbers <- table(
    factor(
      md$sample_id,
      levels = sample_key$sample_id
    )
  )
  
  sample_key$n_cells <- as.integer(
    cell_numbers
  )
  
  sample_key$library_size <- as.numeric(
    Matrix::colSums(pb_counts)
  )
  
  sample_key$patient_id <- sample_key$pb_patient
  
  sample_key$infection_level <- sample_key$pb_level
  
  sample_key$mean_viral_load <- vapply(
    sample_key$sample_id,
    function(s) {
      
      x <- md$pb_viral_load[
        md$sample_id == s
      ]
      
      x <- x[is.finite(x)]
      
      if (length(x) == 0) {
        return(NA_real_)
      }
      
      mean(x)
    },
    numeric(1)
  )
  
  sample_key$median_viral_load <- vapply(
    sample_key$sample_id,
    function(s) {
      
      x <- md$pb_viral_load[
        md$sample_id == s
      ]
      
      x <- x[is.finite(x)]
      
      if (length(x) == 0) {
        return(NA_real_)
      }
      
      stats::median(x)
    },
    numeric(1)
  )
  
  rownames(sample_key) <- sample_key$sample_id
  
  sample_key$passed_qc <-
    sample_key$n_cells >= min_cells_per_sample &
    sample_key$library_size > 0
  
  sample_info_all <- sample_key
  
  sample_key <- sample_key[
    sample_key$passed_qc,
    ,
    drop = FALSE
  ]
  
  pb_counts <- pb_counts[
    ,
    sample_key$sample_id,
    drop = FALSE
  ]
  
  if (ncol(pb_counts) == 0) {
    stop(
      "No pseudobulk samples passed QC."
    )
  }
  
  message(
    "Pseudobulk samples retained: ",
    ncol(pb_counts)
  )
  
  message(
    "Patients retained: ",
    length(unique(sample_key$patient_id))
  )
  
  # ============================================================
  # 11. Gene filtering helper
  # ============================================================
  
  filter_genes <- function(mat) {

    # Remove viral genes before DESeq2 normalization, testing, GSEA,
    # and PCA, rather than merely hiding them in volcano plots.
    if (length(excluded_viral_genes) > 0) {
      mat <- mat[!rownames(mat) %in% excluded_viral_genes, , drop = FALSE]
    }

    keep_gene <- Matrix::rowSums(
      mat >= min_count
    ) >= min_samples
    
    mat[
      keep_gene,
      ,
      drop = FALSE
    ]
  }
  
  # ============================================================
  # 12. DESeq2 helper
  # ============================================================
  
  run_deseq <- function(
    count_mat,
    sample_md,
    reference,
    target
  ) {
    
    sample_md$patient_id <- factor(
      sample_md$patient_id
    )
    
    sample_md$infection_level <- factor(
      sample_md$infection_level,
      levels = c(reference, target)
    )
    
    count_mat <- filter_genes(
      count_mat
    )
    
    if (nrow(count_mat) < 2) {
      
      warning(
        "Too few genes after filtering."
      )
      
      return(NULL)
    }
    
    count_mat <- as.matrix(
      count_mat
    )
    
    storage.mode(count_mat) <- "integer"
    
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = count_mat,
      colData = sample_md,
      design = ~ patient_id + infection_level
    )
    
    # --------------------------------------------------------
    # Size-factor estimation
    # --------------------------------------------------------
    
    sf_use <- sf_type
    
    if (sf_use == "ratio") {
      
      genes_all_positive <- rowSums(
        count_mat > 0
      ) == ncol(count_mat)
      
      if (!any(genes_all_positive)) {
        
        message(
          "No genes positive in all samples; ",
          "using poscounts size-factor estimation."
        )
        
        sf_use <- "poscounts"
      }
    }
    
    # --------------------------------------------------------
    # Run DESeq2
    # --------------------------------------------------------
    
    dds <- DESeq2::DESeq(
      dds,
      fitType = fit_type,
      sfType = sf_use,
      quiet = TRUE
    )
    
    res <- DESeq2::results(
      dds,
      contrast = c(
        "infection_level",
        target,
        reference
      ),
      alpha = alpha
    )
    
    res_df <- as.data.frame(
      res
    )
    
    res_df$gene <- rownames(
      res_df
    )
    
    res_df <- res_df[
      order(
        res_df$padj,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]
    
    # --------------------------------------------------------
    # DEG classification
    # --------------------------------------------------------
    
    res_df$direction <- "NotSig"
    
    sig <- !is.na(res_df$padj) &
      !is.na(res_df$log2FoldChange) &
      res_df$padj < padj_cut
    
    res_df$direction[
      sig &
        res_df$log2FoldChange >= logfc_cut
    ] <- "Up"
    
    res_df$direction[
      sig &
        res_df$log2FoldChange <= -logfc_cut
    ] <- "Down"
    
    list(
      dds = dds,
      results = res_df,
      size_factor_method = sf_use
    )
  }
  
  # ============================================================
  # 13. Hallmark GSEA helper
  # Consistent with original InfectionDynamic
  # ============================================================
  
  run_gsea_analysis <- function(de) {
    
    if (is.null(de) || nrow(de) == 0) {
      return(NULL)
    }
    
    # --------------------------------------------------------
    # Build ranked gene list
    # --------------------------------------------------------
    
    if (ranks_from == "avg_log2FC") {
      
      ranks <- de$log2FoldChange
      
      names(ranks) <- de$gene
      
    } else {
      
      s <- sign(
        de$log2FoldChange
      )
      
      pa <- pmax(
        de$padj,
        1e-300
      )
      
      ranks <- -log10(pa) * s
      
      names(ranks) <- de$gene
    }
    
    ranks <- ranks[
      is.finite(ranks)
    ]
    
    ranks <- ranks[
      !is.na(names(ranks)) &
        names(ranks) != ""
    ]
    
    ranks <- ranks[
      !duplicated(names(ranks))
    ]
    
    ranks <- sort(
      ranks,
      decreasing = TRUE
    )
    
    if (length(ranks) < 10) {
      
      warning(
        "Too few valid genes for GSEA."
      )
      
      return(NULL)
    }
    
    message(
      "   Running GSEA (Hallmark)..."
    )
    
    gsea_res <- tryCatch({
      
      clusterProfiler::GSEA(
        geneList = ranks,
        TERM2GENE = hallmark,
        verbose = FALSE,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        minGSSize = 10,
        eps = 0
      )
      
    }, error = function(e) {
      
      warning(
        "GSEA error: ",
        conditionMessage(e)
      )
      
      NULL
    })
    
    gsea_res
  }
  
  # ============================================================
  # 14. Hallmark GSEA visualization
  # Consistent with original InfectionDynamic
  # ============================================================
  
  make_gsea_plot <- function(
    gsea_res,
    lvl,
    reference
  ) {
    
    if (is.null(gsea_res) ||
        !inherits(gsea_res, "gseaResult")) {
      return(NULL)
    }
    
    df <- tryCatch(
      as.data.frame(gsea_res),
      error = function(e) NULL
    )
    
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }
    
    if (!"Description" %in% colnames(df) &&
        "ID" %in% colnames(df)) {
      
      df$Description <- df$ID
    }
    
    if (drop_prefix) {
      
      df$Description <- gsub(
        "^HALLMARK_",
        "",
        df$Description
      )
    }
    
    df <- df[
      is.finite(df$NES),
      ,
      drop = FALSE
    ]
    
    # --------------------------------------------------------
    # Upregulated pathways
    # --------------------------------------------------------
    
    up_df <- df[
      df$NES > 0,
      ,
      drop = FALSE
    ]
    
    up_df <- up_df[
      order(
        up_df$NES,
        decreasing = TRUE
      ),
      ,
      drop = FALSE
    ]
    
    up_df <- utils::head(
      up_df,
      top_up
    )
    
    # --------------------------------------------------------
    # Downregulated pathways
    # --------------------------------------------------------
    
    down_df <- df[
      df$NES < 0,
      ,
      drop = FALSE
    ]
    
    down_df <- down_df[
      order(down_df$NES),
      ,
      drop = FALSE
    ]
    
    down_df <- utils::head(
      down_df,
      top_down
    )
    
    # --------------------------------------------------------
    # Combine pathways
    # --------------------------------------------------------
    
    plot_df_list <- list()
    
    if (nrow(up_df) > 0) {
      
      plot_df_list <- c(
        plot_df_list,
        list(
          transform(
            up_df,
            direction = "Up"
          )
        )
      )
    }
    
    if (nrow(down_df) > 0) {
      
      plot_df_list <- c(
        plot_df_list,
        list(
          transform(
            down_df,
            direction = "Down"
          )
        )
      )
    }
    
    if (length(plot_df_list) == 0) {
      return(NULL)
    }
    
    plot_df <- do.call(
      rbind,
      plot_df_list
    )
    
    plot_df$Description <- factor(
      plot_df$Description,
      levels = rev(
        c(
          up_df$Description,
          down_df$Description
        )
      )
    )
    
    # --------------------------------------------------------
    # GSEA barplot
    # --------------------------------------------------------
    
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = Description,
        y = NES,
        fill = factor(
          direction,
          levels = c("Up", "Down")
        )
      )
    ) +
      
      ggplot2::geom_col(
        width = 0.7
      ) +
      
      ggplot2::scale_fill_manual(
        values = c(
          Up = col_up,
          Down = col_down
        )
      ) +
      
      ggplot2::coord_flip() +
      
      ggplot2::geom_hline(
        yintercept = 0,
        linewidth = 0.3,
        linetype = "dashed"
      ) +
      
      ggplot2::geom_text(
        ggplot2::aes(
          label = round(NES, 2)
        ),
        hjust = ifelse(
          plot_df$NES > 0,
          -0.1,
          1.1
        ),
        size = 3.2,
        fontface = "bold"
      ) +
      
      ggplot2::labs(
        title = paste0(
          title_prefix,
          " Level ",
          lvl,
          " vs ",
          reference,
          " | Hallmark GSEA"
        ),
        x = "Pathway",
        y = "Normalized Enrichment Score (NES)"
      ) +
      
      ggplot2::theme_classic(
        base_size = 12
      ) +
      
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        ),
        axis.title.x = ggplot2::element_text(
          size = 12,
          face = "bold"
        ),
        axis.title.y = ggplot2::element_text(
          size = 12,
          face = "bold"
        ),
        axis.text.x = ggplot2::element_text(
          size = 12,
          face = "bold"
        ),
        axis.text.y = ggplot2::element_text(
          size = 12,
          face = "bold"
        ),
        legend.text = ggplot2::element_text(
          size = 12
        ),
        legend.title = ggplot2::element_blank()
      )
    
    p
  }
  
  # ============================================================
  # 15. Paired DESeq2 comparisons
  # ============================================================
  
  message(
    "Running paired DESeq2 comparisons..."
  )
  
  DE <- setNames(
    vector("list", length(case_chr)),
    case_chr
  )
  
  DDS <- DE
  GSEA <- DE
  paired_samples <- DE
  volcano_plots <- DE
  gsea_plots <- DE
  
  comparison_summary <- list()
  
  for (lvl in case_chr) {
    
    message(
      "\n===================================="
    )
    
    message(
      "Level ",
      lvl,
      " vs ",
      ref_chr
    )
    
    # --------------------------------------------------------
    # Identify paired patients
    # --------------------------------------------------------
    
    md_pair <- sample_key[
      sample_key$infection_level %in%
        c(ref_chr, lvl),
      ,
      drop = FALSE
    ]
    
    patient_tab <- table(
      md_pair$patient_id,
      factor(
        md_pair$infection_level,
        levels = c(ref_chr, lvl)
      )
    )
    
    valid_patients <- rownames(
      patient_tab
    )[
      rowSums(
        patient_tab > 0
      ) == 2
    ]
    
    n_patients <- length(
      valid_patients
    )
    
    message(
      "Paired patients: ",
      n_patients
    )
    
    comparison_summary[[lvl]] <- data.frame(
      comparison = paste0(
        "Level_",
        lvl,
        "_vs_",
        ref_chr
      ),
      n_patients = n_patients,
      status = "Skipped",
      n_up = NA_integer_,
      n_down = NA_integer_,
      stringsAsFactors = FALSE
    )
    
    if (n_patients < min_patients) {
      
      warning(
        "Skipping Level ",
        lvl,
        ": only ",
        n_patients,
        " paired patients."
      )
      
      next
    }
    
    # --------------------------------------------------------
    # Retain paired samples
    # --------------------------------------------------------
    
    md_pair <- md_pair[
      md_pair$patient_id %in%
        valid_patients,
      ,
      drop = FALSE
    ]
    
    md_pair <- md_pair[
      order(
        md_pair$patient_id,
        match(
          md_pair$infection_level,
          c(ref_chr, lvl)
        )
      ),
      ,
      drop = FALSE
    ]
    
    pair_counts <- pb_counts[
      ,
      md_pair$sample_id,
      drop = FALSE
    ]
    
    # --------------------------------------------------------
    # Run paired DESeq2
    # --------------------------------------------------------
    
    fit <- tryCatch({
      
      run_deseq(
        count_mat = pair_counts,
        sample_md = md_pair,
        reference = ref_chr,
        target = lvl
      )
      
    }, error = function(e) {
      
      warning(
        "DESeq2 failed for Level ",
        lvl,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    })
    
    if (is.null(fit)) {
      next
    }
    
    DE[[lvl]] <- fit$results
    
    DDS[[lvl]] <- fit$dds
    
    paired_samples[[lvl]] <- md_pair
    
    comparison_summary[[lvl]]$status <- "Completed"
    
    comparison_summary[[lvl]]$n_up <- sum(
      fit$results$direction == "Up"
    )
    
    comparison_summary[[lvl]]$n_down <- sum(
      fit$results$direction == "Down"
    )
    
    message(
      "Significant Up genes: ",
      comparison_summary[[lvl]]$n_up
    )
    
    message(
      "Significant Down genes: ",
      comparison_summary[[lvl]]$n_down
    )
    
    # --------------------------------------------------------
    # Hallmark GSEA
    # --------------------------------------------------------
    
    if (run_gsea) {
      
      GSEA[[lvl]] <- run_gsea_analysis(
        fit$results
      )
      
      if (plot_gsea) {
        
        gsea_plots[[lvl]] <- make_gsea_plot(
          gsea_res = GSEA[[lvl]],
          lvl = lvl,
          reference = ref_chr
        )
      }
    }
  }
  
  comparison_summary <- do.call(
    rbind,
    comparison_summary
  )
  
  rownames(comparison_summary) <- NULL
  
  # ============================================================
  # 16. Volcano plots
  # ============================================================
  
  if (plot_volcano) {
    
    message("Generating volcano plots...")
    
    for (lvl in case_chr) {
      
      de <- DE[[lvl]]
      
      if (is.null(de)) {
        next
      }
      
      de$logP <- -log10(
        pmax(
          de$padj,
          1e-300
        )
      )
      
      de$logP[
        is.na(de$logP)
      ] <- 0
      
      p <- ggplot2::ggplot(
        de,
        ggplot2::aes(
          x = log2FoldChange,
          y = logP,
          color = direction
        )
      ) +
        
        ggplot2::geom_point(
          size = 1,
          alpha = 0.7,
          na.rm = TRUE
        ) +
        
        ggplot2::scale_color_manual(
          values = c(
            Up = col_up,
            Down = col_down,
            NotSig = "grey75"
          )
        ) +
        
        ggplot2::geom_vline(
          xintercept = c(
            -logfc_cut,
            logfc_cut
          ),
          linetype = "dashed",
          linewidth = 0.3
        ) +
        
        ggplot2::geom_hline(
          yintercept = -log10(padj_cut),
          linetype = "dashed",
          linewidth = 0.3
        ) +
        
        ggplot2::labs(
          title = paste0(
            "Level ",
            lvl,
            " vs ",
            ref_chr,
            " | Paired Pseudobulk"
          ),
          x = "log2 Fold Change",
          y = expression(-log[10](FDR))
        ) +
        
        ggplot2::theme_classic(
          base_size = 12
        ) +
        
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            face = "bold",
            hjust = 0.5
          ),
          legend.title = ggplot2::element_blank()
        )
      
      # Label the most significant up/down genes (BH-adjusted P value).
      if (label_volcano) {
        label_candidates <- de[
          de$direction %in% c("Up", "Down") &
            is.finite(de$padj) &
            is.finite(de$log2FoldChange) &
            !is.na(de$gene) & de$gene != "",
          , drop = FALSE
        ]

        pick_labels <- function(direction, n) {
          x <- label_candidates[
            label_candidates$direction == direction,
            , drop = FALSE
          ]
          x <- x[order(x$padj, -abs(x$log2FoldChange)),
                 , drop = FALSE]
          utils::head(x, n)
        }

        label_df <- rbind(
          pick_labels("Up", volcano_top_up),
          pick_labels("Down", volcano_top_down)
        )

        if (nrow(label_df) > 0) {
          p <- p +
            ggrepel::geom_text_repel(
              data = label_df,
              ggplot2::aes(label = gene),
              size = volcano_label_size,
              max.overlaps = Inf,
              box.padding = 0.4,
              point.padding = 0.25,
              min.segment.length = 0,
              segment.color = "grey60",
              seed = 123,
              show.legend = FALSE
            ) +
            ggplot2::coord_cartesian(clip = "off") +
            ggplot2::theme(plot.margin =
              ggplot2::margin(12, 32, 12, 20))
        }
      }

      volcano_plots[[lvl]] <- p
    }
  }
  
  # ============================================================
  # 17. Pseudobulk PCA
  # ============================================================
  
  pca_plot <- NULL
  
  if (plot_pca && ncol(pb_counts) >= 3) {
    
    message("Generating pseudobulk PCA...")
    
    pca_counts <- filter_genes(
      pb_counts
    )
    
    if (nrow(pca_counts) >= 2) {
      
      lib <- Matrix::colSums(
        pca_counts
      )
      
      cpm <- t(
        t(pca_counts) / pmax(lib, 1)
      ) * 1e6
      
      log_cpm <- log2(
        as.matrix(cpm) + 1
      )
      
      gene_var <- apply(
        log_cpm,
        1,
        stats::var
      )
      
      top_genes <- names(
        sort(
          gene_var,
          decreasing = TRUE
        )
      )
      
      top_genes <- head(
        top_genes,
        2000
      )
      
      pca <- stats::prcomp(
        t(
          log_cpm[
            top_genes,
            ,
            drop = FALSE
          ]
        ),
        center = TRUE,
        scale. = FALSE
      )
      
      pca_df <- as.data.frame(
        pca$x[
          ,
          1:2,
          drop = FALSE
        ]
      )
      
      pca_df$sample_id <- rownames(
        pca_df
      )
      
      pca_df$patient_id <- sample_key$patient_id[
        match(
          pca_df$sample_id,
          sample_key$sample_id
        )
      ]
      
      pca_df$infection_level <- sample_key$infection_level[
        match(
          pca_df$sample_id,
          sample_key$sample_id
        )
      ]
      
      variance <- 100 *
        pca$sdev^2 /
        sum(pca$sdev^2)
      
      # Patient shapes
      patient_number <- length(
        unique(pca_df$patient_id)
      )
      
      if (patient_number <= 6) {
        
        pca_plot <- ggplot2::ggplot(
          pca_df,
          ggplot2::aes(
            x = PC1,
            y = PC2,
            color = infection_level,
            shape = patient_id
          )
        ) +
          ggplot2::geom_point(
            size = 3
          )
        
      } else {
        
        pca_plot <- ggplot2::ggplot(
          pca_df,
          ggplot2::aes(
            x = PC1,
            y = PC2,
            color = infection_level
          )
        ) +
          ggplot2::geom_point(
            size = 3
          )
      }
      
      pca_plot <- pca_plot +
        
        ggplot2::labs(
          title = "Pseudobulk PCA",
          x = paste0(
            "PC1 (",
            round(variance[1], 1),
            "%)"
          ),
          y = paste0(
            "PC2 (",
            round(variance[2], 1),
            "%)"
          )
        ) +
        
        ggplot2::theme_classic(
          base_size = 12
        )
    }
  }
  
  # ============================================================
  # 18. Viral-load overlap plot
  # ============================================================
  
  overlap_plot <- NULL
  
  if (use_viral_overlap && plot_overlap) {
    
    message(
      "Generating viral-load overlap plot..."
    )
    
    overlap_df <- md_before_overlap[
      md_before_overlap$pb_infected &
        md_before_overlap$pb_patient %in%
        eligible_patients,
      ,
      drop = FALSE
    ]
    
    overlap_plot <- ggplot2::ggplot(
      overlap_df,
      ggplot2::aes(
        x = pb_patient,
        y = pb_viral_load
      )
    ) +
      
      ggplot2::geom_violin(
        fill = "#4DBBD5",
        alpha = 0.5,
        trim = FALSE
      ) +
      
      ggplot2::geom_boxplot(
        width = 0.15,
        outlier.shape = NA,
        fill = "white"
      ) +
      
      ggplot2::geom_hline(
        yintercept = common_low,
        color = col_down,
        linetype = "dashed",
        linewidth = 0.7
      ) +
      
      ggplot2::geom_hline(
        yintercept = common_high,
        color = col_up,
        linetype = "dashed",
        linewidth = 0.7
      ) +
      
      ggplot2::labs(
        title = "Common Viral-load Interval",
        x = "Patient",
        y = "Viral Load"
      ) +
      
      ggplot2::theme_classic(
        base_size = 12
      ) +
      
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          hjust = 0.5
        )
      )
  }
  
  # ============================================================
  # 19. Infection-level distribution plot
  # ============================================================
  
  level_distribution_plot <- NULL
  
  if (plot_level_distribution) {
    
    level_colors <- c(
      "0" = "#BC80BD",
      "1" = "#4DBBD5",
      "2" = "#44C37F",
      "3" = "#F3A13B",
      "4" = "#E64B35"
    )
    
    requested_levels <- c(
      ref_chr,
      case_chr
    )
    
    if (!all(
      requested_levels %in% names(level_colors)
    )) {
      
      level_colors <- setNames(
        grDevices::hcl.colors(
          length(requested_levels),
          palette = "Dynamic"
        ),
        requested_levels
      )
    }
    
    plot_title <- if (use_viral_overlap) {
      "Reclassified Infection Levels"
    } else {
      "Original Infection Levels"
    }
    
    level_distribution_plot <- ggplot2::ggplot(
      level_distribution,
      ggplot2::aes(
        x = patient_id,
        y = percentage,
        fill = infection_level
      )
    ) +
      
      ggplot2::geom_col(
        width = 0.7,
        color = "white"
      ) +
      
      ggplot2::scale_fill_manual(
        values = level_colors,
        drop = FALSE
      ) +
      
      ggplot2::labs(
        title = plot_title,
        x = "Patient",
        y = "Percentage of Cells (%)",
        fill = "Infection Level"
      ) +
      
      ggplot2::theme_classic(
        base_size = 12
      ) +
      
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          hjust = 0.5
        )
      )
  }
  
  # ============================================================
  # 20. Optional joint DESeq2 model
  # ============================================================
  
  joint_dds <- NULL
  joint_results <- NULL
  
  if (run_joint_model) {
    
    message(
      "Running joint DESeq2 model..."
    )
    
    required_levels <- c(
      ref_chr,
      case_chr
    )
    
    joint_tab <- table(
      sample_key$patient_id,
      factor(
        sample_key$infection_level,
        levels = required_levels
      )
    )
    
    valid_joint <- rownames(
      joint_tab
    )[
      rowSums(
        joint_tab[
          ,
          required_levels,
          drop = FALSE
        ] > 0
      ) == length(required_levels)
    ]
    
    if (length(valid_joint) >= min_patients) {
      
      joint_md <- sample_key[
        sample_key$patient_id %in%
          valid_joint &
          sample_key$infection_level %in%
          required_levels,
        ,
        drop = FALSE
      ]
      
      joint_counts <- pb_counts[
        ,
        joint_md$sample_id,
        drop = FALSE
      ]
      
      joint_counts <- filter_genes(
        joint_counts
      )
      
      joint_dds <- tryCatch({
        
        if (nrow(joint_counts) < 2) {
          stop(
            "Too few genes after filtering."
          )
        }
        
        joint_md$patient_id <- factor(
          joint_md$patient_id
        )
        
        joint_md$infection_level <- factor(
          joint_md$infection_level,
          levels = required_levels
        )
        
        count_mat <- as.matrix(
          joint_counts
        )
        
        storage.mode(count_mat) <- "integer"
        
        dds <- DESeq2::DESeqDataSetFromMatrix(
          countData = count_mat,
          colData = joint_md,
          design = ~ patient_id + infection_level
        )
        
        sf_joint <- sf_type
        
        if (sf_joint == "ratio" &&
            !any(
              rowSums(count_mat > 0) ==
              ncol(count_mat)
            )) {
          
          sf_joint <- "poscounts"
        }
        
        DESeq2::DESeq(
          dds,
          fitType = fit_type,
          sfType = sf_joint,
          quiet = TRUE
        )
        
      }, error = function(e) {
        
        warning(
          "Joint DESeq2 failed: ",
          conditionMessage(e)
        )
        
        NULL
      })
      
      if (!is.null(joint_dds)) {
        
        joint_results <- list()
        
        for (lvl in case_chr) {
          
          res <- DESeq2::results(
            joint_dds,
            contrast = c(
              "infection_level",
              lvl,
              ref_chr
            ),
            alpha = alpha
          )
          
          res_df <- as.data.frame(
            res
          )
          
          res_df$gene <- rownames(
            res_df
          )
          
          joint_results[[lvl]] <- res_df
        }
      }
      
    } else {
      
      warning(
        "Joint model skipped: insufficient ",
        "patients with all requested levels."
      )
    }
  }
  
  # ============================================================
  # 21. Save results
  # ============================================================
  
  if (save_results) {
    
    message("Saving results...")
    
    dir.create(
      outdir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    save_csv <- function(
    x,
    filename
    ) {
      
      if (is.null(x)) {
        return(invisible(NULL))
      }
      
      utils::write.csv(
        x,
        file.path(outdir, filename),
        row.names = FALSE
      )
    }
    
    save_pdf <- function(
    p,
    filename,
    width = 8,
    height = 6
    ) {
      
      if (is.null(p)) {
        return(invisible(NULL))
      }
      
      ggplot2::ggsave(
        filename = file.path(
          outdir,
          filename
        ),
        plot = p,
        width = width,
        height = height,
        useDingbats = FALSE
      )
    }
    
    # --------------------------------------------------------
    # Analysis mode
    # --------------------------------------------------------
    
    save_csv(
      data.frame(
        analysis_mode = analysis_mode,
        use_viral_overlap = use_viral_overlap
      ),
      "Analysis_mode.csv"
    )
    
    # --------------------------------------------------------
    # Overlap-specific results
    # --------------------------------------------------------
    
    if (use_viral_overlap) {
      
      save_csv(
        patient_ranges,
        "Patient_viral_load_ranges.csv"
      )
      
      save_csv(
        overlap_interval,
        "Common_viral_load_interval.csv"
      )
      
      save_csv(
        retention_summary,
        "Patient_cell_retention.csv"
      )
      
      save_csv(
        level_boundaries,
        "New_infection_level_boundaries.csv"
      )
    }
    
    # --------------------------------------------------------
    # Infection-level distribution
    # --------------------------------------------------------
    
    save_csv(
      level_distribution,
      "Infection_level_distribution.csv"
    )
    
    # --------------------------------------------------------
    # Cell metadata
    # --------------------------------------------------------
    
    cell_metadata_export <- data.frame(
      cell_id = rownames(md),
      pb_patient = md$pb_patient,
      pb_viral_load = md$pb_viral_load,
      pb_original_level = md$pb_original_level,
      pb_level = md$pb_level,
      sample_id = md$sample_id,
      stringsAsFactors = FALSE
    )
    
    save_csv(
      cell_metadata_export,
      "Pseudobulk_cell_metadata.csv"
    )
    
    # --------------------------------------------------------
    # Pseudobulk metadata
    # --------------------------------------------------------
    
    save_csv(
      sample_info_all,
      "Pseudobulk_sample_metadata_all.csv"
    )
    
    save_csv(
      sample_key,
      "Pseudobulk_sample_metadata_filtered.csv"
    )
    
    # --------------------------------------------------------
    # Comparison summary
    # --------------------------------------------------------
    
    if (length(excluded_viral_genes) > 0) {
      save_csv(data.frame(gene = excluded_viral_genes),
               "Excluded_viral_genes.csv")
    }

    save_csv(
      comparison_summary,
      "Comparison_summary.csv"
    )
    
    # --------------------------------------------------------
    # Pseudobulk counts
    # --------------------------------------------------------
    
    saveRDS(
      pb_counts,
      file.path(
        outdir,
        "Pseudobulk_counts.rds"
      )
    )
    
    # --------------------------------------------------------
    # Each comparison
    # --------------------------------------------------------
    
    for (lvl in case_chr) {
      
      if (!is.null(DE[[lvl]])) {
        
        save_csv(
          DE[[lvl]],
          paste0(
            "DEG_Level_",
            lvl,
            "_vs_",
            ref_chr,
            ".csv"
          )
        )
      }
      
      if (!is.null(paired_samples[[lvl]])) {
        
        save_csv(
          paired_samples[[lvl]],
          paste0(
            "Paired_samples_Level_",
            lvl,
            "_vs_",
            ref_chr,
            ".csv"
          )
        )
      }
      
      if (!is.null(GSEA[[lvl]])) {
        
        save_csv(
          as.data.frame(GSEA[[lvl]]),
          paste0(
            "GSEA_Level_",
            lvl,
            "_vs_",
            ref_chr,
            ".csv"
          )
        )
      }
      
      save_pdf(
        volcano_plots[[lvl]],
        paste0(
          "Volcano_Level_",
          lvl,
          "_vs_",
          ref_chr,
          ".pdf"
        ),
        width = 7,
        height = 6
      )
      
      save_pdf(
        gsea_plots[[lvl]],
        paste0(
          "GSEA_Level_",
          lvl,
          "_vs_",
          ref_chr,
          "_TopUp",
          top_up,
          "_TopDown",
          top_down,
          ".pdf"
        ),
        width = 8,
        height = 6
      )
    }
    
    # --------------------------------------------------------
    # PCA
    # --------------------------------------------------------
    
    save_pdf(
      pca_plot,
      "Pseudobulk_PCA.pdf",
      width = 7,
      height = 6
    )
    
    # --------------------------------------------------------
    # Overlap plot
    # --------------------------------------------------------
    
    if (use_viral_overlap) {
      
      save_pdf(
        overlap_plot,
        "Common_viral_load_overlap.pdf",
        width = 8,
        height = 6
      )
    }
    
    # --------------------------------------------------------
    # Level distribution
    # --------------------------------------------------------
    
    save_pdf(
      level_distribution_plot,
      "Infection_level_distribution.pdf",
      width = 8,
      height = 6
    )
    
    # --------------------------------------------------------
    # Joint model
    # --------------------------------------------------------
    
    if (!is.null(joint_results)) {
      
      for (lvl in names(joint_results)) {
        
        save_csv(
          joint_results[[lvl]],
          paste0(
            "Joint_DEG_Level_",
            lvl,
            "_vs_",
            ref_chr,
            ".csv"
          )
        )
      }
    }
  }
  
  # ============================================================
  # 22. Return results
  # ============================================================
  
  message(
    "\nInfectionPseudobulk completed."
  )
  
  return(
    list(
      
      # Analysis mode
      analysis_mode = analysis_mode,
      use_viral_overlap = use_viral_overlap,
      
      # Viral-load overlap
      patient_ranges = patient_ranges,
      overlap_interval = overlap_interval,
      retention_summary = retention_summary,
      
      # Infection levels
      level_boundaries = level_boundaries,
      level_distribution = level_distribution,
      cell_metadata = md,
      
      # Pseudobulk
      pseudobulk_counts = pb_counts,
      sample_metadata = sample_key,
      sample_metadata_all = sample_info_all,
      
      # Genes excluded from host-expression analysis
      excluded_viral_genes = excluded_viral_genes,

      # Differential expression
      DE = DE,
      DDS = DDS,
      
      # Hallmark GSEA
      GSEA = GSEA,
      
      # Patient pairing
      paired_samples = paired_samples,
      comparison_summary = comparison_summary,
      
      # Visualization
      volcano_plots = volcano_plots,
      gsea_plots = gsea_plots,
      pca_plot = pca_plot,
      overlap_plot = overlap_plot,
      level_distribution_plot = level_distribution_plot,
      
      # Joint model
      joint_dds = joint_dds,
      joint_results = joint_results
    )
  )
}