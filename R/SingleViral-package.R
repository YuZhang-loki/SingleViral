#' SingleViral: A unified toolkit for viral infection analysis and visualization
#'
#' The **SingleViral** package provides a comprehensive workflow for quantifying viral
#' infection levels at the single-cell level, classifying infection stages,
#' identifying infection-associated gene signatures, and exploring host–virus
#' interactions through communication and functional enrichment analyses.
#'
#' ## Core Functions
#'
#' - **ViralloadScoring()**
#'   Compute viral load scores from viral gene expression counts using robust
#'   normalization, scaling, and winsorization procedures.
#'
#' - **TargetInfectionLevel()**
#'   Classify cells into multi-level infection stages (e.g., uninfected, low, moderate,
#'   high) based on viral load scores, using quantile or k-means thresholding.
#'
#' - **InfectionVariableGenes()**
#'   Identify differentially expressed genes between infected and uninfected cells,
#'   and visualize results via publication-quality volcano plots.
#'
#' - **InfectionDynamic()**
#'   Perform Hallmark GSEA or custom GMT-based enrichment analysis across infection
#'   levels to characterize pathway-level dynamics.
#'
#' - **CellVariableInfection()**
#'   Summarize infection proportions, viral load averages, and infection-level
#'   distributions across cell groups or timepoints.
#'
#' - **CellCycleByInfection()**
#'   Quantify and visualize cell-cycle phases (G1, S, G2M) in relation to infection
#'   levels using Seurat’s CellCycleScoring.
#'
#' - **CellChatByCellInfection()**
#'   Integrate CellChat analysis by infection status, generating cell–cell interaction
#'   networks, overall and per-pathway circular visualizations, and communication tables.
#'
#' ## Design Philosophy
#'
#' - Modular: each function can be used independently or as part of a full workflow
#' - Reproducible: integrates robust normalization, parameter transparency, and
#'   clear plotting functions
#' - Publication-ready: outputs ggplot2 figures and csv results formatted for
#'   downstream reporting or paper figures
#'
#' ## Typical Workflow
#'
#' 1. `ViralloadScoring()` → compute per-cell viral load
#' 2. `TargetInfectionLevel()` → classify infection stages
#' 3. `InfectionVariableGenes()` → find infection-associated DEGs
#' 4. `InfectionDynamic()` → run GSEA across levels
#' 5. `CellVariableInfection()` → summarize infection ratio per cell type
#' 6. `CellCycleByInfection()` → examine infection impact on cell cycle
#' 7. `CellChatByCellInfection()` → analyze communication networks
#'
#' ## Author
#' Developed by **[Your Name]**, 2025
#' Institute of Virology & Systems Biology
#'
#' ## Dependencies
#' This package imports functions from:
#' \pkg{Seurat}, \pkg{ggplot2}, \pkg{dplyr}, \pkg{Matrix},
#' \pkg{clusterProfiler}, \pkg{CellChat}, \pkg{future},
#' \pkg{ggrepel}, \pkg{readr}, and \pkg{scales}.
#'
#' @keywords internal
"_PACKAGE"
#' @aliases SingleViral-package
#' @keywords package
#' @import Seurat ggplot2 dplyr Matrix clusterProfiler future ggrepel readr scales
NULL
