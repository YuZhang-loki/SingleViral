#' Viral gene dynamics across infection levels
#'
#' This function filters cells by infection level and timepoint, extracts
#' the expression matrix of user-defined viral genes, computes the mean
#' expression of each viral gene within each infection level, transforms
#' the matrix as \code{log2(avg + 1)}, and visualizes it as a heatmap.
#' Optionally, the averaged expression matrix is written to CSV and the
#' heatmap is saved as a PDF.
#'
#' @param obj A \code{Seurat} object.
#' @param level_col Name of the metadata column storing infection levels
#'   (e.g. \code{"infection_level"}).
#' @param time_col Name of the metadata column storing timepoints
#'   (e.g. \code{"timepoint"}).
#' @param time_include Optional character vector of timepoints to keep.
#'   Cells whose \code{time_col} is not in this vector will be discarded.
#'   If \code{NULL}, all timepoints are used.
#' @param level_include Optional vector of infection levels to keep
#'   (character or numeric). If \code{NULL}, all levels present in
#'   \code{level_col} are used.
#' @param viral_genes Character vector of viral gene symbols to extract.
#'   Only genes present in the expression matrix rownames are used.
#' @param assay Name of the assay to use (e.g. \code{"RNA"}). If
#'   \code{NULL}, \code{DefaultAssay(obj)} is used.
#' @param slot Which expression slot to use, one of \code{"data"} or
#'   \code{"counts"}.
#' @param save_dir Optional output directory. If not \code{NULL}, a CSV
#'   file (\code{viral_expression_log2avg_by_level.csv}) containing the
#'   log2 averaged matrix and a PDF file
#'   (\code{viral_expression_heatmap.pdf}) with the heatmap are saved.
#' @param angle_col Optional angle (in degrees) for column labels in
#'   \code{pheatmap} (passed to \code{angle_col}).
#' @param fontsize Base font size for the heatmap.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{avg_expr}{A numeric matrix of dimension
#'   \code{length(viral_genes_found) x n_levels} containing
#'   \code{log2(avg + 1)} expression values (rows = genes, columns = infection
#'   levels).}
#'   \item{heatmap}{The \code{pheatmap} object returned by
#'   \code{pheatmap::pheatmap}.}
#' }
#'
#' @importFrom Seurat DefaultAssay Assays GetAssayData
#' @importFrom Matrix rowMeans
#' @importFrom pheatmap pheatmap
#' @export
VirusGeneDynamic <- function(
    obj,
    level_col     = "infection_level",
    time_col      = "timepoint",
    time_include  = NULL,
    level_include = NULL,
    viral_genes,
    assay         = NULL,
    slot          = c("data", "counts"),
    save_dir      = NULL,
    angle_col     = NULL,
    fontsize      = 7
){
  stopifnot(inherits(obj, "Seurat"))
  slot <- match.arg(slot)
  
  md <- obj@meta.data
  if (!level_col %in% colnames(md)) stop("Column '", level_col, "' not found.")
  if (!time_col  %in% colnames(md)) stop("Column '", time_col, "' not found.")
  
  # -----------------------------
  # 1. Filter by timepoint and infection level
  # -----------------------------
  md_use <- md
  
  if (!is.null(time_include)) {
    md_use <- md_use[as.character(md_use[[time_col]]) %in%
                       as.character(time_include), , drop = FALSE]
  }
  if (!is.null(level_include)) {
    md_use <- md_use[as.character(md_use[[level_col]]) %in%
                       as.character(level_include), , drop = FALSE]
  }
  
  if (nrow(md_use) == 0) {
    stop("No cells left after filtering by time_include / level_include.")
  }
  
  cells_use <- rownames(md_use)
  
  # Order infection levels as an ordered factor
  md_use[[level_col]] <- factor(
    as.character(md_use[[level_col]]),
    levels = sort(unique(as.character(md_use[[level_col]])))
  )
  
  # -----------------------------
  # 2. Extract viral gene expression matrix
  # -----------------------------
  if (is.null(assay)) assay <- Seurat::DefaultAssay(obj)
  if (!assay %in% Seurat::Assays(obj)) {
    stop("Assay '", assay, "' not found in object.")
  }
  
  mat_all <- Seurat::GetAssayData(obj, assay = assay, slot = slot)
  
  genes_use <- intersect(viral_genes, rownames(mat_all))
  if (length(genes_use) == 0) {
    warning("No viral genes found in object.")
    return(list(avg_expr = NULL, heatmap = NULL))
  }
  
  mat_sub <- mat_all[genes_use, cells_use, drop = FALSE]
  
  # -----------------------------
  # 3. Compute mean expression per infection level
  # -----------------------------
  lev_vec <- md_use[colnames(mat_sub), level_col]
  lev_vec <- factor(lev_vec, levels = levels(md_use[[level_col]]))
  
  avg_mat <- sapply(levels(lev_vec), function(lv) {
    if (sum(lev_vec == lv) == 0) {
      rep(NA_real_, length(genes_use))
    } else {
      Matrix::rowMeans(mat_sub[, lev_vec == lv, drop = FALSE])
    }
  })
  rownames(avg_mat) <- genes_use
  colnames(avg_mat) <- levels(lev_vec)
  
  # log2(avg + 1)
  avg_log_mat <- log2(avg_mat + 1)
  
  # -----------------------------
  # 4. Draw heatmap
  # -----------------------------
  suppressPackageStartupMessages(library(pheatmap))
  hm <- pheatmap::pheatmap(
    avg_log_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize     = fontsize,
    angle_col    = angle_col,
    main         = "Viral Gene Expression (log2 avg + 1)"
  )
  dev.off()
  # -----------------------------
  # 5. Save outputs
  # -----------------------------
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    
    utils::write.csv(
      avg_log_mat,
      file  = file.path(save_dir, "viral_expression_log2avg_by_level.csv"),
      quote = FALSE
    )
    
    grDevices::pdf(
      file.path(save_dir, "viral_expression_heatmap.pdf"),
      width  = 8,
      height = 8
    )
    pheatmap::pheatmap(
      avg_log_mat,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      fontsize     = fontsize,
      angle_col    = angle_col,
      main         = "Viral Gene Expression (log2 avg + 1)"
    )
    grDevices::dev.off()
  }
  
  return(list(avg_expr = avg_log_mat, heatmap = hm))
}