#' Viral load scoring with robust baseline and optional gene weighting
#'
#' @description
#' Compute per-cell viral-load scores from raw counts in a Seurat object.
#' The pipeline applies CPM normalization + log1p, subtracts a per-gene
#' baseline (quantile or negative cells; `q_baseline=0` forces zero),
#' divides by a robust scale (MAD or IQR), truncates negatives to 0,
#' optionally winsorizes large values, and then aggregates per cell
#' (unweighted mean and weighted mean with optional gene weights/boost).
#'
#' @param obj A Seurat object (viral genes must be in the current assay rows).
#' @param viral_genes Character vector of viral gene symbols (must match rownames).
#' @param neg_cells Optional character vector of cell barcodes regarded as
#'   uninfected/negative; used to estimate per-gene baselines via row medians.
#'   If `NULL`, the baseline is the per-gene `q_baseline` quantile.
#' @param q_baseline Numeric in \[0,1\]. Per-gene quantile used as baseline
#'   when `neg_cells` is `NULL`. **If set to 0, the baseline is forced to 0.**
#' @param scale_type Robust scale type per gene: `"MAD"` (median absolute deviation)
#'   or `"IQR"` (interquartile range).
#' @param winsor_q Optional numeric in (0,1]; winsorization quantile per gene
#'   applied after non-negative truncation. Use `NULL` to disable winsorization.
#' @param epsilon Small positive floor for the robust scale to avoid division by zero.
#' @param important_genes Optional character vector; these genes will be
#'   additionally up-weighted by `boost`.
#' @param boost Numeric multiplier applied to `important_genes` weights (default 2).
#' @param gene_weights Optional named numeric vector of per-gene weights
#'   (names must match `viral_genes` present in the object). These are multiplied
#'   with the base weight (1) and with `boost` for `important_genes`.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{Z} Sparse/regular matrix of per-gene standardized, non-negative,
#'         (optionally winsorized) values (rows = genes, cols = cells).
#'   \item \code{weights} Named numeric vector of gene weights used.
#'   \item \code{viral_load_unweighted} Numeric vector: per-cell mean of \code{Z}.
#'   \item \code{viral_load_weighted} Numeric vector: per-cell weighted mean of \code{Z}.
#' }
#'
#' @examples
#' \dontrun{
#' # Minimal example
#' res <- ViralLoadScoring(
#'   obj = seurat_obj,
#'   viral_genes = c("U27","U41","U11"),
#'   q_baseline = 0,              # force zero baseline
#'   scale_type = "MAD",
#'   winsor_q = 0.995,
#'   important_genes = c("U27","U41"),
#'   boost = 2
#' )
#' seurat_obj$viral_load <- scale(res$viral_load_weighted)
#' }
#'
#' @export
#' @importFrom Seurat GetAssayData
#' @importFrom matrixStats rowMedians
#' @importFrom stats quantile mad
ViralLoadScoring <- function(obj, viral_genes, neg_cells = NULL,
                             q_baseline = 0, scale_type = c("MAD","IQR"),
                             winsor_q = 0.995, epsilon = 0.05,
                             important_genes = NULL, boost = 2, gene_weights = NULL) {
  
  scale_type <- match.arg(scale_type)
  genes <- intersect(viral_genes, rownames(obj))
  if (length(genes) == 0) stop("No viral genes found in object.")
  
  ## 1) CPM + log1p
  V_counts <- GetAssayData(obj, layer = "counts")[genes, , drop = FALSE]
  libsize  <- Matrix::colSums(V_counts) + 1
  CPM      <- t(t(V_counts) / libsize * 1e6)
  X        <- log1p(CPM)  # rows = genes, cols = cells
  
  ## 2) Baseline bg
  if (isTRUE(all.equal(q_baseline, 0))) {
    bg <- rep(0, nrow(X)); names(bg) <- rownames(X)
  } else if (!is.null(neg_cells)) {
    neg_cells <- intersect(colnames(X), neg_cells)
    if (length(neg_cells) == 0) stop("None of 'neg_cells' found among object cells.")
    bg <- matrixStats::rowMedians(X[, neg_cells, drop = FALSE])
  } else {
    bg <- apply(X, 1, function(v) stats::quantile(v, probs = q_baseline, names = FALSE, type = 7))
  }
  
  ## 3) Robust scale
  if (scale_type == "MAD") {
    sg <- apply(X, 1, stats::mad)
  } else {
    sg <- apply(X, 1, function(v) diff(stats::quantile(v, probs = c(0.25, 0.75), names = FALSE)))
  }
  sg <- pmax(sg, epsilon)
  
  ## 4) Standardize -> non-negative -> winsor
  X_center <- X - bg
  Z <- X_center / sg
  Z@x <- pmax(Z@x, 0)
  
  if (!is.null(winsor_q)) {
    caps <- apply(Z, 1, function(v) stats::quantile(v, probs = winsor_q, names = FALSE))
    for (g in seq_len(nrow(Z))) {
      idx <- which(Z[g, ] > caps[g])
      if (length(idx) > 0) Z[g, idx] <- caps[g]
    }
  }
  
  neg_cells_used <- character(0)
  if (!is.null(neg_cells)) {
    neg_cells_used <- intersect(colnames(Z), neg_cells)
    if (length(neg_cells_used) > 0) {
      Z[, neg_cells_used] <- 0  # 这样两种得分都会为 0
    }
  }
  
  ## 6) Gene weighting
  w <- rep(1, length(genes)); names(w) <- genes
  if (!is.null(gene_weights)) {
    if (is.null(names(gene_weights))) stop("'gene_weights' must be a named numeric vector.")
    hit <- intersect(names(gene_weights), genes)
    if (length(hit) > 0) w[hit] <- w[hit] * as.numeric(gene_weights[hit])
  }
  if (!is.null(important_genes)) {
    hit2 <- intersect(important_genes, genes)
    if (length(hit2) > 0) w[hit2] <- w[hit2] * boost
  }
  w[!is.finite(w)] <- 0
  if (sum(w) == 0) w[] <- 1
  w <- w / sum(w)
  
  ## 7) Scores
  viral_load_unweighted <- Matrix::colMeans(Z)
  viral_load_weighted   <- as.numeric(Matrix::t(Z) %*% w)

  if (length(neg_cells_used) > 0) {
    viral_load_unweighted[neg_cells_used] <- 0
    viral_load_weighted[neg_cells_used]   <- 0
  }
  
  list(
    Z = Z,
    weights = w,
    viral_load_unweighted = viral_load_unweighted,
    viral_load_weighted   = viral_load_weighted,
    neg_cells_used = neg_cells_used
  )
}
