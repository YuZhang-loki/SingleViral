#' Viral load scoring with robust baseline, gene weighting and visualization
#'
#' @description
#' Compute per-cell viral load scores from raw counts in a Seurat object.
#' The pipeline applies CPM normalization + log1p, subtracts a per-gene
#' baseline (quantile or negative cells; \code{q_baseline = 0} forces zero),
#' divides by a robust scale (MAD or IQR), truncates negatives to 0,
#' optionally winsorizes large values, and then aggregates per cell
#' (unweighted mean and weighted mean with optional gene weights/boost).
#' Optionally, the chosen score vector is written to \code{obj$viral_load}
#' for downstream visualization within the function.
#'
#' @param obj A Seurat object; viral genes must be present in the current assay
#'   (i.e. row names of \code{GetAssayData(obj)}).
#' @param viral_genes Character vector of viral gene symbols. Must match
#'   row names of the current assay; genes not found are silently ignored.
#' @param neg_cells Optional character vector of cell barcodes regarded as
#'   uninfected/negative. If provided, per-gene baselines are estimated as
#'   row medians over these cells. If \code{NULL}, the baseline is the
#'   per-gene \code{q_baseline} quantile.
#' @param q_baseline Numeric in \code{[0,1]}. Per-gene quantile used as
#'   baseline when \code{neg_cells} is \code{NULL}. If set to \code{0},
#'   the baseline is forced to \code{0} for all genes.
#' @param scale_type Robust scale type per gene: \code{"MAD"} (median absolute
#'   deviation) or \code{"IQR"} (interquartile range).
#' @param winsor_q Optional numeric in \code{(0,1]}; winsorization quantile
#'   per gene applied after non-negative truncation. Use \code{NULL} to
#'   disable winsorization.
#' @param epsilon Small positive floor for the robust scale to avoid division
#'   by zero; scales smaller than \code{epsilon} are replaced by \code{epsilon}.
#' @param important_genes Optional character vector; these genes will be
#'   additionally up-weighted by \code{boost} in the weighted score.
#' @param boost Numeric multiplier applied to \code{important_genes} weights
#'   (default \code{2}).
#' @param gene_weights Optional named numeric vector of per-gene weights
#'   (names must match \code{viral_genes} present in the object). These are
#'   multiplied with the base weight (\code{1}) and with \code{boost} for
#'   \code{important_genes}.
#' @param plot Logical. If \code{TRUE}, generate FeaturePlot-based
#'   visualizations of the viral load score and return the ggplot objects.
#' @param write_mode Character, one of \code{"weighted"} or
#'   \code{"unweighted"}. Controls which score vector is written to
#'   \code{obj$viral_load} and used for internal visualization:
#'   \itemize{
#'     \item \code{"weighted"}: use \code{viral_load_weighted}
#'     \item \code{"unweighted"}: use \code{viral_load_unweighted}
#'   }
#' @param group.col Optional character scalar. Column name in
#'   \code{obj@meta.data} used as \code{split.by} in \code{FeaturePlot}
#'   for grouped visualization. If \code{NULL}, only the global plot is drawn.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{Z}: Matrix of per-gene standardized, non-negative,
#'         (optionally winsorized) values (rows = genes, cols = cells).
#'   \item \code{weights}: Named numeric vector of gene weights used.
#'   \item \code{viral_load_unweighted}: Numeric vector, per-cell mean of \code{Z}.
#'   \item \code{viral_load_weighted}: Numeric vector, per-cell weighted mean of \code{Z}.
#'   \item \code{neg_cells_used}: Character vector of negative cells that were set to 0.
#'   \item \code{write_mode}: The mode used to write scores to \code{obj$viral_load}.
#'   \item \code{group.col}: The grouping column used for split visualization (if any).
#'   \item \code{plots}: If \code{plot = TRUE}, a list of ggplot objects:
#'         \code{global} (overall FeaturePlot) and \code{split}
#'         (split-by-\code{group.col} FeaturePlot or \code{NULL}).
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with weighted scores and split by timepoint
#' res <- ViralLoadScoring(
#'   obj          = seurat_obj,
#'   viral_genes  = c("U27","U41","U11"),
#'   q_baseline   = 0,              # force zero baseline
#'   scale_type   = "MAD",
#'   winsor_q     = 0.995,
#'   important_genes = c("U27","U41"),
#'   boost        = 2,
#'   write_mode   = "weighted",
#'   group.col    = "timepoint",
#'   plot         = TRUE
#' )
#'
#' # Access numeric scores
#' head(res$viral_load_weighted)
#' head(res$viral_load_unweighted)
#'
#' # Use the returned ggplot objects
#' res$plots$global
#' res$plots$split
#' }
#'
#' @export
#' @importFrom Seurat GetAssayData FeaturePlot
#' @importFrom matrixStats rowMedians
#' @importFrom stats quantile mad
ViralLoadScoring <- function(
    obj,
    viral_genes,
    neg_cells = NULL,
    q_baseline = 0,
    scale_type = c("MAD","IQR"),
    winsor_q = 0.995,
    epsilon = 0.05,
    important_genes = NULL,
    boost = 2,
    gene_weights = NULL,
    plot = TRUE,
    write_mode = c("weighted", "unweighted"),
    group.col = NULL
) {
  write_mode <- match.arg(write_mode)
  scale_type <- match.arg(scale_type)
  
  genes <- intersect(viral_genes, rownames(obj))
  if (length(genes) == 0) stop("No viral genes found in object.")
  
  ## 1) CPM + log1p
  V_counts <- GetAssayData(obj, layer = "counts")[genes, , drop = FALSE]
  libsize  <- Matrix::colSums(V_counts) + 1
  CPM      <- t(t(V_counts) / libsize * 1e6)
  X        <- log1p(CPM)
  
  ## 2) Baseline
  if (isTRUE(all.equal(q_baseline, 0))) {
    bg <- rep(0, nrow(X)); names(bg) <- rownames(X)
  } else if (!is.null(neg_cells)) {
    neg_cells <- intersect(colnames(X), neg_cells)
    if (length(neg_cells) == 0) stop("None of 'neg_cells' found among object cells.")
    bg <- matrixStats::rowMedians(X[, neg_cells, drop = FALSE])
  } else {
    bg <- apply(X, 1, function(v)
      stats::quantile(v, probs = q_baseline, names = FALSE)
    )
  }
  
  ## 3) Robust scale
  if (scale_type == "MAD") {
    sg <- apply(X, 1, stats::mad)
  } else {
    sg <- apply(X, 1, function(v)
      diff(stats::quantile(v, probs = c(0.25, 0.75), names = FALSE))
    )
  }
  sg <- pmax(sg, epsilon)
  
  ## 4) Z-score + truncation + winsorization
  Z <- (X - bg) / sg
  Z@x <- pmax(Z@x, 0)
  
  if (!is.null(winsor_q)) {
    caps <- apply(Z, 1, function(v)
      stats::quantile(v, probs = winsor_q, names = FALSE)
    )
    for (g in seq_len(nrow(Z))) {
      idx <- which(Z[g, ] > caps[g])
      if (length(idx) > 0) Z[g, idx] <- caps[g]
    }
  }
  
  ## 5) Negative cells reset
  neg_cells_used <- character(0)
  if (!is.null(neg_cells)) {
    neg_cells_used <- intersect(colnames(Z), neg_cells)
    if (length(neg_cells_used) > 0) {
      Z[, neg_cells_used] <- 0
    }
  }
  
  ## 6) Gene weighting
  w <- rep(1, length(genes)); names(w) <- genes
  
  if (!is.null(gene_weights)) {
    if (is.null(names(gene_weights))) {
      stop("'gene_weights' must be a named numeric vector.")
    }
    hit <- intersect(names(gene_weights), genes)
    if (length(hit) > 0) {
      w[hit] <- w[hit] * gene_weights[hit]
    }
  }
  
  if (!is.null(important_genes)) {
    hit2 <- intersect(important_genes, genes)
    if (length(hit2) > 0) {
      w[hit2] <- w[hit2] * boost
    }
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
  
  ## 8) Write chosen score into obj$viral_load for internal plotting
  if (write_mode == "weighted") {
    obj$viral_load <- viral_load_weighted
  } else {
    obj$viral_load <- viral_load_unweighted
  }
  
  ## 9) Visualization with FeaturePlot
  plots <- NULL
  if (isTRUE(plot)) {
    p1 <- FeaturePlot(
      obj,
      features = "viral_load",
      cols = c("grey", "red")
    ) +
      ggtitle("Viral load scoring") +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
      )
    
    if (!is.null(group.col)) {
      if (!group.col %in% colnames(obj@meta.data)) {
        stop(paste0("group.col '", group.col, "' not found in obj@meta.data"))
      }
      p2 <- FeaturePlot(
        obj,
        features = "viral_load",
        cols = c("grey", "red"),
        split.by = group.col
      )
    } else {
      p2 <- NULL
    }
    
    plots <- list(global = p1, split = p2)
  }
  
  list(
    Z                      = Z,
    weights                = w,
    viral_load_unweighted  = viral_load_unweighted,
    viral_load_weighted    = viral_load_weighted,
    neg_cells_used         = neg_cells_used,
    write_mode             = write_mode,
    group.col              = group.col,
    plots                  = plots
  )
}