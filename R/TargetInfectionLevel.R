#' Discretize viral-load into infection levels and visualize with DimPlot
#'
#' @description
#' Assigns per-cell infection levels using a viral-load column in a Seurat object.
#' Cells with (absolute) score <= `zero_tol` are forced to level 0 ("Uninfected").
#' Positive cells are partitioned into `n_levels - 1` bins using either quantiles
#' or 1D k-means on positive scores only. Both the numeric level and a label are
#' written back to `meta.data`. Optionally, DimPlot visualizations are generated.
#'
#' @param obj A **Seurat** object.
#' @param score_col Character. Column in `obj@meta.data` that stores viral-load scores.
#' @param n_levels Integer (>=2). Total number of levels **including** level 0.
#'   Level 0 is reserved for uninfected cells.
#' @param method `"quantile"` or `"kmeans"`. Binning strategy for **positive** cells.
#' @param label_col Character. Name of the label column to write.
#' @param level_col Character. Name of the numeric level column to write.
#' @param labels Optional character vector of length `n_levels`. Custom labels where
#'   `labels[1]` corresponds to level 0 ("Uninfected").
#' @param zero_tol Numeric. Scores with `abs(score) <= zero_tol` are treated as zero.
#' @param kmeans_nstart Integer. `nstart` for k-means when `method = "kmeans"`.
#' @param plot Logical. Whether to generate DimPlot visualizations.
#' @param level_cols Named character vector. Custom colors for `infection_level`.
#'   If `NULL`, Seurat default colors are used.
#' @param label_cols Named character vector. Custom colors for `infection_label`.
#'   If `NULL`, Seurat default colors are used.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{object} Updated Seurat object.
#'   \item \code{thresholds} Numeric vector of bin thresholds for positives.
#'   \item \code{table} Frequency table of assigned levels.
#'   \item \code{levels} Integer vector of per-cell infection levels.
#'   \item \code{labels} Character vector of per-cell infection labels.
#'   \item \code{plots} A list with \code{level_plot} and \code{label_plot}
#'         (only returned when \code{plot = TRUE}).
#' }
#'
#' @examples
#' \dontrun{
#' res <- TargetInfectionLevel(
#'   obj       = data,
#'   score_col = "viral_load",
#'   n_levels  = 5,
#'   labels    = c("Uninfected", "Low infection", "Moderate infection",
#'                 "Moderate to high", "High infection"),
#'   plot = TRUE
#' )
#'
#' res$plots$level_plot
#' res$plots$label_plot
#' }
#'
#' @export
#' @importFrom Seurat DimPlot
TargetInfectionLevel <- function(obj,
                                 score_col     = "viral_load",
                                 n_levels      = 4,
                                 method        = c("quantile", "kmeans"),
                                 label_col     = "infection_label",
                                 level_col     = "infection_level",
                                 labels        = NULL,
                                 zero_tol      = 1e-12,
                                 kmeans_nstart = 10,
                                 plot          = TRUE,
                                 level_cols    = NULL,
                                 label_cols    = NULL) {
  
  stopifnot(inherits(obj, "Seurat"))
  method <- match.arg(method)
  if (n_levels < 2) stop("`n_levels` must be >= 2.")
  
  ## 1. Check whether the viral load column exists
  if (!(score_col %in% colnames(obj[[]]))) {
    stop(sprintf("Column '%s' not found in meta.data.", score_col))
  }
  
  scores <- obj[[score_col]][, 1]
  names(scores) <- colnames(obj)
  
  ## 2. Define uninfected cells (level 0)
  zero_mask <- abs(scores) <= zero_tol
  lvl <- rep(NA_integer_, length(scores))
  names(lvl) <- names(scores)
  lvl[zero_mask] <- 0
  
  pos_scores <- scores[!zero_mask]
  pos_idx <- names(pos_scores)
  
  ## 3. Stratify positive cells into 1...(n_levels-1)
  n_pos_bins <- n_levels - 1
  thresholds <- NULL
  
  if (length(pos_scores) > 0) {
    
    if (method == "quantile") {
      ord <- order(pos_scores)
      bin <- integer(length(pos_scores))
      
      base_size <- floor(length(pos_scores) / n_pos_bins)
      sizes <- rep(base_size, n_pos_bins)
      remainder <- length(pos_scores) - base_size * n_pos_bins
      if (remainder > 0) sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1
      
      start <- 1
      for (j in seq_len(n_pos_bins)) {
        end <- start + sizes[j] - 1
        idx_j <- ord[start:end]
        bin[idx_j] <- j
        start <- end + 1
      }
      
      lvl[pos_idx] <- bin
      thresholds <- sapply(seq_len(n_pos_bins), function(j)
        max(pos_scores[bin == j]))
      
    } else {  ## kmeans
      km <- kmeans(matrix(pos_scores), centers = n_pos_bins, nstart = kmeans_nstart)
      ord <- order(km$centers)
      map <- setNames(seq_along(ord), ord)
      lvl[pos_idx] <- map[as.character(km$cluster)]
      thresholds <- sort(as.numeric(km$centers))
    }
  }
  
  lvl <- as.integer(lvl)
  
  ## 4. Generate labels
  if (is.null(labels)) {
    labels <- c("Uninfected", paste0("Level_", 1:(n_levels - 1)))
  }
  lab_vec <- labels[pmax(1, lvl + 1)]
  
  ## 5. Write back to Seurat object
  obj[[level_col]] <- lvl
  obj[[label_col]] <- lab_vec
  
  ## 6. DimPlot visualization
  plots <- NULL
  if (isTRUE(plot)) {
    
    if (is.null(level_cols)) {
      p1 <- DimPlot(
        obj,
        group.by = level_col,
        repel = TRUE,
        label = TRUE
      )
    } else {
      p1 <- DimPlot(
        obj,
        group.by = level_col,
        repel = TRUE,
        label = TRUE,
        cols = level_cols
      )
    }
    
    p1 <- p1 +
      ggtitle("Infection level") +
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
    
    if (is.null(label_cols)) {
      p2 <- DimPlot(
        obj,
        group.by = label_col,
        repel = TRUE,
        label = TRUE
      )
    } else {
      p2 <- DimPlot(
        obj,
        group.by = label_col,
        repel = TRUE,
        label = TRUE,
        cols = label_cols
      )
    }
    
    p2 <- p2 +
      ggtitle("Infection label") +
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
    
    plots <- list(level_plot = p1, label_plot = p2)
  }
  
  list(
    thresholds = thresholds,
    table      = table(lvl),
    levels     = lvl,
    labels     = lab_vec,
    plots      = plots
  )
}
