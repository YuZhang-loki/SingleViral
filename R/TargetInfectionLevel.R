#' Discretize viral-load into infection levels and write to Seurat meta.data
#'
#' @description
#' Assigns per-cell infection levels using a viral-load column in a Seurat object.
#' Cells with (absolute) score <= `zero_tol` are forced to level 0 ("Uninfected").
#' Positive cells are partitioned into `n_levels - 1` bins using either quantiles
#' or 1D k-means on positive scores only. Both the numeric level and a label are
#' written back to `meta.data`.
#'
#' @param obj A **Seurat** object.
#' @param score_col Character. Column in `obj@meta.data` that stores viral-load scores.
#' @param n_levels Integer (>=2). Total number of levels **including** level 0.
#'   Level 0 is reserved for uninfected cells.
#' @param method `"quantile"` or `"kmeans"`. Binning strategy for **positive** cells.
#' @param label_col Character. Name of the label column to write (e.g., `"infection_label"`).
#' @param level_col Character. Name of the numeric level column to write (e.g., `"infection_level"`).
#' @param labels Optional character vector of length `n_levels`. Custom labels where
#'   `labels[1]` corresponds to level 0. If `NULL`, defaults to
#'   `c("Uninfected", paste0("Level_", 1:(n_levels - 1)))`.
#' @param zero_tol Numeric. Scores with `abs(score) <= zero_tol` are treated as zero.
#' @param kmeans_nstart Integer. `nstart` for k-means when `method = "kmeans"`.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{object} Updated Seurat object.
#'   \item \code{thresholds} Numeric vector of bin edges/centers used for positives.
#'   \item \code{table} Frequency table of assigned levels.
#'   \item \code{levels} Integer vector of per-cell levels (aligned with columns of \code{obj}).
#'   \item \code{labels} Character vector of per-cell level labels.
#' }
#'
#' @examples
#' \dontrun{
#' res <- TargetInfectionLevel(
#'   obj       = data,
#'   score_col = "viral_load",
#'   n_levels  = 5,
#'   method    = "quantile",
#'   labels    = c("Uninfected", "Low infection", "Moderate infection",
#'                 "Moderate to high", "High infection")
#' )
#' table(res$object$infection_level)
#' }
#'
#' @export
TargetInfectionLevel <- function(obj,
                                 score_col = "viral_load",
                                 n_levels  = 3,
                                 method    = c("quantile", "kmeans"),
                                 label_col = "infection_label",
                                 level_col = "infection_level",
                                 labels    = NULL,
                                 zero_tol  = 1e-12,
                                 kmeans_nstart = 10) {
  stopifnot(inherits(obj, "Seurat"))
  method <- match.arg(method)
  if (n_levels < 2) stop("`n_levels` must be >= 2.")

  # --- check score column in meta.data ---
  if (!(score_col %in% colnames(obj[[]]))) {
    stop(sprintf("Column '%s' not found in meta.data. Please write viral-load first.", score_col))
  }

  # 1) fetch scores
  scores <- obj[[score_col]][, 1]
  names(scores) <- colnames(obj)

  # 2) zero (uninfected) handling
  zero_mask <- abs(scores) <= zero_tol
  lvl <- rep(NA_integer_, length(scores))
  names(lvl) <- names(scores)
  lvl[zero_mask] <- 0

  pos_scores <- scores[!zero_mask]
  pos_idx    <- names(pos_scores)

  if (length(pos_scores) == 0) {
    # all cells are zero/uninfected
    lvl[is.na(lvl)] <- 0
    obj[[level_col]] <- lvl
    if (is.null(labels)) {
      labels <- c("Uninfected", paste0("Level_", 1:(n_levels - 1)))
    } else if (length(labels) != n_levels) {
      stop("`labels` must have length equal to `n_levels`.")
    }
    lab_vec <- labels[pmax(1, lvl + 1)]
    obj[[label_col]] <- lab_vec
    return(list(object = obj, thresholds = 0, table = table(lvl),
                levels = lvl, labels = lab_vec))
  }

  # 3) bin positives into 1..(n_levels-1)
  n_pos_bins <- n_levels - 1

  if (method == "quantile") {
    probs <- seq(0, 1, length.out = n_pos_bins + 1)
    brks  <- unique(quantile(pos_scores, probs = probs, names = FALSE, type = 7))
    # if too few unique breaks, fallback to equal-width range
    if (length(brks) <= 2) {
      brks <- seq(min(pos_scores), max(pos_scores), length.out = n_pos_bins + 1)
    }
    bins <- cut(pos_scores, breaks = brks, include.lowest = TRUE, labels = FALSE)
    lvl[pos_idx] <- bins
    thresholds <- brks[-1]
  } else {
    km  <- kmeans(matrix(pos_scores), centers = n_pos_bins, nstart = kmeans_nstart)
    ord <- order(km$centers)                   # low -> high
    map <- setNames(seq_along(ord), ord)       # cluster -> level index (1..K)
    lvl[pos_idx] <- map[as.character(km$cluster)]
    centers <- sort(as.numeric(km$centers))
    # midpoints between adjacent centers can be treated as thresholds
    thresholds <- c(0, (head(centers, -1) + tail(centers, -1)) / 2)
  }

  lvl <- as.integer(lvl)

  # 4) labels
  if (is.null(labels)) {
    labels <- c("Uninfected", paste0("Level_", 1:(n_levels - 1)))
  } else if (length(labels) != n_levels) {
    stop("`labels` must have length equal to `n_levels`.")
  }
  lab_vec <- labels[pmax(1, lvl + 1)]  # level 0 -> labels[1]

  # 5) write back
  obj[[level_col]] <- lvl
  obj[[label_col]] <- lab_vec

  list(
    object     = obj,
    thresholds = thresholds,
    table      = table(lvl),
    levels     = lvl,
    labels     = lab_vec
  )
}
