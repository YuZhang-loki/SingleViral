#' CellChat analysis by infection status (derived from viral load), with overview
#' circle plots, heatmaps and per-pathway circle plots
#'
#' @description
#' Builds an infection-status label per cell using a viral-load column and a
#' threshold, appends it as \code{<cell_type>_<inf/uninf>} to a new grouping column,
#' runs CellChat on this grouping, exports the communication table, and saves
#' overview circle plots (counts/weights). Optionally draws (1) two global
#' heatmaps for interaction count and interaction strength, and (2) per-pathway
#' circle plots (top-K by total weight) in a combined PDF and individual PDFs.
#'
#' @param obj A \strong{Seurat} object.
#' @param cell_type Character. Meta-data column for cell type labels (e.g. \code{"manual_label"}).
#' @param score_col Character. Meta-data column with viral-load scores (numeric).
#' @param group_col Character. Meta-data column for groups.
#' @param level_col Character. Infection level column (only used when \code{level_include} is not \code{NULL}).
#' @param level_include Optional vector of infection levels to keep (e.g., \code{c(0,1,2,3,4)}). \code{NULL} keeps all.
#'
#' @param infection_threshold Numeric. Cells with \code{viral_load > infection_threshold}
#'   are labeled "inf", otherwise "uninf".
#' @param new_col Character. Name of the new grouping column to write (default \code{"cell_infection"}).
#' @param suffix_sep Character separator between \code{cell_type} and infection suffix.
#'
#' @param celltype_include Optional character vector of cell types to retain.
#' @param group_include Optional character vector of groups to retain.
#'
#' @param organism \code{"human"} or \code{"mouse"}; choose CellChat DB accordingly.
#' @param db_category Optional character to subset the CellChat DB (\code{CellChat::subsetDB}).
#' @param assay Optional assay for the Seurat object. If \code{NULL}, keep default.
#' @param min_cells Integer. Minimum cells per group for \code{CellChat::filterCommunication}.
#' @param raw_use Logical. Passed to \code{CellChat::computeCommunProb}.
#'
#' @param save_dir Output directory for results (tables/plots).
#' @param save_plots Logical. Whether to save the circle/heatmap plots.
#' @param width,height Numeric. Size (inches) for overview circle plots PDF.
#' @param width_merge,height_merge Numeric. Size for the per-pathways combined PDF.
#' @param label_size,label_size_merge Numeric. Node label cex for overview / merged plots.
#'
#' @param mfrow Integer vector of length 2 for combined per-pathway layout (rows, cols).
#' @param per_pathway_plots Logical. Whether to draw per-pathway circle plots.
#' @param per_pathway_top Integer. Number of top pathways (by total weight) to plot.
#' @param per_pathway_folder Character. Subfolder name for individual per-pathway PDFs.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{object_with_label}: Seurat object with infection columns written
#'   \item \code{cellchat}: \code{CellChat} object after computation
#'   \item \code{interactions}: data.frame of communication results
#'   \item \code{outdir}: output directory used
#' }
#'
#' @examples
#' \dontrun{
#' res <- CellChatByCellInfection(
#'   obj         = data,
#'   cell_type   = "manual_label",
#'   score_col   = "viral_load",
#'   group_col = "timepoint",
#'   infection_threshold = 0,
#'   save_dir    = "cellchat_cell_infection",
#'   per_pathway_plots = TRUE,
#'   per_pathway_top   = 12
#' )
#' }
#'
#' @seealso \code{\link[CellChat]{createCellChat}}, \code{\link[CellChat]{subsetDB}},
#'   \code{\link[CellChat]{computeCommunProb}}, \code{\link[CellChat]{netVisual_circle}}
#'
#' @importFrom utils write.csv
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics par
#' @export
CellChatByCellInfection <- function(
    obj,
    # — Meta-data columns — #
    cell_type     = "manual_label",
    score_col     = "viral_load",
    group_col     = "timepoint",
    level_col     = "infection_level",
    level_include    = NULL,

    # — Infection labeling — #
    infection_threshold = 0,
    new_col       = "cell_infection",
    suffix_sep    = "_",

    # — Optional filters — #
    celltype_include    = NULL,
    group_include   = NULL,

    # — CellChat config — #
    organism      = c("human","mouse"),
    db_category   = NULL,
    assay         = NULL,
    min_cells     = 5,
    raw_use       = TRUE,

    # — Output & plots — #
    save_dir      = "cellchat_cell_infection",
    save_plots    = TRUE,
    width         = 8, height = 6,
    width_merge   = 16, height_merge = 12,
    label_size    = 1,  label_size_merge = 1,

    # — Per-pathway plotting — #
    mfrow = c(4, 3),
    per_pathway_plots  = TRUE,
    per_pathway_top    = 12,
    per_pathway_folder = "pathways"
){
  # Check optional dependencies only when this function is called.
  optional_pkgs <- c("CellChat")
  missing_optional <- optional_pkgs[
    !vapply(optional_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_optional) > 0L) {
    stop(
      "Optional package(s) required for CellChatByCellInfection(): ",
      paste(missing_optional, collapse = ", "),
      ". Please install them before using this function.",
      call. = FALSE
    )
  }

  if (!inherits(obj, "Seurat")) stop("'obj' must be a Seurat object.")

  md <- obj@meta.data
  for (cn in c(cell_type, score_col, group_col)) {
    if (!cn %in% colnames(md)) stop(sprintf("'%s' not found in meta.data.", cn))
  }

  # set assay if requested
  if (!is.null(assay)) Seurat::DefaultAssay(obj) <- assay

  # — filter by timepoints and labels —
  keep <- rep(TRUE, nrow(md))
  if (!is.null(celltype_include))  keep <- keep & (md[[cell_type]] %in% celltype_include)
  if (!is.null(group_include)) keep <- keep & (as.character(md[[group_col]]) %in% as.character(group_include))
  if (!all(keep)) obj <- subset(obj, cells = rownames(md)[keep])

  # — optional filter by infection levels —
  if (!is.null(level_include)) {
    md2 <- obj@meta.data
    if (!level_col %in% colnames(md2)) stop(sprintf("'%s' not found in meta.data.", level_col))
    keep_lvl <- md2[[level_col]] %in% level_include
    obj <- subset(obj, cells = rownames(md2)[keep_lvl])
  }

  # — build infection status from viral load —
  md <- obj@meta.data
  if (!is.numeric(md[[score_col]])) stop(sprintf("'%s' must be numeric.", score_col))
  infection_status <- ifelse(md[[score_col]] > infection_threshold, "inf", "uninf")

  # write to meta.data
  obj@meta.data$infection_threshold <- infection_threshold
  obj@meta.data$infection_status    <- infection_status
  obj@meta.data[[new_col]]          <- paste0(md[[cell_type]], suffix_sep, infection_status)

  # — choose CellChat DB —
  organism  <- match.arg(organism)
  CellChatDB <- if (identical(organism, "human")) CellChat::CellChatDB.human else CellChat::CellChatDB.mouse
  if (!is.null(db_category)) {
    CellChatDB <- CellChat::subsetDB(CellChatDB, search = db_category)
  }

  # — create & run CellChat —
  cellchat <- CellChat::createCellChat(object = obj, group.by = new_col)
  cellchat@DB <- CellChatDB
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, raw.use = raw_use)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)

  # — output dir & export table —
  outdir <- file.path(save_dir, paste0("by_", new_col))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  df.net <- CellChat::subsetCommunication(cellchat)
  utils::write.csv(df.net, file = file.path(outdir, "cellcell_interactions.csv"), row.names = FALSE)

  # — overview circle plots & heatmaps —
  groupSize <- as.numeric(table(cellchat@idents))
  if (save_plots) {
    ## circle: number of interactions
    grDevices::pdf(file.path(outdir, "circle_counts.pdf"), width = width, height = height)
    CellChat::netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                               weight.scale = TRUE, label.edge = FALSE,
                               title.name = "Number of interactions",
                               vertex.label.cex = label_size)
    grDevices::dev.off()

    ## circle: interaction weights/strength
    grDevices::pdf(file.path(outdir, "circle_weights.pdf"), width = width, height = height)
    CellChat::netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                               weight.scale = TRUE, label.edge = FALSE,
                               title.name = "Interaction weights/strength",
                               vertex.label.cex = label_size)
    grDevices::dev.off()

    ## heatmap of interaction counts
    if (is.matrix(cellchat@net$count) && nrow(cellchat@net$count) > 0 &&
        ncol(cellchat@net$count) > 0) {
      grDevices::pdf(file.path(outdir, "heatmap_counts.pdf"),
                     width = width, height = height)
      pheatmap::pheatmap(
        cellchat@net$count,
        border_color   = "black",
        cluster_cols   = FALSE,
        cluster_rows   = FALSE,
        fontsize       = 10,
        display_numbers = TRUE,
        number_color   = "black",
        number_format  = "%.0f",
        main           = "Number of interactions"
      )
      grDevices::dev.off()
    }

    ##  heatmap of interaction weights
    if (is.matrix(cellchat@net$weight) && nrow(cellchat@net$weight) > 0 &&
        ncol(cellchat@net$weight) > 0) {
      grDevices::pdf(file.path(outdir, "heatmap_weights.pdf"),
                     width = width, height = height)
      pheatmap::pheatmap(
        cellchat@net$weight,
        border_color   = "black",
        cluster_cols   = FALSE,
        cluster_rows   = FALSE,
        fontsize       = 10,
        display_numbers = TRUE,
        number_color   = "black",
        main           = "Interaction strength"
      )
      grDevices::dev.off()
    }
  }

  # — per-pathway circle plots —
  if (save_plots && per_pathway_plots) {
    mat <- cellchat@net$weight
    if (is.matrix(mat) && nrow(mat) > 0 && ncol(mat) > 0) {
      ord  <- order(rowSums(mat, na.rm = TRUE), decreasing = TRUE)
      pick <- ord[seq_len(min(per_pathway_top, length(ord)))]

      grDevices::pdf(file.path(outdir, sprintf("per_pathway_top%d.pdf", length(pick))),
                     width = width_merge, height = height_merge)
      graphics::par(mfrow = mfrow, xpd = TRUE, mar = c(1,1,1,1))

      for (k in pick) {
        if (sum(mat[k, ], na.rm = TRUE) <= 0) next
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[k, ] <- mat[k, ]
        CellChat::netVisual_circle(
          mat2,
          vertex.weight    = groupSize,
          weight.scale     = TRUE,
          arrow.width      = 1,
          arrow.size       = 0.05,
          edge.weight.max  = max(mat, na.rm = TRUE),
          title.name       = rownames(mat)[k],
          vertex.label.cex = label_size_merge
        )
      }
      grDevices::dev.off()

      pw_dir <- file.path(outdir, per_pathway_folder)
      dir.create(pw_dir, recursive = TRUE, showWarnings = FALSE)
      for (k in pick) {
        if (sum(mat[k, ], na.rm = TRUE) <= 0) next
        fn <- paste0(gsub("[^A-Za-z0-9._-]+", "_", rownames(mat)[k]), ".pdf")
        grDevices::pdf(file.path(pw_dir, fn), width = width, height = height)
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[k, ] <- mat[k, ]
        CellChat::netVisual_circle(
          mat2,
          vertex.weight    = groupSize,
          weight.scale     = TRUE,
          arrow.width      = 1,
          arrow.size       = 0.05,
          edge.weight.max  = max(mat, na.rm = TRUE),
          title.name       = rownames(mat)[k],
          vertex.label.cex = label_size
        )
        grDevices::dev.off()
      }
    } else {
      message("[Skip] per-pathway plots: weight matrix not available.")
    }
  }

  saveRDS(cellchat, file.path(outdir, "cellchat_cell_infection.rds"))
  message("Done. Outputs saved in: ", normalizePath(outdir, mustWork = FALSE))

  list(
    object_with_label = obj,
    cellchat = cellchat,
    interactions = df.net,
    outdir = outdir
  )
}
