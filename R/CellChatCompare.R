#' Compare Cell–Cell Communication Between Two Merged Condition Groups Using CellChat
#'
#' @description
#' `CellChatCompare` compares cell–cell communication between **two condition groups**
#' created by merging values of a metadata column (`compare_col`) into:
#' \itemize{
#'   \item Group A: `compare_1`
#'   \item Group B: `compare_2`
#' }
#'
#' The function:
#' \itemize{
#'   \item Optionally filters the Seurat object by `group_include`, `celltype_include`, and `level_include`
#'   \item Builds two merged groups (A/B) according to `compare_col` and (`compare_1`, `compare_2`)
#'   \item **Pre-checks cell numbers per cell type within each group**:
#'     \itemize{
#'       \item If a cell type has fewer than `min_cells_per_group` cells in that group, it is removed
#'     }
#'   \item Runs CellChat for each group independently
#'   \item Merges the two CellChat objects for comparison
#'   \item Saves and returns comparison plots:
#'     \itemize{
#'       \item Interaction count/weight barplots
#'       \item Global heatmaps (count & weight)
#'       \item Signaling-role heatmaps (A vs B)
#'       \item Circle plots for each group
#'     }
#' }
#'
#' @param obj A \code{Seurat} object.
#' @param compare_col Character. Metadata column used to define conditions to merge
#'   (e.g. \code{"timepoint"}, \code{"infection_level"}, \code{"treatment"}).
#'   If \code{NULL}, defaults to \code{group_col}.
#' @param compare_1 Vector. Values in \code{compare_col} to be merged into **Group A**.
#' @param compare_2 Vector. Values in \code{compare_col} to be merged into **Group B**.
#' @param group_col Character. Metadata column used for optional pre-filtering (default \code{"timepoint"}).
#' @param group_include Optional vector. Values of \code{group_col} to keep. Default \code{NULL}.
#' @param cell_type Character. Metadata column defining sender/receiver cell types (default \code{"manual_label"}).
#' @param celltype_include Optional vector. Cell types to keep before analysis. Default \code{NULL}.
#' @param level_col Character. Metadata column defining infection levels (default \code{"infection_level"}).
#' @param level_include Optional vector. Infection levels to keep before analysis. Default \code{NULL}.
#' @param assay Character. Seurat assay to use (default \code{"RNA"}).
#' @param slot Character. Expression slot used by CellChat (default \code{"data"}).
#' @param organism Character. \code{"human"} or \code{"mouse"} (selects CellChat database).
#' @param min_cells Integer. Minimum total cells required to run CellChat for a group (A or B).
#' @param min_cells_per_group Integer. Minimum cells per \code{cell_type} within each group.
#'   Cell types failing this threshold in a given group are dropped for that group.
#' @param save_dir Character. Directory to store outputs.
#' @param group_names Character length-2. Names for the two merged groups (default \code{c("GroupA","GroupB")}).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{cellchat_list}: list of CellChat objects for the two merged groups
#'   \item \code{merged}: merged CellChat object for comparison
#'   \item \code{plot_compare}: barplot (count & weight)
#'   \item \code{plot_heatmap}: global heatmap list (count & weight)
#'   \item \code{plot_heatmap_signal}: signaling-role heatmap (optional)
#'   \item \code{circle_plot_function}: function that draws circle plots for both groups
#'   \item \code{kept_cells}: named list of cell barcodes used for each group
#'   \item \code{dropped_celltypes}: named list of dropped cell types per group (by min_cells_per_group)
#' }
#'
#' @import Seurat
#' @import ggplot2
#' @importFrom graphics par
#'
#' @export
CellChatCompare <- function(
    obj,
    compare_col       = NULL,
    compare_1,
    compare_2,
    group_col         = "timepoint",
    group_include     = NULL,
    cell_type         = "manual_label",
    celltype_include  = NULL,
    level_col         = "infection_level",
    level_include     = NULL,
    assay             = "RNA",
    slot              = "data",
    organism          = c("human", "mouse"),
    min_cells         = 10,
    min_cells_per_group = 10,
    save_dir          = "CellChat_two_group_compare",
    group_names       = c("GroupA", "GroupB")
){
  # Check optional dependencies only when this function is called.
  optional_pkgs <- c("CellChat", "ComplexHeatmap")
  missing_optional <- optional_pkgs[
    !vapply(optional_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_optional) > 0L) {
    stop(
      "Optional package(s) required for CellChatCompare(): ",
      paste(missing_optional, collapse = ", "),
      ". Please install them before using this function.",
      call. = FALSE
    )
  }

  stopifnot(inherits(obj, "Seurat"))
  organism <- match.arg(organism)
  
  if (length(group_names) != 2) stop("'group_names' must be length 2.")
  gA_name <- group_names[[1]]
  gB_name <- group_names[[2]]
  
  md <- obj@meta.data
  
  # ---------- basic checks ----------
  if (!group_col %in% colnames(md)) stop("Column '", group_col, "' not found in meta.data.")
  if (!cell_type %in% colnames(md)) stop("Column '", cell_type, "' not found in meta.data.")
  if (!level_col %in% colnames(md)) stop("Column '", level_col, "' not found in meta.data.")
  if (is.null(compare_col)) compare_col <- group_col
  if (!compare_col %in% colnames(md)) stop("Column '", compare_col, "' not found in meta.data.")
  
  # ---------- 1) global filtering ----------
  keep <- rep(TRUE, nrow(md))
  if (!is.null(group_include)) {
    keep <- keep & as.character(md[[group_col]]) %in% as.character(group_include)
  }
  if (!is.null(celltype_include)) {
    keep <- keep & as.character(md[[cell_type]]) %in% as.character(celltype_include)
  }
  if (!is.null(level_include)) {
    keep <- keep & as.character(md[[level_col]]) %in% as.character(level_include)
  }
  cells_keep <- rownames(md)[keep]
  if (length(cells_keep) == 0) {
    stop("No cells left after filtering by group_include / celltype_include / level_include.")
  }
  
  obj_sub <- subset(obj, cells = cells_keep)
  md_sub  <- obj_sub@meta.data
  
  # ---------- 2) build merged groups A/B ----------
  compare_chr <- as.character(md_sub[[compare_col]])
  A_vals <- as.character(compare_1)
  B_vals <- as.character(compare_2)
  
  A_cells_all <- rownames(md_sub)[compare_chr %in% A_vals]
  B_cells_all <- rownames(md_sub)[compare_chr %in% B_vals]
  
  if (length(A_cells_all) == 0) stop("No cells matched compare_1 in column '", compare_col, "'.")
  if (length(B_cells_all) == 0) stop("No cells matched compare_2 in column '", compare_col, "'.")
  
  # ---------- helper: filter cell types by min_cells_per_group within one group ----------
  filter_celltypes_in_group <- function(cells_vec, group_label) {
    meta_g <- md_sub[cells_vec, , drop = FALSE]
    ct <- as.character(meta_g[[cell_type]])
    tab <- table(ct)
    keep_ct <- names(tab)[as.integer(tab) >= min_cells_per_group]
    drop_ct <- setdiff(names(tab), keep_ct)
    
    if (length(keep_ct) == 0) {
      return(list(cells = character(0), dropped = names(tab), kept_ct = character(0), tab = tab))
    }
    cells2 <- rownames(meta_g)[ct %in% keep_ct]
    list(cells = cells2, dropped = drop_ct, kept_ct = keep_ct, tab = tab)
  }
  
  A_f <- filter_celltypes_in_group(A_cells_all, gA_name)
  B_f <- filter_celltypes_in_group(B_cells_all, gB_name)
  
  # If either group loses all cells after per-celltype filtering, stop with a clear message
  if (length(A_f$cells) < min_cells) {
    stop(gA_name, " has only ", length(A_f$cells),
         " cells after min_cells_per_group filtering (< min_cells=", min_cells, ").")
  }
  if (length(B_f$cells) < min_cells) {
    stop(gB_name, " has only ", length(B_f$cells),
         " cells after min_cells_per_group filtering (< min_cells=", min_cells, ").")
  }
  
  # subset to each group
  data_A <- subset(obj_sub, cells = A_f$cells)
  data_B <- subset(obj_sub, cells = B_f$cells)
  
  # choose CellChat DB
  if (organism == "human") {
    CellChatDB <- CellChat::CellChatDB.human
  } else {
    CellChatDB <- CellChat::CellChatDB.mouse
  }
  
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  # ---------- helper: run CellChat for one Seurat subset ----------
  run_one_cellchat <- function(data_cond, cond_name) {
    meta_cond <- data_cond@meta.data
    meta_cond[[cell_type]] <- droplevels(factor(meta_cond[[cell_type]]))
    data_cond@meta.data <- meta_cond
    Seurat::Idents(data_cond) <- meta_cond[[cell_type]]
    
    Seurat::DefaultAssay(data_cond) <- assay
    data.input <- Seurat::GetAssayData(data_cond, assay = assay, slot = slot)
    if (ncol(data.input) == 0) stop("Expression matrix is empty for ", cond_name, ".")
    
    cellchat_cond <- CellChat::createCellChat(object = data_cond, group.by = cell_type)
    cellchat_cond@DB <- CellChatDB
    
    cellchat_cond <- CellChat::subsetData(cellchat_cond)
    cellchat_cond <- CellChat::identifyOverExpressedGenes(cellchat_cond)
    cellchat_cond <- CellChat::identifyOverExpressedInteractions(cellchat_cond)
    cellchat_cond <- CellChat::computeCommunProb(cellchat_cond, raw.use = TRUE)
    cellchat_cond <- CellChat::filterCommunication(cellchat_cond, min.cells = min_cells)
    cellchat_cond <- CellChat::computeCommunProbPathway(cellchat_cond)
    cellchat_cond <- CellChat::aggregateNet(cellchat_cond)
    cellchat_cond <- CellChat::netAnalysis_computeCentrality(cellchat_cond)
    
    saveRDS(cellchat_cond, file = file.path(save_dir, paste0("CellChat_", compare_col, "_", cond_name, ".rds")))
    cellchat_cond
  }
  
  message("\n==============================")
  message(">> Building CellChat for merged ", compare_col, ": ", gA_name,
          " (", paste(A_vals, collapse = ","), ")")
  message("==============================\n")
  cellchat_A <- run_one_cellchat(data_A, gA_name)
  
  message("\n==============================")
  message(">> Building CellChat for merged ", compare_col, ": ", gB_name,
          " (", paste(B_vals, collapse = ","), ")")
  message("==============================\n")
  cellchat_B <- run_one_cellchat(data_B, gB_name)
  
  cellchat_list <- list()
  cellchat_list[[gA_name]] <- cellchat_A
  cellchat_list[[gB_name]] <- cellchat_B
  
  # ---------- merge & visualize ----------
  object.list <- cellchat_list
  cellchat_merged <- CellChat::mergeCellChat(object.list, add.names = names(object.list))
  
  # compareInteractions barplots: count & weight
  gg1 <- CellChat::compareInteractions(cellchat_merged, show.legend = FALSE, group = 1:2)
  gg2 <- CellChat::compareInteractions(cellchat_merged, show.legend = FALSE, group = 1:2, measure = "weight")
  p_compare <- gg1 + gg2
  
  # global heatmaps
  h1 <- CellChat::netVisual_heatmap(cellchat_merged)
  h2 <- CellChat::netVisual_heatmap(cellchat_merged, measure = "weight")
  p_heatmap <- h1 + h2
  
  # signaling-role heatmaps (A vs B)
  p_heatmap_signal <- NULL
  pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
  if (length(pathway.union) > 0) {
    ht1 <- CellChat::netAnalysis_signalingRole_heatmap(
      object.list[[1]],
      pattern       = "all",
      signaling     = pathway.union,
      title         = names(object.list)[1],
      color.heatmap = "OrRd"
    )
    ht2 <- CellChat::netAnalysis_signalingRole_heatmap(
      object.list[[2]],
      pattern       = "all",
      signaling     = pathway.union,
      title         = names(object.list)[2],
      color.heatmap = "OrRd"
    )
    p_heatmap_signal <- ht1 + ht2
  }
  
  # circle plots
  circle_plot_fun <- function() {
    groupSize_list <- lapply(object.list, function(x) as.numeric(table(x@idents)))
    
    weight.max.count <- CellChat::getMaxWeight(object.list, attribute = c("idents", "count"))
    par(mfrow = c(1, length(object.list)), xpd = TRUE)
    for (k in seq_along(object.list)) {
      CellChat::netVisual_circle(
        object.list[[k]]@net$count,
        vertex.weight   = groupSize_list[[k]],
        weight.scale    = TRUE,
        label.edge      = FALSE,
        edge.weight.max = weight.max.count[2],
        title.name      = paste0("Number of interactions - ", names(object.list)[k])
      )
    }
    
    weight.max.str <- CellChat::getMaxWeight(object.list, attribute = c("idents", "weight"))
    par(mfrow = c(1, length(object.list)), xpd = TRUE)
    for (k in seq_along(object.list)) {
      CellChat::netVisual_circle(
        object.list[[k]]@net$weight,
        vertex.weight   = groupSize_list[[k]],
        weight.scale    = TRUE,
        label.edge      = FALSE,
        edge.weight.max = weight.max.str[2],
        title.name      = paste0("Interaction weights - ", names(object.list)[k])
      )
    }
  }
  
  # ---------- save plots ----------
  ggplot2::ggsave(
    file.path(save_dir, paste0("CellChat_compare_barplot_", compare_col, "_", gA_name, "_vs_", gB_name, ".pdf")),
    p_compare, width = 10, height = 5
  )
  
  pdf(file.path(save_dir, paste0("CellChat_compare_heatmap_", compare_col, "_", gA_name, "_vs_", gB_name, ".pdf")),
      width = 10, height = 5)
  ComplexHeatmap::draw(p_heatmap, heatmap_legend_side = "right")
  dev.off()
  
  if (!is.null(p_heatmap_signal)) {
    pdf(file.path(save_dir, paste0("CellChat_compare_signaling_heatmap_", compare_col, "_", gA_name, "_vs_", gB_name, ".pdf")),
        width = 10, height = 5)
    ComplexHeatmap::draw(p_heatmap_signal, heatmap_legend_side = "right")
    dev.off()
  }
  
  pdf(file.path(save_dir, paste0("CellChat_compare_circle_", compare_col, "_", gA_name, "_vs_", gB_name, ".pdf")),
      width = 8, height = 4)
  circle_plot_fun()
  dev.off()
  
  invisible(list(
    cellchat_list        = cellchat_list,
    merged               = cellchat_merged,
    plot_compare         = p_compare,
    plot_heatmap         = p_heatmap,
    plot_heatmap_signal  = p_heatmap_signal,
    circle_plot_function = circle_plot_fun,
    kept_cells = list(
      GroupA = A_f$cells,
      GroupB = B_f$cells
    ),
    dropped_celltypes = list(
      GroupA = A_f$dropped,
      GroupB = B_f$dropped
    )
  ))
}