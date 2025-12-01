#' Compare Cell–Cell Communication Across Multiple Conditions Using CellChat
#'
#' @description
#' `CellChatCompare` performs multi-condition comparison of cell–cell
#' communication. It:
#' \itemize{
#'   \item Filters Seurat object (timepoint / cell type / infection level)
#'   \item Splits cells by a metadata column (`meta_col`)
#'   \item Runs CellChat independently for each condition
#'   \item Merges CellChat objects for comparison
#'   \item Generates publication-ready visualizations:
#'     \itemize{
#'       \item Interaction count / strength barplots
#'       \item Global CellChat heatmaps
#'       \item Signaling-role heatmaps (per-pair comparison)
#'       \item Circular network plots for each condition
#'     }
#' }
#'
#' @param obj A Seurat object.
#' @param meta_col Metadata column used for splitting conditions (e.g. `"timepoint"`,
#'   `"infection_level"`, `"treatment"`). Default: `time_col`.
#' @param meta_include Character vector of values of `meta_col` to include.
#'   If `NULL`, all values in `meta_col` are used.
#' @param time_col Metadata column defining timepoint (default `"timepoint"`).
#' @param timepoints Timepoints to keep. Default `NULL`.
#' @param group_col Metadata column defining sender/receiver cell type
#'   (default `"manual_label"`).
#' @param group_include Vector of cell types to keep. Default `NULL`.
#' @param level_col Metadata column defining infection levels (default `"infection_level"`).
#' @param level_include Infection levels to keep. Default `NULL`.
#' @param assay Seurat assay to use (default `"RNA"`).
#' @param slot Expression slot (e.g. `"data"`, `"counts"`). Default `"data"`.
#' @param organism `"human"` or `"mouse"` (for CellChat database).
#' @param min_cells Minimum cells required to run CellChat for a condition.
#' @param workers Number of parallel workers.
#' @param save_dir Directory to store output plots.
#' @param vi Integer specifying which two conditions to compare in signaling-role
#'   heatmap (`vi` vs `vi + 1`).
#'
#' @return A list including:
#' \itemize{
#'   \item `cellchat_list` — list of CellChat objects per condition
#'   \item `merged` — merged CellChat object
#'   \item `plot_compare` — barplot object
#'   \item `plot_heatmap` — CellChat heatmap list
#'   \item `plot_heatmap_signal` — signaling-role heatmap (optional)
#'   \item `circle_plot_function` — function used to generate circle plots
#' }
#'
#' @details
#' This function automates multi-condition communication analysis with CellChat.
#' It supports time-course, infection-level, treatment-group, and cluster-level
#' comparisons.
#'
#' The function generates:
#' \itemize{
#'   \item Interaction-count comparison barplots
#'   \item Interaction-strength barplots
#'   \item Global communication heatmaps
#'   \item Signaling role heatmaps (per selected pair)
#'   \item Circular communication network plots (per condition, NOT merged)
#' }
#'
#' @examples
#' \dontrun{
#' res <- CellChatCompare(
#'   obj          = seurat_object,
#'   meta_col     = "infection_status",
#'   meta_include = c("uninfected", "infected"),
#'   group_col    = "manual_label",
#'   vi           = 1
#' )
#'
#' # draw circle plots
#' res$circle_plot_function()
#' }
#'
#' @import Seurat
#' @import CellChat
#' @import ComplexHeatmap
#' @import ggplot2
#' @importFrom graphics par
#' @importFrom utils write.csv
#'
#' @export
CellChatCompare <- function(
    obj,
    meta_col      = NULL,        # e.g. "timepoint", "infection_level", "group_id" ...
    meta_include  = NULL,        # which values of meta_col to compare
    time_col      = "timepoint",
    timepoints    = NULL,        # pre-filter by time
    group_col     = "manual_label",
    group_include = NULL,        # pre-filter by cell types
    level_col     = "infection_level",
    level_include = NULL,        # pre-filter by infection levels
    assay         = "RNA",
    slot          = "data",
    organism      = c("human", "mouse"),
    min_cells     = 10,
    workers       = 4,
    save_dir      = "CellChat_timepoint_compare",
    vi            = 1            # which pair in object.list to use for signalingRole heatmap
    ){
      stopifnot(inherits(obj, "Seurat"))
      organism <- match.arg(organism)
      
      md <- obj@meta.data
      
      ## ---------- 0. basic column checks ----------
      if (!time_col %in% colnames(md))
        stop("Column '", time_col, "' not found in meta.data.")
      if (!group_col %in% colnames(md))
        stop("Column '", group_col, "' not found in meta.data.")
      if (!level_col %in% colnames(md))
        stop("Column '", level_col, "' not found in meta.data.")
      
      ## ---------- 1. global filtering by time / group / infection level ----------
      keep <- rep(TRUE, nrow(md))
      if (!is.null(timepoints)) {
        keep <- keep & as.character(md[[time_col]]) %in% as.character(timepoints)
      }
      if (!is.null(group_include)) {
        keep <- keep & as.character(md[[group_col]]) %in% as.character(group_include)
      }
      if (!is.null(level_include)) {
        keep <- keep & as.character(md[[level_col]]) %in% as.character(level_include)
      }
      
      cells_keep <- rownames(md)[keep]
      if (length(cells_keep) == 0)
        stop("No cells left after filtering by timepoints / group_include / level_include.")
      
      obj_sub <- subset(obj, cells = cells_keep)
      md_sub  <- obj_sub@meta.data
      
      ## ---------- 2. split by meta_col / meta_include & run CellChat ----------
      if (is.null(meta_col)) {
        meta_col <- time_col
      }
      if (!meta_col %in% colnames(md_sub))
        stop("Column '", meta_col, "' not found in filtered meta.data.")
      
      if (is.null(meta_include)) {
        meta_include <- sort(unique(as.character(md_sub[[meta_col]])))
      } else {
        meta_include <- as.character(meta_include)
      }
      
      # choose CellChat DB
      if (organism == "human") {
        CellChatDB <- CellChatDB.human
      } else {
        CellChatDB <- CellChatDB.mouse
      }
      
      if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
      
      # optional: enable parallel
      # future::plan(future::multisession, workers = workers)
      
      cellchat_list <- list()
      
      for (cond in meta_include) {
        message("\n==============================")
        message(">> Building CellChat for ", meta_col, " = ", cond)
        message("==============================\n")
        
        cond_cells <- rownames(md_sub)[as.character(md_sub[[meta_col]]) == cond]
        if (length(cond_cells) < min_cells) {
          warning("Condition ", cond, " has only ", length(cond_cells),
                  " cells (< min_cells). Skipped.")
          next
        }
        
        data_cond <- subset(obj_sub, cells = cond_cells)
        meta_cond <- data_cond@meta.data
        
        # define sender/receiver groups
        meta_cond[[group_col]] <- droplevels(factor(meta_cond[[group_col]]))
        data_cond@meta.data <- meta_cond
        Seurat::Idents(data_cond) <- meta_cond[[group_col]]
        
        # expression matrix
        Seurat::DefaultAssay(data_cond) <- assay
        data.input <- Seurat::GetAssayData(data_cond, assay = assay, slot = slot)
        if (ncol(data.input) == 0) {
          warning("Expression matrix is empty for condition ", cond, ". Skipped.")
          next
        }
        
        # create CellChat object
        cellchat_cond <- CellChat::createCellChat(
          object   = data_cond,
          group.by = group_col
        )
        cellchat_cond@DB <- CellChatDB
        
        # standard CellChat workflow
        cellchat_cond <- CellChat::subsetData(cellchat_cond)
        cellchat_cond <- CellChat::identifyOverExpressedGenes(cellchat_cond)
        cellchat_cond <- CellChat::identifyOverExpressedInteractions(cellchat_cond)
        cellchat_cond <- CellChat::computeCommunProb(cellchat_cond, raw.use = TRUE)
        cellchat_cond <- CellChat::filterCommunication(cellchat_cond, min.cells = min_cells)
        cellchat_cond <- CellChat::computeCommunProbPathway(cellchat_cond)
        cellchat_cond <- CellChat::aggregateNet(cellchat_cond)
        # centrality (for signalingRole heatmap)
        cellchat_cond <- CellChat::netAnalysis_computeCentrality(cellchat_cond)
        
        cellchat_list[[cond]] <- cellchat_cond
        
        saveRDS(
          cellchat_cond,
          file = file.path(save_dir, paste0("CellChat_", meta_col, "_", cond, ".rds"))
        )
      }
      
      ## ---------- 3. require at least two conditions ----------
      cellchat_list <- cellchat_list[!vapply(cellchat_list, is.null, logical(1))]
      if (length(cellchat_list) < 2) {
        stop("Less than two valid conditions left in meta_col for comparison.")
      }
      
      ## ---------- 4. merge & visualize ----------
      object.list     <- cellchat_list
      cellchat_merged <- CellChat::mergeCellChat(
        object.list,
        add.names = names(object.list)
      )
      
      
      # 4.2 compareInteractions barplots: interaction count & weight
      gg1 <- CellChat::compareInteractions(
        cellchat_merged,
        show.legend = FALSE,
        group       = 1:length(object.list)
      )
      gg2 <- CellChat::compareInteractions(
        cellchat_merged,
        show.legend = FALSE,
        group       = 1:length(object.list),
        measure     = "weight"
      )
      p_compare <- gg1 + gg2
      
      # 4.3 global heatmaps (count & weight)
      h1 <- CellChat::netVisual_heatmap(cellchat_merged)
      h2 <- CellChat::netVisual_heatmap(cellchat_merged, measure = "weight")
      p_heatmap <- h1 + h2   # HeatmapList object
      
      # 4.4 signaling-role heatmaps，使用 object.list[[vi]] 和 [[vi+1]]
      p_heatmap_signal <- NULL
      if (vi < length(object.list)) {
        i <- vi
        pathway.union <- union(object.list[[i]]@netP$pathways,
                               object.list[[i + 1]]@netP$pathways)
        ht1 <- CellChat::netAnalysis_signalingRole_heatmap(
          object.list[[i]],
          pattern       = "all",
          signaling     = pathway.union,
          title         = names(object.list)[i],
          color.heatmap = "OrRd"
        )
        ht2 <- CellChat::netAnalysis_signalingRole_heatmap(
          object.list[[i + 1]],
          pattern       = "all",
          signaling     = pathway.union,
          title         = names(object.list)[i + 1],
          color.heatmap = "OrRd"
        )
        p_heatmap_signal <- ht1 + ht2
      } else {
        warning("`vi` is too large; cannot form a pair (vi, vi+1) for signalingRole heatmap.")
      }
      
      # 4.5 circle plots for each condition
      circle_plot_fun <- function() {
        # node size: groupSize per condition
        groupSize_list <- lapply(object.list, function(x) as.numeric(table(x@idents)))
        
        ## --- number of interactions ---
        weight.max.count <- CellChat::getMaxWeight(
          object.list, attribute = c("idents", "count")
        )
        par(mfrow = c(1, length(object.list)), xpd = TRUE)
        for (k in seq_along(object.list)) {
          CellChat::netVisual_circle(
            object.list[[k]]@net$count,
            vertex.weight  = groupSize_list[[k]],
            weight.scale   = TRUE,
            label.edge     = FALSE,
            edge.weight.max = weight.max.count[2],
            title.name     = paste0("Number of interactions - ", names(object.list)[k])
          )
        }
        
        ## --- interaction strength (weight) ---
        weight.max.str <- CellChat::getMaxWeight(
          object.list, attribute = c("idents", "weight")
        )
        par(mfrow = c(1, length(object.list)), xpd = TRUE)
        for (k in seq_along(object.list)) {
          CellChat::netVisual_circle(
            object.list[[k]]@net$weight,
            vertex.weight  = groupSize_list[[k]],
            weight.scale   = TRUE,
            label.edge     = FALSE,
            edge.weight.max = weight.max.str[2],
            title.name     = paste0("Interaction weights - ", names(object.list)[k])
          )
        }
      }
      
      ## ---------- 5. save all plots ----------
      # barplot
      ggplot2::ggsave(
        file.path(save_dir, paste0("CellChat_compare_barplot_", meta_col, ".pdf")),
        p_compare, width = 10, height = 5
      )
      
      # global heatmaps
      pdf(file.path(save_dir, paste0("CellChat_compare_heatmap_", meta_col, ".pdf")),
          width = 10, height = 5)
      ComplexHeatmap::draw(p_heatmap, heatmap_legend_side = "right")
      dev.off()
      
      # signaling-role heatmaps (if available)
      if (!is.null(p_heatmap_signal)) {
        pdf(file.path(save_dir, paste0("CellChat_compare_signaling_heatmap_", meta_col, ".pdf")),
            width = 10, height = 5)
        ComplexHeatmap::draw(p_heatmap_signal, heatmap_legend_side = "right")
        dev.off()
      }
      
      # circle plots for each condition
      pdf(file.path(save_dir, paste0("CellChat_compare_circle_", meta_col, ".pdf")),
          width = 4 * length(object.list), height = 4)
      circle_plot_fun()
      dev.off()
      
      invisible(list(
        cellchat_list        = cellchat_list,
        merged               = cellchat_merged,
        plot_compare         = p_compare,
        plot_heatmap         = p_heatmap,
        plot_heatmap_signal  = p_heatmap_signal,
        circle_plot_function = circle_plot_fun
      ))
    }
