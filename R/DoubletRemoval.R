#' Perform quality control, clustering, and doublet removal for a single Seurat object
#'
#' @description
#' This function performs a standard single-cell RNA-seq preprocessing workflow
#' on a single \code{Seurat} object, including:
#' \itemize{
#'   \item mitochondrial content calculation (species-aware),
#'   \item quality control filtering,
#'   \item normalization, feature selection, dimensionality reduction and clustering,
#'   \item automatic doublet detection and removal using \pkg{DoubletFinder}.
#' }
#'
#' The calculation of mitochondrial gene percentage depends on the specified
#' \code{species}:
#' \itemize{
#'   \item \code{"human"}: mitochondrial genes are detected using the \code{"^MT-"} prefix;
#'   \item \code{"mouse"}: mitochondrial genes are loaded from the built-in
#'         \file{mt_genes.txt} file bundled in the \pkg{SingleViral} package;
#'   \item \code{"other"}: mitochondrial genes must be provided explicitly via
#'         the \code{mt.genes} argument.
#' }
#'
#' @param data A \code{Seurat} object containing raw UMI counts.
#' @param mt.genes Character vector of mitochondrial gene names. Required only
#'   when \code{species = "other"}.
#' @param species Character. Species type used to determine how mitochondrial
#'   genes are defined. One of \code{"human"}, \code{"mouse"}, or \code{"other"}.
#' @param min.features Minimum number of detected features per cell.
#' @param max.features Maximum number of detected features per cell.
#' @param max.counts Maximum total UMI counts per cell.
#' @param max.percent.mt Maximum allowed percentage of mitochondrial reads.
#' @param nfeatures Number of variable features to select.
#' @param dims_use Integer vector specifying PCA dimensions to use.
#' @param resolution Resolution parameter for graph-based clustering.
#' @param sct Logical. Whether SCTransform normalization was used upstream
#'   (passed to \pkg{DoubletFinder}).
#' @param doublet_rate Expected doublet rate (fraction of total cells).
#' @param pN Numeric. Proportion of artificial doublets used in DoubletFinder.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return
#' A filtered \code{Seurat} object with doublets removed. The returned object
#' contains updated clustering, UMAP embeddings, and QC metadata.
#'
#' @details
#' This function automatically detects the appropriate \code{DoubletFinder}
#' API (v3 or legacy) to ensure compatibility across package versions.
#'
#' @seealso
#' \code{\link[DoubletFinder]{doubletFinder}},
#' \code{\link[Seurat]{PercentageFeatureSet}}
#'
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA
#'   FindNeighbors FindClusters RunUMAP PercentageFeatureSet
#' @importFrom DoubletFinder paramSweep summarizeSweep find.pK
#'
#' @export
ProcessOneSeurat <- function(
    data,
    mt.genes = NULL,
    species = c("human", "mouse", "other"),
    min.features = 200,
    max.features = 6000,
    max.counts = 20000,
    max.percent.mt = 60,
    nfeatures = 2000,
    dims_use = 1:20,
    resolution = 0.5,
    sct = TRUE,
    doublet_rate = 0.075,
    pN = 0.25,
    verbose = TRUE
) {
  stopifnot(inherits(data, "Seurat"))
  species <- match.arg(species)

  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat not installed.")
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) stop("DoubletFinder not installed.")

  # -------------------- percent.mt (species-aware) --------------------
  if (species == "human") {
    # Seurat default for human
    data[["percent.mt"]] <- Seurat::PercentageFeatureSet(data, pattern = "^MT-")

  } else {
    # mouse: read from SingleViral built-in mt_genes.txt
    if (species == "mouse") {
      mt_file <- system.file("mt_genes.txt", package = "SingleViral")
      if (mt_file == "") stop("Cannot find built-in 'mt_genes.txt' in package 'SingleViral'.")
      mt.genes <- readLines(mt_file)
    }

    # other: must provide mt.genes
    if (species == "other") {
      if (is.null(mt.genes) || length(mt.genes) == 0) {
        stop("For species = 'other', you must provide 'mt.genes'.")
      }
    }

    mt_use <- intersect(mt.genes, rownames(data))
    if (length(mt_use) == 0) {
      warning("No mt.genes found in object. percent.mt set to 0.")
      data[["percent.mt"]] <- 0
    } else {
      data[["percent.mt"]] <- Seurat::PercentageFeatureSet(data, features = mt_use)
    }
  }

  # -------------------- QC --------------------
  data <- subset(
    data,
    subset = nFeature_RNA > min.features &
      nFeature_RNA < max.features &
      nCount_RNA < max.counts &
      percent.mt < max.percent.mt
  )

  # -------------------- Standard workflow --------------------
  data <- Seurat::NormalizeData(data, verbose = FALSE)
  data <- Seurat::FindVariableFeatures(data, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  data <- Seurat::ScaleData(data, features = rownames(data), verbose = FALSE)
  data <- Seurat::RunPCA(data, features = Seurat::VariableFeatures(data), verbose = FALSE)
  data <- Seurat::FindNeighbors(data, dims = dims_use, verbose = FALSE)
  data <- Seurat::FindClusters(data, resolution = resolution, verbose = FALSE)
  data <- Seurat::RunUMAP(data, dims = dims_use, verbose = FALSE)

  # -------------------- DoubletFinder sweep (v3/old compatible) --------------------
  if (exists("paramSweep_v3", where = asNamespace("DoubletFinder"), inherits = FALSE)) {
    sweep.res.list <- DoubletFinder::paramSweep_v3(data, PCs = dims_use, sct = sct)
    sweep.stats    <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn          <- DoubletFinder::find.pK(sweep.stats)
  } else {
    sweep.res.list <- DoubletFinder::paramSweep(data, PCs = dims_use, sct = sct)
    sweep.stats    <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn          <- DoubletFinder::find.pK(sweep.stats)
  }

  pK_bcmvn <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  homotypic.prop <- DoubletFinder::modelHomotypic(data$seurat_clusters)
  nExp_poi <- round(doublet_rate * ncol(data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

  if (exists("doubletFinder_v3", where = asNamespace("DoubletFinder"), inherits = FALSE)) {
    data <- DoubletFinder::doubletFinder_v3(
      data,
      PCs = dims_use, pN = pN, pK = pK_bcmvn,
      nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = sct
    )
  } else {
    data <- DoubletFinder::doubletFinder(
      data,
      PCs = dims_use, pN = pN, pK = pK_bcmvn,
      nExp = nExp_poi.adj, reuse.pANN = NULL, sct = sct
    )
  }

  df_cols <- grep("^DF\\.classifications", colnames(data@meta.data), value = TRUE)
  if (length(df_cols) == 0) stop("DoubletFinder output DF.classifications_* not found.")
  df_col <- df_cols[length(df_cols)]

  n_before <- ncol(data)
  meta <- data@meta.data
  keep_cells <- rownames(meta)[meta[[df_col]] == "Singlet"]
  data <- subset(data, cells = keep_cells)
  n_after <- ncol(data)

  if (verbose) {
    message("species=", species,
            " | pK=", pK_bcmvn,
            " | homotypic=", round(homotypic.prop, 4),
            " | nExp=", nExp_poi, " -> adj ", nExp_poi.adj,
            " | cells ", n_before, " -> ", n_after,
            " (", df_col, ")")
  }

  return(data)
}

#' Batch-process multiple Seurat objects in place
#'
#' @description
#' This function applies \code{\link{ProcessOneSeurat}} to multiple
#' \link[Seurat]{Seurat} objects stored in an environment, and overwrites
#' each object with its processed version.
#'
#' The function is designed for large projects with many samples
#' (e.g. multiple 10x libraries) and minimizes memory usage by processing
#' objects sequentially.
#'
#' @param obj_names Character vector of object names to process.
#'   Each name must correspond to a \link[Seurat]{Seurat} object in \code{envir}.
#' @param mt.genes Character vector of mitochondrial gene names used for QC.
#' @param envir Environment containing the Seurat objects.
#'   Defaults to \code{.GlobalEnv}.
#' @param ... Additional parameters passed directly to
#'   \code{\link{ProcessOneSeurat}}.
#'
#' @return
#' Invisibly returns \code{NULL}.
#' Each Seurat object specified in \code{obj_names} is replaced in \code{envir}
#' by its processed version.
#'
#' @details
#' This function modifies objects \strong{in place}.
#' Use with caution if the original objects need to be preserved.
#'
#' Garbage collection (\code{gc()}) is invoked after each sample to reduce
#' memory usage.
#'
#' @examples
#' \dontrun{
#' obj_names <- ls(pattern = "^data_")
#'
#' BatchProcessSeurat_InPlace(
#'   obj_names = obj_names,
#'   mt.genes  = mt_genes,
#'   min.features = 300,
#'   max.percent.mt = 50,
#'   resolution = 0.6
#' )
#' }
#'
#' @seealso
#' \link{ProcessOneSeurat}
#'
#' @export
BatchProcessSeurat_InPlace <- function(
    obj_names,
    mt.genes,
    envir = .GlobalEnv,
    ...
) {
  stopifnot(is.character(obj_names), length(obj_names) > 0)
  
  for (nm in obj_names) {
    message("\n==== ", nm, " ====")
    obj <- get(nm, envir = envir)
    obj2 <- ProcessOneSeurat(obj, mt.genes = mt.genes, ...)
    assign(nm, obj2, envir = envir)  
    rm(obj, obj2)
    gc()
  }
  
  invisible(NULL)  
}
