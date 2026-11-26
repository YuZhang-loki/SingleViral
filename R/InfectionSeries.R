#' Infection-level progression clustering for viral and host genes with ggplot2 visualization
#'
#' @description
#' This function performs an infection-progression style clustering analysis by treating
#' discrete \code{infection_level} values (e.g. 0/1/2/3/4) as an ordered "process axis".
#' Viral genes and host genes are clustered separately using Mfuzz soft clustering.
#'
#' Workflow:
#' \enumerate{
#'   \item Subset cells by \code{infection_level} (optional).
#'   \item Treat each infection level as one bin (NO quantile binning).
#'   \item Compute mean expression per gene within each infection level.
#'   \item Filter low-information genes (expressed in too few levels).
#'   \item Z-score per gene across levels.
#'   \item Mfuzz clustering separately for viral and host genes.
#'   \item Publication-ready ggplot2 trend plots for each cluster (viral & host).
#'   \item Export cluster gene lists (ranked by membership).
#' }
#'
#' @param obj A \strong{Seurat} object.
#' @param level_col Character. Column name in \code{obj@meta.data} storing infection levels.
#' @param level_include Optional vector. Keep only these infection levels (e.g. c(0,1,2,3,4)).
#' @param assay Character. Seurat assay to use (default \code{"RNA"}).
#' @param slot Character. Expression slot (default \code{"data"}).
#'
#' @param virus_genes Character vector of viral gene symbols (must match assay feature names).
#' @param host_genes Optional character vector of host gene symbols.
#'   If \code{NULL}, host genes are defined as \code{setdiff(all_genes, virus_genes)}.
#'
#' @param min_cells_per_level Integer. Minimum number of cells required in each infection level.
#' @param min_expr_levels Integer. Keep genes that are expressed (>0) in at least this many levels.
#'
#' @param drop_mt Logical. Remove mitochondrial genes from host gene set.
#' @param mt_pattern Character regex for mitochondrial genes (default \code{"^MT-"}).
#' @param drop_ribo Logical. Remove ribosomal genes from host gene set.
#' @param ribo_pattern Character regex for ribosomal genes (default \code{"^RP[SL]"}).
#'
#' @param k_virus Integer. Number of Mfuzz clusters for viral genes.
#' @param k_host Integer. Number of Mfuzz clusters for host genes.
#' @param m_virus Numeric. Mfuzz fuzzifier for viral genes.
#'   If \code{NULL}, estimated automatically by \code{Mfuzz::mestimate()}.
#' @param m_host Numeric. Mfuzz fuzzifier for host genes.
#'   If \code{NULL}, estimated automatically by \code{Mfuzz::mestimate()}.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @param export_clusters Logical. Whether to export per-cluster gene lists.
#' @param export_dir Character. Output directory for cluster gene lists.
#' @param membership_cut Numeric. Export only genes with membership >= this cutoff.
#'
#' @param plot_dir Optional character. If not NULL, save plots to this folder as PDF.
#' @param plot_prefix Character. Prefix used in plot filenames.
#' @param show_gene_individual Logical. Whether to overlay thin gene-level lines.
#' @param max_genes_individual Integer. Maximum number of genes drawn per cluster.
#'
#' @param verbose Logical. Print progress messages.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{bins}: data.frame of infection levels and cell counts used
#'   \item \code{expr_bin}: raw mean expression (genes x levels)
#'   \item \code{expr_virus_z}, \code{expr_host_z}: scaled matrices used for clustering
#'   \item \code{virus}, \code{host}: Mfuzz results (clusters, membership, m, mfuzz_obj)
#'   \item \code{plots}: ggplot objects (viral_cluster_plot, host_cluster_plot)
#'   \item \code{cluster_tables}: list of per-cluster gene tables (virus, host)
#' }
#'
#' @examples
#' \dontrun{
#' res <- InfectionLevelClusterVH(
#'   obj = data,
#'   level_col = "infection_level",
#'   virus_genes = readLines("virus_gene.txt"),
#'   min_cells_per_level = 20,
#'   k_virus = 6,
#'   k_host = 12,
#'   plot_dir = "plots"
#' )
#' res$plots$viral_cluster_plot
#' }
#'
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom Matrix rowMeans
#' @importFrom dplyr filter arrange group_by summarise mutate select slice_head ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon facet_wrap
#'   theme_classic labs ggsave
#' @importFrom readr write_csv
#'
#' @export
InfectionSeries <- function(
    obj,
    level_col = "infection_level",
    level_include = NULL,
    assay = "RNA",
    slot  = "data",
    virus_genes,
    host_genes = NULL,
    min_cells_per_level = 10,
    min_expr_levels = 2,
    drop_mt = TRUE,
    mt_pattern = "^MT-",
    drop_ribo = FALSE,
    ribo_pattern = "^RP[SL]",
    k_virus = 6,
    k_host  = 12,
    m_virus = NULL,
    m_host  = NULL,
    seed = 123,
    export_clusters = TRUE,
    export_dir = "infections_series_clusters",
    membership_cut = 0,
    plot_dir = NULL,
    plot_prefix = "infections_series_level",
    show_gene_individual = TRUE,
    max_genes_individual = 200,
    verbose = TRUE
) {
  # Check optional dependencies only when this function is called.
  optional_pkgs <- c("Biobase", "Mfuzz")
  missing_optional <- optional_pkgs[
    !vapply(optional_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_optional) > 0L) {
    stop(
      "Optional package(s) required for InfectionSeries(): ",
      paste(missing_optional, collapse = ", "),
      ". Please install them before using this function.",
      call. = FALSE
    )
  }

  ## ------------------------------- ##
  ## 0) Basic checks
  ## ------------------------------- ##
  stopifnot(inherits(obj, "Seurat"))                              
  if (!level_col %in% colnames(obj@meta.data))                     
    stop("Column '", level_col, "' not found in meta.data.")
  
  if (missing(virus_genes) || length(virus_genes) == 0)           
    stop("'virus_genes' is required and cannot be empty.")
  
  
  pkgs_need <- c("Seurat", "Matrix", "dplyr", "tidyr", "ggplot2", "readr")
  miss_pkg <- pkgs_need[!vapply(pkgs_need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_pkg) > 0) stop("Missing packages: ", paste(miss_pkg, collapse = ", "))
  
  ## ------------------------------- ##
  ## 1) Prepare infection_level and subset cells
  ## ------------------------------- ##
  md <- obj@meta.data                                                
  lvl_raw <- md[[level_col]]                                         
  lvl_chr <- as.character(lvl_raw)                                   
  lvl_num <- suppressWarnings(as.numeric(lvl_chr))                   
  
  
  if (all(is.na(lvl_num))) stop("infection_level cannot be coerced to numeric. Please check: ", level_col)
  
  
  keep <- rep(TRUE, nrow(md))                                        
  if (!is.null(level_include)) {                                     
    keep <- keep & (lvl_chr %in% as.character(level_include))        
  }
  
  
  keep <- keep & !is.na(lvl_num)                                     
  cells_use <- rownames(md)[keep]                                    
  if (length(cells_use) == 0) stop("No cells left after level_include / NA filtering.")
  
  
  obj_sub <- subset(obj, cells = cells_use)                          
  md2 <- obj_sub@meta.data                                           
  lvl2 <- suppressWarnings(as.numeric(as.character(md2[[level_col]])))
  if (any(is.na(lvl2))) stop("Some infection_level values are NA after subsetting. Please check metadata.")
  
  ## ------------------------------- ##
  ## 2) Drop infection levels with too few cells
  ## ------------------------------- ##
  level_counts <- as.data.frame(table(lvl2), stringsAsFactors = FALSE)
  colnames(level_counts) <- c("level", "n_cells")                    
  level_counts$level_num <- as.numeric(level_counts$level)           
  level_counts <- level_counts[order(level_counts$level_num), ]      
  
  
  good_levels <- level_counts$level[level_counts$n_cells >= min_cells_per_level]
  if (length(good_levels) < 2) {
    stop("Less than 2 infection levels have >= min_cells_per_level. Lower min_cells_per_level.")
  }
  
  
  keep2 <- as.character(md2[[level_col]]) %in% as.character(good_levels) 
  cells_use2 <- rownames(md2)[keep2]                                    
  obj_sub <- subset(obj_sub, cells = cells_use2)                        
  md2 <- obj_sub@meta.data                                               
  lvl2 <- suppressWarnings(as.numeric(as.character(md2[[level_col]])))   
  lvl_fac <- factor(lvl2, levels = sort(unique(lvl2)), ordered = TRUE)   
  md2$.__inf_level <- lvl_fac                                            
  
  
  bins <- as.data.frame(table(md2$.__inf_level), stringsAsFactors = FALSE)
  colnames(bins) <- c("level", "n_cells")                                 
  bins$level_num <- as.numeric(as.character(bins$level))                  
  
  if (verbose) {
    message("Infection levels used (after filtering):")
    print(bins)
  }
  
  ## ------------------------------- ##
  ## 3) Compute mean expression per infection level
  ## ------------------------------- ##
  Seurat::DefaultAssay(obj_sub) <- assay                                   
  expr <- Seurat::GetAssayData(obj_sub, assay = assay, slot = slot)        
  
  
  level_levels <- levels(md2$.__inf_level)                                 
  expr_bin <- sapply(level_levels, function(lv) {                          
    cells_lv <- rownames(md2)[md2$.__inf_level == lv]                      
    Matrix::rowMeans(expr[, cells_lv, drop = FALSE])                       
  })
  
  
  if (is.null(dim(expr_bin))) expr_bin <- matrix(expr_bin, ncol = 1)
  
  colnames(expr_bin) <- level_levels                                       
  all_genes <- rownames(expr_bin)                                          
  
  ## ------------------------------- ##
  ## 4) Define viral genes and host genes
  ## ------------------------------- ##
  virus_genes_use <- intersect(virus_genes, all_genes)                     
  if (length(virus_genes_use) < 2)
    stop("Too few virus genes matched assay features. Check virus_genes and assay rownames.")
  
  
  if (is.null(host_genes)) {
    host_genes_use <- setdiff(all_genes, virus_genes_use)                  
  } else {
    host_genes_use <- intersect(host_genes, all_genes)                    
  }
  
  
  if (drop_mt)   host_genes_use <- host_genes_use[!grepl(mt_pattern, host_genes_use)]
  if (drop_ribo) host_genes_use <- host_genes_use[!grepl(ribo_pattern, host_genes_use)]
  
  ## ------------------------------- ##
  ## 5) Filter low-information genes across levels
  ## ------------------------------- ##
  keep_by_levels <- function(mat, min_levels = 2) {                        
    nz <- rowSums(mat > 0) >= min_levels                                   
    mat[nz, , drop = FALSE]                                                
  }
  
  mat_v <- keep_by_levels(expr_bin[virus_genes_use, , drop = FALSE], min_expr_levels) 
  mat_h <- keep_by_levels(expr_bin[host_genes_use,  , drop = FALSE], min_expr_levels) 
  
  if (nrow(mat_v) < 2) stop("Too few viral genes remain after filtering.")
  if (nrow(mat_h) < 2) stop("Too few host genes remain after filtering.")
  
  ## ------------------------------- ##
  ## 6) Z-score scaling per gene across infection levels
  ## ------------------------------- ##
  zscore <- function(x) {                                                  
    mu <- rowMeans(x)                                                      
    sdv <- apply(x, 1, stats::sd)                                           
    sdv[sdv == 0] <- 1                                                     
    (x - mu) / sdv                                                         
  }
  
  expr_virus_z <- zscore(as.matrix(mat_v))                                 
  expr_host_z  <- zscore(as.matrix(mat_h))                                  
  
  ## ------------------------------- ##
  ## 7) Mfuzz clustering
  ## ------------------------------- ##
  set.seed(seed)                                                            
  
  run_mfuzz <- function(mat_z, k, m = NULL) {                                
    eset <- Biobase::ExpressionSet(assayData = mat_z)                        
    if (is.null(m)) m <- Mfuzz::mestimate(eset)                              
    cl <- Mfuzz::mfuzz(eset, c = k, m = m)                                   
    list(
      mfuzz_obj  = cl,                                                      
      m          = m,                                                       
      clusters   = cl$cluster,                                               
      membership = cl$membership                                            
    )
  }
  
  if (verbose) message("Clustering viral genes with Mfuzz ...")
  virus_res <- run_mfuzz(expr_virus_z, k_virus, m_virus)
  
  if (verbose) message("Clustering host genes with Mfuzz ...")
  host_res  <- run_mfuzz(expr_host_z,  k_host,  m_host)
  
  ## ------------------------------- ##
  ## 8) Export per-cluster gene lists (ranked by membership)
  ## ------------------------------- ##
  make_cluster_tables <- function(res) {                                     
    cl  <- res$clusters                                                     
    mem <- res$membership                                                  
    
    out <- vector("list", length = max(cl))                                 
    names(out) <- as.character(seq_len(max(cl)))                             
    
    for (k in seq_len(max(cl))) {                                           
      genes_k <- names(cl)[cl == k]                                         
      if (length(genes_k) == 0) {                                           
        out[[as.character(k)]] <- data.frame()                              
        next
      }
      df <- data.frame(
        gene       = genes_k,
        cluster    = k,
        membership = mem[genes_k, k],
        stringsAsFactors = FALSE
      )
      df <- df |>
        dplyr::filter(membership >= membership_cut) |>
        dplyr::arrange(dplyr::desc(membership))
      out[[as.character(k)]] <- df
    }
    out
  }
  
  cluster_tables <- list(
    virus = make_cluster_tables(virus_res),
    host  = make_cluster_tables(host_res)
  )
  
  if (export_clusters) {                                                     
    if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)    
    dir_v <- file.path(export_dir, "virus")                                  
    dir_h <- file.path(export_dir, "host")                                   
    dir.create(dir_v, showWarnings = FALSE, recursive = TRUE)
    dir.create(dir_h, showWarnings = FALSE, recursive = TRUE)
    
    
    for (k in names(cluster_tables$virus)) {
      df <- cluster_tables$virus[[k]]
      if (nrow(df) == 0) next
      readr::write_csv(df, file.path(dir_v, sprintf("virus_cluster_%02d_genes.csv", as.integer(k))))
    }
    
    
    for (k in names(cluster_tables$host)) {
      df <- cluster_tables$host[[k]]
      if (nrow(df) == 0) next
      readr::write_csv(df, file.path(dir_h, sprintf("host_cluster_%02d_genes.csv", as.integer(k))))
    }
    
    if (verbose) message("✅ Cluster gene lists exported to: ", normalizePath(export_dir))
  }
  
  ## ------------------------------- ##
  ## 9) Prepare ggplot2 publication-style cluster trend plots
  ## ------------------------------- ##
  
  prep_plot_df <- function(mat_z, res, kind = c("virus", "host")) {
    kind <- match.arg(kind)
    df <- as.data.frame(mat_z)                                               
    df$gene <- rownames(df)                                                  
    df_long <- df |>
      tidyr::pivot_longer(                                                   
        cols = -gene,
        names_to = "level",
        values_to = "z"
      )
    df_long$level_num <- as.numeric(as.character(df_long$level))             
    df_long$cluster <- res$clusters[df_long$gene]                            
    df_long$membership <- NA_real_                                           
    
    
    cl <- res$clusters
    mem <- res$membership
    df_long$membership <- mem[cbind(df_long$gene, as.character(df_long$cluster))]
    
    
    df_sum <- df_long |>
      dplyr::group_by(cluster, level_num) |>
      dplyr::summarise(
        mean_z = mean(z, na.rm = TRUE),
        se_z   = stats::sd(z, na.rm = TRUE) / sqrt(dplyr::n()),
        .groups = "drop"
      )
    
    list(long = df_long, summary = df_sum)
  }
  
  pv <- prep_plot_df(expr_virus_z, virus_res, "virus")
  ph <- prep_plot_df(expr_host_z,  host_res,  "host")
  
  
  plot_cluster_trends <- function(plot_pack, title) {
    df_long <- plot_pack$long
    df_sum  <- plot_pack$summary
    
    col_gene  <- "#F2A541"  
    col_trend <- "#B22222"  
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(
        data = df_sum,
        ggplot2::aes(x = level_num, ymin = mean_z - se_z, ymax = mean_z + se_z),
        fill  = col_trend,
        alpha = 0.18
      ) +
      ggplot2::geom_line(
        data = df_sum,
        ggplot2::aes(x = level_num, y = mean_z),
        color = col_trend,
        linewidth = 1.25
      ) +
      ggplot2::geom_point(
        data = df_sum,
        ggplot2::aes(x = level_num, y = mean_z),
        color = col_trend,
        size = 2.1
      ) +
      ggplot2::facet_wrap(~cluster, scales = "free_y") +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::labs(
        title = title,
        x = "Infection level",
        y = "Z-scored mean expression"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(face = "bold", size = 12)
      )
    
    if (isTRUE(show_gene_individual)) {
      df_long2 <- df_long |>
        dplyr::group_by(cluster, gene) |>
        dplyr::summarise(
          membership = max(membership, na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::group_by(cluster) |>
        dplyr::arrange(dplyr::desc(membership)) |>
        dplyr::slice_head(n = max_genes_individual) |>
        dplyr::ungroup() |>
        dplyr::select(cluster, gene)
      
      df_draw <- dplyr::inner_join(
        df_long,
        df_long2,
        by = c("cluster", "gene")
      )
      
      p <- p +
        ggplot2::geom_line(
          data = df_draw,
          ggplot2::aes(x = level_num, y = z, group = gene),
          color = col_gene,
          alpha = 0.1,          
          linewidth = 0.35
        )
    }
    
    p
  }
  
  p_virus <- plot_cluster_trends(pv, "Viral gene expression dynamics along infection levels")
  p_host  <- plot_cluster_trends(ph, "Host gene expression dynamics along infection levels")
  
  
  if (!is.null(plot_dir)) {
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(
      filename = file.path(plot_dir, paste0(plot_prefix, "_virus_clusters.pdf")),
      plot = p_virus, width = 11, height = 7, useDingbats = FALSE
    )
    ggplot2::ggsave(
      filename = file.path(plot_dir, paste0(plot_prefix, "_host_clusters.pdf")),
      plot = p_host, width = 11, height = 7, useDingbats = FALSE
    )
    if (verbose) message("✅ Plots saved to: ", normalizePath(plot_dir))
  }
  
  ## ------------------------------- ##
  ## 10) Return everything
  ## ------------------------------- ##
  invisible(list(
    bins = bins,
    expr_bin = expr_bin,
    expr_virus_z = expr_virus_z,
    expr_host_z  = expr_host_z,
    virus = virus_res,
    host  = host_res,
    plots = list(
      viral_cluster_plot = p_virus,
      host_cluster_plot  = p_host
    ),
    cluster_tables = cluster_tables
  ))
}

