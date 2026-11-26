library(tidyverse)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(dplyr)
library(ggrepel)
library(Seurat)
library(SingleR)
library(SingleViral)
base.dir <- "/home/tang/HIV_GSE266455/matrix"
sample.names <- list.dirs(base.dir, full.names = FALSE, recursive = FALSE)
sample.names <- sample.names[sample.names != ""]
for (sid in sample.names) {
  
  data.dir <- file.path(base.dir, sid, "outs", "filtered_feature_bc_matrix")
  cat("正在处理：", sid, "  路径：", data.dir, "\n")
  
  counts <- Read10X(data.dir = data.dir)
  
  obj <- CreateSeuratObject(
    counts      = counts,
    project     = sid,
    min.cells   = 3,
    min.features = 200
  )
  
  obj$sample <- sid   # 给 meta.data 加一个样本名（之后很有用）
  
  # 在全局环境中生成和你习惯一样的对象名：data_2dpi_blood3 等
  assign(paste0("data_", sid), obj, envir = .GlobalEnv)
}

data <- merge(
  data_8026_r5,
  y = c(
    data_8026_r6,data_8027_r1,data_8027_r2,
    data_8027_r5,data_8027_r6,data_145_r1,data_145_r2
  ),
  add.cell.ids = c("8026_r5","8026_r6","8027_r1","8027_r2",
                   "8027_r5","8027_r6","145_r1","145_r2"),
  project = "hiv"
)


data$patient <- gsub("_r.*", "", data$orig.ident)

data <- JoinLayers(data)

data[["percent.mt"]] <- PercentageFeatureSet(data,pattern = "^MT-")

data <- subset(data,subset=nFeature_RNA>200 & percent.mt<60)

data <- NormalizeData(data)

data <- FindVariableFeatures(data,selection.method = "vst",
                             nfeatures = 2000)

all.genes <- rownames(data)
data <- ScaleData(data,features =all.genes)

data <- RunPCA(data,features = VariableFeatures(object = data))

data <- FindNeighbors(data,dims = 1:20)
data <- FindClusters(data,resolution = 0.5)
data <- RunUMAP(data,dims=1:20)  

data@meta.data$manual_label <- "unclassified cells"

CD4Tcells <- WhichCells(object = data,expression = MS4A1==0 & MZB1==0)

data$manual_label[CD4Tcells] <- "CD4+Tcells"

data <- subset(data,manual_label %in% c("CD4+Tcells"))

data <- NormalizeData(data)

data <- FindVariableFeatures(data,selection.method = "vst",
                             nfeatures = 2000)

all.genes <- rownames(data)
data <- ScaleData(data,features =all.genes)

data <- RunPCA(data,features = VariableFeatures(object = data))

data <- FindNeighbors(data,dims = 1:20)
data <- FindClusters(data,resolution = 0.5)
data <- RunUMAP(data,dims=1:20)  

DimPlot(data,reduction = "umap",label = TRUE)
DimPlot(data,reduction = "umap", group.by = "patient", label = TRUE,repel = TRUE)

virus_gene <- c(
  "gag","pol","vif","vpr","tat","rev","vpu","env","vef"
)


##1.ViralLoadScoring

result <- ViralLoadScoring(
  obj = data,
  viral_genes = virus_gene,
  neg_cells = NULL,
  background_method = "none",
  group.col = "patient",
  plot = TRUE
)

data <- result$obj

result$plots$global
result$plots$split

##2.TargetInfectionLevel

custom_level_cols <- c(
  "2" = "#E64B35",
  "1" = "#F3A13B",
  "0" = "#BC80BD"
)
custom_label_cols <- c(
  "High infection"      = "#E64B35",
  "Low infection"    = "#F3A13B",
  "Uninfected"         = "#BC80BD"
)

res_1 <- TargetInfectionLevel(
  obj = data,
  score_col = "viral_load",
  method = "quantile",
  n_levels = 3,
  labels = c("Uninfected","Low infection","High infection"),
  plot = TRUE,
  level_cols = custom_level_cols,
  label_cols = custom_label_cols
) 

res_1$plots$level_plot #4*3
res_1$plots$label_plot #6*4
data$infection_level <- res_1$levels
data$infection_label <- res_1$labels

##3.InfectionPseudobulk

res_pb <- InfectionPseudobulk(
  obj = data,
  
  patient_col = "patient",
  level_col = "infection_level",
  viral_col = "viral_load",
  
  ref_level = 0,
  case_levels = 1:2,
  
  use_viral_overlap = FALSE,
  
  min_cells_per_sample = 10,
  min_patients = 2,
  
  run_gsea = TRUE,
  species = "human",
  
  plot_volcano = TRUE,
  label_volcano = TRUE,
  
  volcano_top_up = 15,
  volcano_top_down = 15,
  volcano_label_size = 3.5,
  
  plot_pca = TRUE,
  plot_gsea = TRUE,
  
  save_results = TRUE,
  outdir = "Pseudobulk_results",
  padj_cut = 0.1,
  logfc_cut = 0,
  viral_genes = virus_gene
)

# 查看分析模式
res_pb$analysis_mode

# 查看配对患者数量及分析状态
res_pb$comparison_summary

# 查看每位患者的感染级别分布
res_pb$level_distribution

# 查看Pseudobulk样本信息
res_pb$sample_metadata

# 查看Level 1差异表达结果
head(res_pb$DE[["1"]])

# 查看Level 1 GSEA结果
res_pb$GSEA[["1"]]

# 绘制Level 1火山图
res_pb$volcano_plots[["2"]] #5.5*4

# 绘制Level 1 GSEA图
res_pb$gsea_plots[["2"]]#10.5*6.5

# 绘制PCA
res_pb$pca_plot