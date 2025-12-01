library(tidyverse)
library(ggrepel)
library(SingleR)
library(ggplot2)
library(Seurat)
library(SingleViral)
data$timepoint <- sapply(strsplit(colnames(data), "_"), `[`, 1)
data$timepoint <- gsub("[0-9]+$", "", sapply(strsplit(colnames(data), "_"), `[`, 1))
data@meta.data$timepoint <- factor(data@meta.data$timepoint,levels =
                                     c("control","24hpi","48hpi",
                                       "72hpi","96hpi"))
data[["percent.mt"]] <- PercentageFeatureSet(data,pattern = "^MT-")
VlnPlot(data,features = c("nFeature_RNA","nCount_RNA","percent.mt"),
        ncol = 3)
plot1 <- FeatureScatter(data,feature1 = "nFeature_RNA",
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(data,feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA")
plot1+plot2
data <- subset(data,subset=nFeature_RNA>200&percent.mt<50)
data <- NormalizeData(data)
data <- FindVariableFeatures(data,selection.method = "vst",
                             nfeatures = 2000)
all.genes <- rownames(data)
data <- ScaleData(data,features =all.genes)
data <- RunPCA(data,features = VariableFeatures(object = data))
VizDimLoadings(data,dims = 1:12,reduction = "pca")
DimHeatmap(data,dims = 1:12,cells = 500,balanced = TRUE)
data <- JackStraw(data,num.replicate = 100)
data <- ScoreJackStraw(data,dims = 1:20)
data <- FindNeighbors(data,dims = 1:20)
data <- FindClusters(data,resolution = 0.5)
data <- RunUMAP(data,dims=1:20)
DimPlot(data,reduction = "umap",label = TRUE)
DimPlot(data,reduction = "umap", group.by = "timepoint", label = TRUE, pt.size = 0.5)


data <- JoinLayers(data)
virus_gene <- readLines("virus_gene.txt")
control_cells <- rownames(data@meta.data)[data$timepoint == "control"]

clusters <- data$seurat_clusters
data_singler <- GetAssayData(data,slot = "data")
humancell <- celldex::HumanPrimaryCellAtlasData()
single3 <- SingleR(test = data_singler,ref = humancell,
                   labels = humancell$label.main)
single3$labels[single3$labels != "T_cells"&single3$labels !="B_cell"&
                 single3$labels!="NK_cell"] <- "stem/progenitor_cells"
data@meta.data$manual_label <- single3$labels
data_Tcell <- subset(data,manual_label %in% c("T_cells"))
dnTcells <- WhichCells(object = data_Tcell,expression = CD8A == 0 & CD4 == 0)
dpTcells <- WhichCells(object = data_Tcell,expression = CD8A > 0 & CD4 > 0)
CD4Tcells <- WhichCells(object = data_Tcell,expression = CD4 > 0 & CD8A == 0)
CD8Tcells <- WhichCells(object = data_Tcell,expression = CD8A > 0 & CD4 == 0)
data$manual_label[CD8Tcells] <- "CD8+Tcells"
data$manual_label[CD4Tcells] <- "CD4+Tcells"
data$manual_label[dpTcells] <- "DPTcells"
data$manual_label[dnTcells] <- "DNTcells"
data1 <- subset(data,manual_label %in% c("CD4+Tcells","CD4+CD8+Tcells"))
Th1_cells <- WhichCells(data1, expression = TBX21 > 0 )
Th2_cells <- WhichCells(data1, expression = GATA3 > 0 )
Th17_cells <- WhichCells(data1, expression = RORC > 0)
Tfh <- WhichCells(data1, expression = BCL6 > 0)
Treg_cells <- WhichCells(data1,expression = IL2RA > 0 & FOXP3 > 0)
data$manual_label[Treg_cells] <- "Treg"
data$manual_label[Th1_cells] <- "Th1"
data$manual_label[Th2_cells] <- "Th2"
data$manual_label[Th17_cells] <- "Th17"
data$manual_label[Tfh] <- "Tfh"
data@meta.data$manual_label <- factor(data@meta.data$manual_label,levels =
                                        c("CD4+Tcells","CD8+Tcells","DPTcells",
                                          "DNTcells","Th1","Th2","Tfh","Treg","Th17","NK_cell","B_cell","stem/progenitor_cells"))

new_colors <- c(
  "#E64B35", # B_cell  深红
  "#4DBBD5", # CD4-CD8-Tcells  湖蓝
  "#00A087", # CD4+CD8+Tcells  青绿
  "#3C5488", # CD4+Tcells  靛蓝
  "#F39B7F", # CD8+Tcells  淡橙红
  "#8491B4", # exhausted_Tcells  灰蓝
  "#91D1C2", # NK_cell  浅青绿
  "#DC0000", # stem/progenitor_cells  鲜红
  "#7E6148", # Th1_cells  棕色
  "#B09C85", # Th17_cells  沙色
  "#9A77CF", # Th2_cells  紫色
  "#F7C530"  # Treg_cells  金黄
)
DimPlot(data, group.by = "manual_label", cols = new_colors)+
  ggtitle("Cell Annotation")