library(Seurat)
library(SingleR)
library(SingleViral)
library(CellChat)
library(SoupX)
library(DoubletFinder)

preanalysis(base.dir = "/home/tang/CART_HHV6_reactivation/matrix", species = "human")

data$timepoint <- sapply(strsplit(colnames(data), "_"), `[`, 1)
data$timepoint <- gsub("[0-9]+$", "", sapply(strsplit(colnames(data), "_"), `[`, 1))
data$timepoint <- factor(data$timepoint,levels =c("control","2dpi","3dpi","5dpi","14dpi"))

data <- JoinLayers(data)

data <- NormalizeData(data)

data <- FindVariableFeatures(data,selection.method = "vst",
                             nfeatures = 2000)
all.genes <- rownames(data)
data <- ScaleData(data,features =all.genes)
data <- RunPCA(data,features = VariableFeatures(object = data))
data <- FindNeighbors(data,dims = 1:20)
data <- FindClusters(data,resolution = 0.5)
data <- RunUMAP(data,dims=1:20)  

mousecell <- celldex::ImmGenData()
clusters <- data$seurat_clusters
data_singler <- GetAssayData(data,layer = "data")
single3 <- SingleR(test = data_singler,ref = mousecell,
                   labels = mousecell$label.main)
plotDeltaDistribution(single3, ncol = 3)
single3$pruned.labels[is.na(single3$pruned.labels)] <- "low confidence"
single3$pruned.labels[single3$pruned.labels == "B cells, pro"] <- "B cells"

single3$pruned.labels[single3$pruned.labels == "Microglia"|single3$pruned.labels == "Mast cells"|
                        single3$pruned.labels == "Stem cells"] <- "low confidence"
data@meta.data$manual_label <- single3$pruned.labels
DimPlot(data,reduction = "umap",group.by = c("manual_label"),
        label = TRUE,repel = TRUE)
table(data$manual_label)

FeaturePlot(data, features = c("Epcam","Col1a1","Sftpa1","Pecam1","Ptprc"))
FeaturePlot(data, features = c("Ager","Rtkn2","Foxj1","Lamp3"))
Lymphoid_cells <- c("B cells","ILC","NKT","T cells","Tgd","NK cells")
Myeloid_cells <- c("Basophils","DC","Eosinophils","Macrophages","Monocytes","Neutrophils")
Stromal_cells <- c("Endothelial cells","Epithelial cells","Fibroblasts","Stromal cells","AT1","AT2")
data_Lymphoid <-  subset(data,manual_label %in% Lymphoid_cells)
Lymphoid_low_confidence <- WhichCells(data_Lymphoid, expression = Epcam > 0 | Col1a1 > 0 |
                                        Sftpb >0 | Pecam1 >0)
data$manual_label[Lymphoid_low_confidence] <- "low confidence"

data_myeloid <-  subset(data,manual_label %in% Myeloid_cells)
myeloid_low_confidence <- WhichCells(data_myeloid, expression = Epcam > 0 | Col1a1 > 0 |
                                       Sftpb >0 | Pecam1 >0)
data$manual_label[myeloid_low_confidence] <- "low confidence"

data_Stromal <-  subset(data,manual_label %in% Stromal_cells)
stromal_low_confidence <- WhichCells(data_Stromal,expression = Ptprc > 0)
data$manual_label[stromal_low_confidence] <- "low confidence"

DimPlot(data,reduction = "umap",label = TRUE,repel = TRUE)
DimPlot(data,reduction = "umap",group.by = c("manual_label"),
        label = TRUE,repel = TRUE)
data_at1 <-  subset(data,manual_label %in% c("Epithelial cells")&
                      seurat_clusters %in% c("15"))
data_at2 <-  subset(data,manual_label %in% c("Epithelial cells")&
                      seurat_clusters %in% c("9"))
data_stromal <-  subset(data,manual_label %in% c("Stromal cells"))
at1 <- WhichCells(object = data_at1)
at2 <- WhichCells(object = data_at2)
smooth <- WhichCells(data_stromal, expression = Acta2 > 1 & Myh11 > 1 &
                       Tagln >1 )
data$manual_label[at1] <- "AT1"
data$manual_label[at2] <- "AT2"
data$manual_label[smooth] <- "smooth muscle"
data <-  subset(data,!manual_label %in% c("low confidence"))
#重新标准化
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
table(data$manual_label)
cell_colors <- c(
  
  # --- 肺泡上皮 ---
  "AT1" = "#E64B35",        # 砖红
  "AT2" = "#F39B7F",        # 珊瑚橙
  "Epithelial cells" = "#FFD166",  # 金黄
  # --- 内皮 / 间质 ---
  "Endothelial cells" = "#4DBBD5",  # 湖蓝
  "Fibroblasts"       = "#00A087",  # 青绿
  "Stromal cells"     = "#9A7DCF",  # 淡紫
  "smooth muscle"    = "#E18727",  # 橘棕
  # --- 免疫：髓系 ---
  "Monocytes"   = "#91D1C2",  # 深蓝
  "Macrophages" = "#2E9FDF",  # 天蓝
  "Neutrophils" = "#DC0000",  # 鲜红
  "DC"          = "#FC4E07",  # 淡蓝
  # --- 免疫：淋巴系 ---
  "B cells"   = "#B09C85",   # 沙色
  "T cells"   = "#Fad",   # 粉紫
  "NK cells"  = "#0072B2",   # 薄荷绿
  "NKT"       = "#E64B35",   # 棕色
  "ILC"       = "#B0E0E6",  
  "Tgd"       = "#FFFF99",   
  # --- 粒细胞 ---
  "Basophils"   = "#BC3C29", # 深砖红
  "Eosinophils" = "#7CAE00" # 草绿
)
level <- names(sort(table(data@meta.data$manual_label), decreasing = TRUE))

data@meta.data$manual_label <- factor(data@meta.data$manual_label,levels = level)
DimPlot(data,reduction = "umap",group.by = c("manual_label"),
        label = TRUE,repel = TRUE,cols = cell_colors)+
  ggtitle("Cell Annotation")


##1.ViralLoadScoring
virus_gene <- c("ORF1ab","S","M","E","N","ORF3a","ORF6","ORF7a",
                "ORF7b","ORF8","ORF10")
control_cells <- rownames(data@meta.data[data@meta.data$timepoint == "control", ])

res <- ViralLoadScoring(obj = data, viral_genes = virus_gene, 
                        neg_cells = control_cells,write_mode = "unweighted",
                        group.col = "timepoint",plot = TRUE,  
                        q_baseline = 0, scale_type = "MAD",
                        winsor_q = 0.995)
res$plots$global
res$plots$split
data$viral_load <- res$viral_load_unweighted

##2.TargetInfectionLevel
custom_level_cols <- c(
  "4" = "#E64B35",
  "3" = "#F3A13B",
  "2" = "#44C37F",
  "1" = "#4DBBD5",
  "0" = "#BC80BD"
)
custom_label_cols <- c(
  "High infection"      = "#E64B35",
  "Moderate to high"    = "#F3A13B",
  "Moderate infection" = "#44C37F",
  "Low infection"      = "#4DBBD5",
  "Uninfected"         = "#BC80BD"
)

res_1 <- TargetInfectionLevel(
  obj = data,
  score_col = "viral_load",
  method = "kmeans",
  n_levels = 5,
  labels = c("Uninfected", "Low infection", "Moderate infection",
             "Moderate to high", "High infection"),
  plot = TRUE,
  level_cols = custom_level_cols,
  label_cols = custom_label_cols
) 

res_1$plots$level_plot
res_1$plots$label_plot
data$infection_level <- res_1$levels
data$infection_label <- res_1$labels

##3.InfectionDynamic

res_gsea <- InfectionDynamic(
  obj        = data,
  species    = "mouse",
  level_col  = "infection_level",
  ref_level  = 0,
  min.pct    = 0.1,
  logfc.threshold = 0,
  ranks_from = "avg_log2FC",
  top_up     = 10,
  top_down   = 10,
  col_up     = "#D74C4C",
  col_down   = "#3B6FB3",
  drop_prefix = TRUE,
  save_plot  = TRUE,
  outdir     = "GSEA_results",
  group_col = "timepoint",
  group_include = c("2dpi","3dpi","5dpi"),
  cell_type = "manual_label",
  celltype_include = Myeloid_cells,
  title_prefix = "Myeloid_cells"
)
res_gsea$plots[["1"]] # 10.5*8
res_gsea$plots[["2"]] 
res_gsea$plots[["3"]] 
res_gsea$plots[["4"]] 

##4.CellVariableInfection
res_2 <- CellVariableInfection(
  obj       = data,
  group_col = "manual_label",  
  score_col = "viral_load",
  level_col = "infection_level",
  zero_tol  = 0,
  palette   = NULL,            
  lvl_colors= c("0"="#7FB3D5","1"="turquoise","2"="#FFFF99","3"="#fdb462","4"="#D74C4C"),
  save_dir  = "infection_summary_plots",
  show_stack_count_label = FALSE,
  show_stack_percent_label = FALSE,
  split_ncol = 5
)

res_2$plots$infection_proportion
res_2$plots$mean_viral_load
res_2$plots$infection_level_counts_plot
res_2$plots$infection_level_percent_plot
res_2$plots$umap_viral_load_split
res_2$plots$umap_infection_level_split

##5.VirusGeneDynamic

res_3 <- VirusGeneDynamic(
  obj           = data,
  meta_col     = "manual_label",
  group_col      = "timepoint",
  group_include  = c("2dpi","3dpi","5dpi"),
  meta_include = NULL,
  viral_genes   = virus_gene,
  assay         = "RNA",
  layer          = "data",
  save_dir      = "infection_viral_summary",
  angle_col     = 315,
  fontsize = 7
)
res_3$avg_expr
res_3$heatmap

##6.CellCycleByInfection

lvl_colors <- c(
  "#7FB3D5", 
  "#5DADE2", 
  "#48C9B0", 
  "#F5B041",  
  "#D35400"   
)

res_cc <- CellCycleByInfection(
  obj        = data,
  species    = "mouse",
  level_col  = "infection_level",
  save_dir   = "results/cell_cycle",
  title_phase  = "Stromal cells:Cell-cycle phase (%) by infection level",
  title_violin = "Stromal cells:Cell-cycle scores by infection level",
  lvl_colors = lvl_colors,
  cell_type = "manual_label",
  celltype_include = Stromal_cells,
  show_points  = FALSE
)

res_cc$violin_plot
res_cc$phase_plot

##7.InfectionVariableGenes

res_4 <- InfectionVariableGenes(
  obj            = data,
  group_col      = "timepoint",
  cell_type      = "manual_label",
  level_col      = "infection_level",
  group_include  = c("2dpi","3dpi","5dpi"),
  celltype_include = Myeloid_cells,
  case_levels    = 1:4,
  control_levels = 0,
  virus_genes    = virus_gene,
  save_dir       = "SARScov2/Myeloid_cells_DEG",
  save_csv       = "Myeloid_cells_DEG.csv",
  save_plot      = "Myeloid_cells_Volcano.pdf",
  plot_title     = "Myeloid_cells(2-5dpi): level 1-4 vs level 0",
  top_n_pval     = 45,
  min.pct = 0
)

##8.CellChatCompare

data$infection_status <- ifelse(data$viral_load > 0, "infected", "uninfected")

##8.CellChatCompare

data$infection_status <- ifelse(data$viral_load > 0, "infected", "uninfected")

res_cc <- CellChatCompare(
  obj          = data,
  compare_col     ="timepoint",
  compare_1 =    c("control"),
  compare_2 =   c("2dpi","3dpi","5dpi"),
  group_col     = "infection_level",
  group_include   = c("0","1","2","3","4"),
  cell_type    = "manual_label",
  celltype_include = Myeloid_cells,  
  assay        = "RNA",
  slot         = "data",
  organism     = "human",
  min_cells    = 5 ,
  save_dir     = "CellChat_TP_compare",
  group_names = c("control","2-5dpi")
)