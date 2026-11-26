library(tidyverse)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(dplyr)
library(ggrepel)
library(Seurat)
library(SingleR)
library(SingleViral)
umi <- read.delim(file.choose(),
                  stringsAsFactors = F,row.names = 1)
umi_u1102 <- read.delim(file.choose(),
                        stringsAsFactors = F,row.names = 1)

gene.df <- bitr (row.names(umi),
                 fromType = "GENEID",
                 toType = c("SYMBOL"),
                 OrgDb = EnsDb.Hsapiens.v86)
colnames(gene.df)=c("gene","symbol")
umi <- mutate(umi,gene=rownames(umi))
umi <- merge(umi,gene.df)

umi=umi[!duplicated(umi$symbol),]
umi <- na.omit(umi )
rownames(umi) = umi$symbol
umi <- umi[,-386]
umi <- umi[,-1]
rownames(umi)[rownames(umi) == "U2"] <- "U2_hg"
rownames(umi)[rownames(umi) == "U3"] <- "U3_hg"
rownames(umi)[rownames(umi) == "U4"] <- "U4_hg"
rownames(umi)[rownames(umi) == "U6"] <- "U6_hg"
rownames(umi)[rownames(umi) == "U7"] <- "U7_hg"
#矩阵纵向合并
umi2 <- bind_rows(umi, umi_u1102)
umi2[is.na(umi2)]=0

data_c4 <- CreateSeuratObject(counts = umi,
                             project = "control4",
                             min.cells = 3,min.features = 200)

data_96h5 <- CreateSeuratObject(counts = umi2,
                                project = "96h5",
                                min.cells = 3,min.features = 200)

data <- merge(data_c1, y = c(data_c2, data_c3,data_24h1,data_24h2,data_24h3,data_24h4,data_48h1,data_48h2,data_48h3,data_48h4,data_72h1,data_72h2,
                             data_72h3,data_72h4,data_96h1,data_96h2,data_96h3,data_96h4), 
              add.cell.ids = c("control1", "control2", "control3","24hpi1","24hpi2","24hpi3","24hpi4","48hpi1","48hpi2","48hpi3",
                               "48hpi4","72hpi1","72hpi2","72hpi3","72hpi4","96hpi1","96hpi2","96hpi3","96hpi4"), 
              project = "data_HHV6")

data$timepoint <- sapply(strsplit(colnames(data), "_"), `[`, 1)

data$timepoint <- gsub("[0-9]+$", "", sapply(strsplit(colnames(data), "_"), `[`, 1))

data$timepoint <- factor(data$timepoint,levels =c("control","24hpi","72hpi","96hpi"))

data[["percent.mt"]] <- PercentageFeatureSet(data,pattern = "^MT-")

data <- JoinLayers(data)

data <- subset(data,subset=nFeature_RNA>200 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt<60)

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
DimPlot(data,reduction = "umap", group.by = "timepoint", label = TRUE,repel = TRUE)

data_singler <- GetAssayData(data,layer = "data")

humancell <- celldex::HumanPrimaryCellAtlasData()
single3 <- SingleR(test = data_singler,ref = humancell,
                   labels = humancell$label.main)
table(single3$labels)
single3$labels[single3$labels == "Pro-B_cell_CD34-"] <- "B_cell"
single3$labels[single3$labels != "T_cells"&single3$labels !="B_cell"&
                 single3$labels!="NK_cell"] <- "unclassified cells"
data@meta.data$manual_label <- single3$labels
data_Tcell <- subset(data,manual_label %in% c("T_cells"))
dnTcells <- WhichCells(object = data_Tcell,expression = CD8A == 0 & CD4 == 0)
dpTcells <- WhichCells(object = data_Tcell,expression = CD8A > 0 & CD4 > 0)
CD4Tcells <- WhichCells(object = data_Tcell,expression = CD4 > 0 & CD8A == 0)
CD8Tcells <- WhichCells(object = data_Tcell,expression = CD8A > 0 & CD4 == 0)
data$manual_label[CD8Tcells] <- "CD8+Tcells"
data$manual_label[CD4Tcells] <- "naive CD4+Tcells"
data$manual_label[dpTcells] <- "DPTcells"
data$manual_label[dnTcells] <- "DNTcells"
data1 <- subset(data,manual_label %in% c("naive CD4+Tcells","DNTcells","DPTcells"))
Th1_cells <- WhichCells(data1, expression = TBX21 > 0 )
Th2_cells <- WhichCells(data1, expression = GATA3 > 0 )
Th17_cells <- WhichCells(data1, expression = RORC > 0)
Th22_cells <- WhichCells(data1,expression =  AHR > 0)
Tfh <- WhichCells(data1, expression = BCL6 > 0)
Treg_cells <- WhichCells(data1,expression = IL2RA > 0 & FOXP3 > 0)
data$manual_label[Th2_cells] <- "Th2"
data$manual_label[Th22_cells] <- "Th22"
data$manual_label[Tfh] <- "Tfh"
data$manual_label[Th1_cells] <- "Th1"
data$manual_label[Treg_cells] <- "Treg"
data$manual_label[Th17_cells] <- "Th17"

data@meta.data$manual_label <- factor(data@meta.data$manual_label,levels =
                                        c("naive CD4+Tcells","CD8+Tcells","DPTcells",
                                          "DNTcells","Th1","Th2","Tfh","Th17","Th22","Treg","NK_cell","B_cell","unclassified cells"))

new_colors <- c(
  "#E64B35",
  "#4DBBD5", 
  "#00A087", 
  "#Fad", 
  "#F39B7F", 
  "#FFFF99",
  "#91D1C2", 
  "#DC0000",
  "#2E9FDF",
  "#8491B4",
  "#9A77CF", 
  "#F7C530",
  "#3C5488"
)

DimPlot(data, group.by = "manual_label", cols = new_colors)+
  ggtitle("CBMC HHV-6B")

virus_gene <- c(
  "DR4","DR3","DR5","DR6","DR8",
  "U2","U2R","U3","U4","U7","U6","U10","U9",
  "U11","U11R","U13","U12","U13R","U14","U14R",
  "U15","U17","U17.5","U18","U18R",
  "U19","U20","U20.5","U21","U22","U23","U23R",
  "U25R","U24","U25","U26","U27","U27R","U28",
  "U29","U30","U30R","U31R","U31","U32","U33",
  "U35R","U34","U35","U36","U36R","U37","U38",
  "U38R","U39","U39R","U40","U40R","U41","U41R",
  "U46","U42","U43","U44","U45","U47","U48",
  "U48R","U49","U50","U50R","U51","U52","U52R",
  "U53","U53R","U54","U54R","U55","U57","U56",
  "U58","U59","U60/66","U62","U63","U64","U65",
  "U76R","U66R","U67","U68","U69","U70","U71","U72",
  "U73","U73R","U74","U75","U76","U77",
  "U77/79","U78","U79",
  "U81","U81.5","U82","U83","U83R","U84",
  "U85R","U85","IE2","U88","IE1","U91",
  "U94","U95","U95R","gQ2","gQ1","U100R"
)


##1.ViralLoadScoring
control_cells <- rownames(data@meta.data[data@meta.data$timepoint == "control", ])

result <- ViralLoadScoring(
  obj = data,
  viral_genes = virus_gene,
  neg_cells = control_cells,
  background_method = "negative_reference",
  group.col = "timepoint",
  plot = TRUE
)

data <- result$obj

result$plots$viral_load
result$plots$global
result$plots$split

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

res_1$plots$level_plot #4*3
res_1$plots$label_plot #6*4
data$infection_level <- res_1$levels
data$infection_label <- res_1$labels

##3.InfectionDynamic

res_gsea <- InfectionDynamic(
  obj        = data,
  species    = "human",
  level_col  = "infection_level",
  ref_level  = 0,
  min.pct    = 0.1,
  logfc.threshold = 0.1,
  ranks_from = "avg_log2FC",
  top_up     = 15,
  top_down   = 15,
  col_up     = "#D74C4C",
  col_down   = "#3B6FB3",
  drop_prefix = TRUE,
  save_plot  = TRUE,
  outdir     = "GSEA_results",
  group_col = "timepoint",
  group_include = c("24hpi","48hpi","72hpi","96hpi"),
  cell_type = "manual_label",
  celltype_include = c("naive CD4+Tcells","CD8+Tcells","DPTcells",
                     "DNTcells","Th1","Th2","Th17","Th22","Tfh","Treg"),
  title_prefix = NULL
)
res_gsea$plots[["1"]] # 12*8
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
  celltype_include = c("naive CD4+Tcells","CD8+Tcells","DPTcells",
                       "DNTcells","Th1","Th2","Tfh","Th17","Th22","Treg","NK_cell","B_cell")
)

res_2$plots$mean_viral_load_violin #8*4
res_2$plots$infection_proportion
res_2$plots$mean_viral_load
res_2$plots$infection_level_counts_plot
res_2$plots$infection_level_percent_plot
res_2$plots$umap_viral_load_split
res_2$plots$umap_infection_level_split

##5.VirusGeneDynamic
res_3 <- VirusGeneDynamic(
  obj           = data,
  meta_col     = "infection_level",
  group_col      = "timepoint",
  group_include  =c("24hpi","72hpi","96hpi"),
  meta_include = c("1","2","3","4"),
  viral_genes   = virus_gene,
  assay         = "RNA",
  layer          = "data",
  save_dir      = "infection_viral_summary",
  angle_col     = 315,
  fontsize = 7
)
res_3$avg_expr
res_3$heatmap #9*9

res_3 <- VirusGeneDynamic(
  obj           = data,
  meta_col     = "manual_label",
  group_col      = "timepoint",
  group_include  = c("24hpi","48hpi","72hpi","96hpi"),
  meta_include = c("naive CD4+Tcells","CD8+Tcells","DPTcells",
                   "DNTcells","Th1","Th2","Tfh","Th17","Th22","Treg","NK_cell","B_cell"),
  viral_genes   = virus_gene,
  assay         = "RNA",
  layer          = "data",
  save_dir      = "infection_viral_summary",
  angle_col     = 315,
  fontsize = 7
)

res_3$avg_expr
res_3$heatmap #9*9

##6.CellCycleByInfection

lvl_colors <- c("#7FB3D5", "#5DADE2", "#48C9B0", "#F5B041", "#D35400")

res_cc <- CellCycleByInfection(
  obj        = data,
  species    = "human",
  level_col  = "infection_level",
  save_dir   = "results/cell_cycle",
  title_phase  = "Cell-cycle phase (%) by infection level",
  title_violin = "Cell-cycle scores by infection level",
  lvl_colors = lvl_colors,
  cell_type = "manual_label",
  celltype_include = NULL,
  show_points  = FALSE
)

res_cc$violin_plot#8*4
res_cc$phase_plot

##7.InfectionVariableGenes

result <- InfectionVariableGenes(
  obj = data,
  group_include = c("24hpi", "48hpi", "72hpi", "96hpi"),
  celltype_include =  c("naive CD4+Tcells","CD8+Tcells","DPTcells",
                        "DNTcells","Th1","Th2","Th17","Th22","Tfh","Treg"),
  trend_top_n = 9,
  virus_genes = virus_gene,
  save_dir = "./results"
)

result$correlation_plot

result$top_positive
result$trend_plot_positive

result$top_negative
result$trend_plot_negative

result$significant_combined

##8.InfectionSeries
res <- InfectionSeries(
  obj = data,
  level_col = "infection_level",
  virus_genes = virus_gene,
  min_cells_per_level = 20,
  k_virus = 6,
  k_host  = 12,
  membership_cut = 0.3,
  export_dir = "infections_series_clusters",
  plot_dir = "infections_series_plots",
  show_gene_individual = TRUE
)

res$plots$viral_cluster_plot
res$plots$host_cluster_plot
res$cluster_tables

##9.GeneViralLoadCorrelation
res1 <- GeneViralLoadCorrelation(
  obj = data,
  genes = "NCL",
  group_col = "manual_label",
  group_include = c("naive CD4+Tcells","CD8+Tcells","DPTcells",
                    "DNTcells","Th1","Th2","Th22","Tfh","Treg"),
  viral_col = "viral_load",
  level_col = "infection_level",
  method = "pearson",
  show_cell_scatter = TRUE,
  show_cell_points = TRUE
)
res1$plot
res1$cor_table

