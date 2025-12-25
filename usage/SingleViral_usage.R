library(Seurat)
library(SingleR)
library(SingleViral)
library(CellChat)
library(SoupX)
library(DoubletFinder)

preanalysis(base.dir = "/home/tang/file/matrix_1", species = "mouse") #Compatible with 10x Genomics single-cell sequencing data

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

##1.ViralLoadScoring
path <- system.file("virus_gene.txt", package = "SingleViral")
virus_gene <- readLines("virus_gene.txt")
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

T_cells <- c("naive CD4+T cells","CD8+T cells",
  "DNT","DPT","Tfh","Th1","Th2")

res_gsea <- InfectionDynamic(
  obj        = data,
  species    = "human",
  level_col  = "infection_level",
  ref_level  = 0,
  min.pct    = 0.1,
  logfc.threshold = 0,
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
  celltype_include = T_cells,
  title_prefix = "T cells"
)
res_gsea$plots[["1"]] # 11*8
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
  show_stack_percent_label = FALSE
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
  meta_col     = "infection_level",
  group_col      = "timepoint",
  group_include  = c("24hpi","48hpi","72hpi","96hpi"),
  meta_include = c("1","2","3","4"),
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

lvl_colors <- c("#7FB3D5", "#5DADE2", "#48C9B0", "#F5B041", "#D35400")

res_cc <- CellCycleByInfection(
  obj        = data,
  species    = "mouse",
  level_col  = "infection_level",
  save_dir   = "results/cell_cycle",
  title_phase  = "Cell-cycle phase (%) by infection level",
  title_violin = "Cell-cycle scores by infection level",
  lvl_colors = lvl_colors,
  cell_type = "manual_label",
  celltype_include = NULL,
  show_points  = FALSE
)

res_cc$violin_plot
res_cc$phase_plot

##7.InfectionVariableGenes

res <- InfectionVariableGenes(
  obj            = data,
  group_col      = "timepoint",
  cell_type      = "manual_label",
  level_col      = "infection_level",
  group_include  = c("24hpi","48hpi","72hpi","96hpi"),
  celltype_include = Stromal_cells,
  case_levels    = 1:4,
  control_levels = 0,
  virus_genes    = virus_gene,
  save_dir       = "SARScov2/Stromal_cells_DEG",
  save_csv       = "Stromal_cells_DEG.csv",
  save_plot      = "Stromal_cells_Volcano.pdf",
  plot_title     = "Stromal_cells(2-5dpi): level 1-4 vs level 0",
  top_n_pval     = 45,
  min.pct = 0
)

##8.CellChatCompare

data$infection_status <- ifelse(data$viral_load > 0, "infected", "uninfected")

res_cc <- CellChatCompare(
  obj          = data,
  compare_col     ="timepoint",
  compare_1 =    c("24hpi"),
  compare_2 =   c("96hpi"),
  group_col     = "infection_level",
  group_include   = c("0","1","2","3","4"),
  cell_type    = "manual_label",
  celltype_include = c("B cells","NK cells","naive CD4+T cells","CD8+T cells",
                       "DNT","DPT","Tfh","Th1","Th2"),  
  assay        = "RNA",
  slot         = "data",
  organism     = "human",
  min_cells    = 5 ,
  save_dir     = "CellChat_TP_compare",
  group_names = c("24hpi","96hpi")
)

##9.CellChatByCellInfection

res_ci <- CellChatByCellInfection(
  obj = data,
  cell_type = "manual_label",
  score_col = "viral_load",
  group_col = "timepoint",
  group_include = c("24hpi","48hpi","72hpi","96hpi"),
  infection_threshold = 0,            
  celltype_include =  c("Treg","B cells","NK cells","naive CD4+T cells","CD8+T cells",
                  "DNT","DPT","Tfh","Th1","Th17","Th2"),
  organism = "human",
  min_cells = 5,
  save_dir = "CellChat_CellInfection",
  per_pathway_top = 20,
  mfrow = c(5,4),
  label_size = 0.5,
  label_size_merge = 0.5
)
res_ci$outdir

##10.RunTraceVirus
#Place the data to be analyzed into the working directory.
#For single-end sequencing data, rename the sample as ***_R1.fastq.gz (e.g., 101_R1.fastq.gz).
#For paired-end sequencing data, rename the files as ***_R1.fastq.gz and ***_R2.fastq.gz (e.g., 101_R1.fastq.gz and 101_R2.fastq.gz).
#After renaming, run the following function.
RunTraceVirus(
  mode             = "virus",#or"host"
  run_wgcna        = FALSE,
  n_thread         = 20,
  genome_fasta     = "hhv6.fa",
  genome_gtf       = "hhv6.gtf",
  host_fasta       = "hg38.fa",
  host_gtf         = "hg38.gtf",
  output_dir       = "./Virus_Track_Output",
  output_dir_host  = "./Host_Virus_Output",
  rsem_index_name  = "hhv6_RSEM",
  star_index_dir   = "hhv6_star_index",
  fastq_dir        = "./",
  paired_end       = FALSE,
  traits_file      = "/path/to/traits.csv",
  work_dir         = "singleviral_tracevirus_test",
  sequence_length  = 149
)
