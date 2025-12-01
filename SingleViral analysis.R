##1.viral load scoring

res <- ViralLoadScoring(obj = data, viral_genes = virus_gene, 
                        neg_cells = control_cells,  
                        q_baseline = 0, scale_type = "MAD",
                        winsor_q = 0.995)
data$viral_load <- res$viral_load_unweighted

FeaturePlot(data, features = "viral_load", cols = c("grey","red"))+
  ggtitle("Viral load scoring")+
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )
FeaturePlot(data, features = "viral_load", cols = c("grey","red"),
            split.by = "timepoint")


##2.TargetInfectionLevel
out <- TargetInfectionLevel(
  obj       = data,
  score_col = "viral_load",
  n_levels  = 5,
  method    = "quantile",
  labels    = c("Uninfected", "Low infection","Mederate infection","Mederate to high", "High infection")  # ⬅️ 自定义标签
)
data <- out$object
table(data$infection_label)
table(data$infection_level)
table(data$timepoint,data$infection_level)

DimPlot(data, group.by = "infection_level", label = TRUE,cols =  c(
  "4" = "#E64B35",
  "3" = "#F3A13B",
  "2" = "#44C37F",
  "1" = "#4DBBD5",
  "0" = "#BC80BD"
))+
  ggtitle("Infection level")+
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

DimPlot(data, group.by = "infection_label", label = TRUE,cols =  c(
  "High infection" = "#E64B35",
  "Mederate to high" = "#F3A13B",
  "Mederate infection" = "#44C37F",
  "Low infection" = "#4DBBD5",
  "Uninfected" = "#BC80BD"
))+
  ggtitle("Infection label")+
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )



##3.InfectionDynamic

res_gsea <- InfectionDynamic(
  obj        = data,
  group_col  = "infection_level",
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
  time_col = "timepoint",
  time_include = c("24hpi","48hpi","72hpi","96hpi")
)
res_gsea$plots[["1"]]
res_gsea$plots[["2"]] 
res_gsea$plots[["3"]] 
res_gsea$plots[["4"]] 


##4.CellVariableInfection

new_colors <- c(
  "CD4+Tcells"             = "#E64B35",
  "CD8+Tcells"             = "#4DBBD5",
  "DNTcells"               = "#00A087",
  "DPTcells"               = "#3C5488",
  "Th1"                    = "#F39B7F",
  "Th2"                    = "#7E6148",
  "Th17"                   = "#BB0021",
  "Tfh"                    = "#91D1C2",
  "Treg"                   = "#DC0000",
  "NK_cell"                = "#7E64A6",
  "B_cell"                 = "#00A1D5",
  "stem/progenitor_cells"  = "#FFC53D"
)

res <- CellVariableInfection(
  obj       = data,
  group_col = "manual_label",   # 或 "timepoint"
  score_col = "viral_load",
  level_col = "infection_level",
  zero_tol  = 0,
  palette   = new_colors,             # 可传入自定义颜色向量
  lvl_colors= c("0"="#7FB3D5","1"="#5DADE2","2"="turquoise","3"="#fdb462","4"="#D74C4C"),
  save_dir  = "infection_summary_plots"
)

res$plots$infection_proportion
res$plots$mean_viral_load
res$plots$infection_level_counts_plot
res$plots$infection_level_percent_plot
data1 <- subset(data,manual_label %in% c("CD4+Tcells","CD8+Tcells","DPTcells",
                                         "DNTcells","Th1","Th2"))
data2 <- subset(data,manual_label %in% c("Tfh","Treg","Th17","NK_cell","B_cell","stem/progenitor_cells"))
FeaturePlot(data1, features = "viral_load", cols = c("grey","red"),
            split.by = "manual_label")
FeaturePlot(data2, features = "viral_load", cols = c("grey","red"),
            split.by = "manual_label")
DimPlot(data, reduction = "umap", group.by = "infection_level", label = FALSE,
        split.by = "manual_label",cols = c(
          "4" = "#E64B35",
          "3" = "#F3A13B",
          "2" = "#44C37F",
          "1" = "#4DBBD5",
          "0" = "#BC80BD"
        ),ncol = 6)+
  ggtitle("Infection level") #15*5

##5.CellCycleByInfection

lvl_colors <- c(
  "#7FB3D5", 
  "#5DADE2", 
  "#48C9B0", 
  "#F5B041",  
  "#D35400"   
)
res_cc <- CellCycleByInfection(
  obj        = data,
  level_col  = "infection_level",
  save_dir   = "results/cell_cycle",
  title_phase  = "Cell-cycle phase (%) by infection level",
  title_violin = "Cell-cycle scores by infection level",
  lvl_colors = lvl_colors,
  show_points  = TRUE
)

res_cc$violin_plot
res_cc$phase_plot


##6.InfectionVariableGenes

res <- InfectionVariableGenes(
  obj            = data,
  timepoint_col  = "timepoint",
  manual_col     = "manual_label",
  level_col      = "infection_level",
  time_values    = c("24hpi","48hpi","72hpi","96hpi"),
  manual_include = c("CD4+Tcells","CD8+Tcells","DPTcells",
                     "DNTcells","Th1","Th2","Tfh","Treg","Th17"),
  case_levels    = 1:4,
  control_levels = 0,
  virus_genes    = virus_gene,
  save_dir       = "1BU1102/Tcells_DEG",
  save_csv       = "Tcells_DEG.csv",
  save_plot      = "Tcells_Volcano.pdf",
  plot_title     = "T cells(24-96hpi): level 1-4 vs level 0",
  top_n_pval     = 45,
  min.pct = 0.5
)


##7.VirusGeneDynamic

res <- VirusGeneDynamic(
  obj           = data,
  level_col     = "infection_level",
  time_col      = "timepoint",
  time_include  = c("24hpi","48hpi","72hpi","96hpi"),
  level_include = c("1","2","3","4"),
  viral_genes   = virus_gene,
  assay         = "RNA",
  slot          = "data",
  save_dir      = "infection_viral_summary",
  angle_col     = 315,
  fontsize = 7
  
)
res$avg_expr
res$heatmap

##8.CellChatCompare

library(CellChat)
data$infection_status <- ifelse(data$viral_load > 0, "infected", "uninfected")

res_cc <- CellChatCompare(
  obj          = data,
  meta_col     ="infection_status",
  meta_include = c("uninfected", "infected"),
  time_col     = "timepoint",
  timepoints   = c("24hpi","48hpi","72hpi","96hpi"),
  group_col    = "manual_label",
  group_include = c("CD4+Tcells","CD8+Tcells","DPTcells",
                    "DNTcells","Th1","Th2","Tfh","NK_cell","B_cell"),  # 可选
  assay        = "RNA",
  slot         = "data",
  organism     = "human",
  min_cells    = 10,
  workers      = 16,
  save_dir     = "CellChat_TP_compare",
  vi           = 1
)

res_cc$plot_heatmap_signal
res_cc$plot_compare
res_cc$plot_heatmap

##9.CellChatByCellInfection

res_ci <- CellChatByCellInfection(
  obj = data,
  label_col = "manual_label",
  score_col = "viral_load",
  timepoint_col = "timepoint",
  time_values = c("24hpi","48hpi","72hpi","96hpi"),
  infection_threshold = 0,            
  labels_use = c("CD4+Tcells","CD8+Tcells","DPTcells",
                 "DNTcells","Th1","Th2","Tfh","Treg","Th17","NK_cell","B_cell"),
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

RunTraceVirus(
  mode             = "virus",#or"host"
  run_wgcna        = FALSE,
  n_thread         = 20,
  genome_fasta     = "hhv6.fa",
  genome_gtf       = "hhv6.gtf",
  host_fasta       = "hg38.fa",
  host_gtf         = "hg38.gtf",
  host_virus_fasta = "host_virus.fa",
  host_virus_gtf   = "host_virus.gtf",
  output_dir       = "./Virus_Track_Output",
  output_dir_host  = "./Host_Virus_Output",
  rsem_index_name  = "hhv6_RSEM",
  star_index_dir   = "hhv6_star_index",
  fastq_dir        = "./",
  paired_end       = FALSE,
  traits_file      = "/path/to/traits.csv",
  work_dir         = "singleviral_tracevirus_test"
)
