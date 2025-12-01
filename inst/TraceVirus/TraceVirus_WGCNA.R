options(timeout = 1000)

required_cran <- c("stringr", "reshape2","ggplot2")
                   
for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
  

required_bioc <- c("AnnotationDbi","GO.db", "preprocessCore", "impute")
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  install.packages("WGCNA", repos = "https://cloud.r-project.org")
}

library(WGCNA)
library(data.table)
library(reshape2)
library(stringr)
enableWGCNAThreads()
data <- read.table(file = "Virus-Host_FPKM_matrix.txt", sep='\t', row.names=1, header=T,
                       quote="", comment="", check.names=F)
m.mad <- apply(data, 1, mad)

dataExpr <- data[which(m.mad > 0), ]

dataExpr <- log2(dataExpr+1)

dataExpr <- t(dataExpr)

sampleTree <- hclust(dist(dataExpr), method = "average")

pdf("Sample clustering tree.pdf", width=12, height=9)

plot(sampleTree, main = "Sample clustering to detect outliers", xlab = "", sub = "")

dev.off()

nSamples <- nrow(dataExpr)

nGenes <- ncol(dataExpr)

softPower <- ifelse(nSamples < 20, 18,
                ifelse(nSamples < 30, 16,
                       ifelse(nSamples < 40, 14, 12)))

net <-  blockwiseModules(dataExpr, power = softPower, maxBlockSize = nGenes,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = "pearson",
                       maxPOutliers= 1, loadTOMs=TRUE,
                       saveTOMFileBase = "virusTOM",
                       verbose = 3)

moduleLabels = net$colors

moduleColors = labels2colors(moduleLabels)

label_color_map <- data.frame(
  Label = sort(unique(moduleLabels)),
  Color = unique(moduleColors[order(moduleLabels)])
)

write.table(label_color_map,
            file = "label_color_map.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
            
pdf("Sample clustering tree.pdf", width=12, height=9)

plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
                                          
dev.off()

MEs = net$MEs

MEs_col = MEs

colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
  
MEs_col = orderMEs(MEs_col)

pdf("Eigengene adjacency heatmap.pdf", width=12, height=9)

plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)

dev.off() 
                     
traitData <- read.csv(file = "traits.csv",row.names = 1)

modTraitCor = cor(MEs_col, traitData, use = "p")

modTraitP = corPvalueStudent(modTraitCor, nSamples) 

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")

dim(textMatrix) = dim(modTraitCor)

pdf("Module-trait relationships.pdf", width=12, height=9)

labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData),
               yLabels = colnames(MEs_col),
               cex.lab = 0.5,
               ySymbols = colnames(MEs_col), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = TRUE,
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
               
dev.off()
dataExpr <- as.data.frame(dataExpr)
modules <- unique(moduleColors)

module_gene_list <- lapply(modules, function(mod) {
  names(dataExpr)[moduleColors == mod]
})

names(module_gene_list) <- modules

for (mod in names(module_gene_list)) {
  genes <- module_gene_list[[mod]]
  fname <- paste0("WGCNA_genes_", mod, ".txt")
  write.table(genes, file = fname, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}


required_bioc <- c("clusterProfiler", "AnnotationDbi", "org.Hs.eg.db","DOSE")

for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DOSE)
library(cowplot)
library(ggplot2)
options(timeout = 1000)
# 设定路径
files <- list.files(pattern = "^WGCNA_genes_.*\\.txt$")

# 循环分析
for (file in files) {
  # 提取模块名称
  mod <- sub("WGCNA_genes_(.*)\\.txt", "\\1", file)
  cat(">>> Processing module:", mod, "\n")
  
  # 读取基因
  genes <- read.table(file)

  # ENSEMBL 转换为 ENTREZID
  gene_entrez <- bitr(genes$V1, fromType = "ENSEMBL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)

  if (nrow(gene_entrez) < 5) {
    cat("Module", mod, "has too few genes after ID conversion. Skipping...\n")
    next
  }

  # GO 富集分析
  ego <- enrichGO(gene         = gene_entrez$ENTREZID,
                  OrgDb        = org.Hs.eg.db,
                  keyType      = "ENTREZID",
                  ont          = "ALL",
                  pvalueCutoff = 0.2,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)

  # 保存 GO 结果
  write.csv(ego, file = paste0("GO_enrichment_", mod, ".csv"))

  # KEGG 富集分析
  ekegg <- enrichKEGG(gene         = gene_entrez$ENTREZID,
                      organism     = "hsa",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.5,
                      qvalueCutoff  = 0.5)

  # 保存 KEGG 结果
  write.csv(ekegg, file = paste0("KEGG_enrichment_", mod, ".csv"))

  # 可选：画图
   if (!is.null(ego) && nrow(ego) > 0) {
   go_plot <- barplot(ego, showCategory = 20, color = "pvalue",
                       label_format = 250, split = "ONTOLOGY",cex.axis=0.05) +
               facet_grid(ONTOLOGY ~ ., scale = "free")

    ggsave(filename = paste0("GO_barplot_", mod, ".pdf"),
           plot = go_plot,
           width = 12, height = 9)
  } else {
    cat("Module", mod, "has no GO result.\n")
  }


 if (!is.null(ekegg) && nrow(ekegg) > 0) {
    kegg_plot <- barplot(ekegg, showCategory = 20, color = "pvalue",
                         label_format = 250,cex.axis=0.05)

    ggsave(filename = paste0("KEGG_barplot_", mod, ".pdf"),
           plot = kegg_plot,
           width = 12, height = 9)
  } else {
    cat("Module", mod, "has no KEGG result.\n")
  }
}           
  

