library(Gviz)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(txdbmaker)

options(ucscChromosomeNames = FALSE)


axisTrack <- GenomeAxisTrack(col="#7e6fb0")
displayPars(axisTrack) <- list(
  cex = 0.8,
  fontcolor = "#7E6FB0"
)
gtf_path1 <- "Nanopore_U1102/canonical-1.gtf"
gtf_path2 <- "Nanopore_U1102/variant-1.gtf"
gtf_path3 <- "Nanopore_U1102/fusion-1.gtf"
txdb1 <- txdbmaker::makeTxDbFromGFF(gtf_path1,format = "gtf")
txdb2 <- txdbmaker::makeTxDbFromGFF(gtf_path2,format = "gtf")
txdb3 <- txdbmaker::makeTxDbFromGFF(gtf_path3,format = "gtf")

grTrack1 <- GeneRegionTrack(
  txdb1,
  format = "gtf",
  chromosome = "NC_001664.4",
  start = 26000,
  end = 48000,
  name = "Canonical",
  transcriptAnnotation = "gene_name",
  collapseTranscripts = FALSE,
  showId = TRUE,
  geneSymbol = TRUE,
  col = NA,
  fill = "#E64B35",
  background.panel = "#ffe",
  background.title = "#9E2F2F",
  stackHeight = 0.35,
)


grTrack2 <- GeneRegionTrack(
  txdb2,
  format = "gtf",
  chromosome = "NC_001664.4",
  start = 26000,
  end = 48000,
  name = "Variant",
  transcriptAnnotation = "gene_name",
  collapseTranscripts = FALSE,
  showId = TRUE,
  geneSymbol = TRUE,
  col = NA,
  fill = "#4DBBD5",
  background.panel = "#dbeeff",
  background.title = "darkblue",
  stackHeight = 0.35
)
grTrack3 <- GeneRegionTrack(
  txdb3,
  format = "gtf",
  chromosome = "NC_001664.4",
  end = 48000,
  name = "Fusion",
  transcriptAnnotation = "gene_name",
  collapseTranscripts = FALSE,
  showId = TRUE,
  geneSymbol = TRUE,
  col = NA,
  fill = "#7e6fb0",
  background.panel = "#FFCCFF",
  background.title = "#990099",
  stackHeight = 0.35
)

displayPars(grTrack1) <- list(cex.title = 1.0,fontsize.feature = 0.8)
displayPars(grTrack2) <- list(cex.title = 1.0,fontsize.feature = 0.8)
displayPars(grTrack3) <- list(cex.title = 1.0,fontsize.feature = 0.8)
bw_file <- "Nanopore_U1102/u1102.bw"

covTrack <- DataTrack(
  range = bw_file,
  genome ="custom",
  chromosome = "NC_001664.4",
  start = 26000,
  end = 48000,
  type = "histogram",
  name = "Coverage",
  col.histogram = "#91D1C2",
  background.title = "#339999",
)
displayPars(covTrack) <- list(cex.title = 1.0)

plotTracks(
  list(covTrack,grTrack1,grTrack2,grTrack3,axisTrack),
  chromosome = "NC_001664.4",
  from = 26000,
  to = 48000,
  transcriptAnnotation = "symbol",
  main = ""
)
