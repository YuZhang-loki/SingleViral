#!/bin/bash

MODE=$1
SUBMODE=$2

if [ -z "$MODE" ]; then
  echo "Usage: bash TraceVirus.sh [virus|host] [WGCNA (optional)]"
  exit 1
fi

source Parameters.txt

# =======================
# VIRUS MODULE
# =======================
if [ "$MODE" == "virus" ]; then
  echo ">>> Running Virus Quantification Module"

mkdir -p $Output_directory
cd $Output_directory

# ===== [1] Build STAR index =====
echo ">>> Checking STAR index..."
if [ -d "$STAR_index_dir" ] && [ -f "$STAR_index_dir/SA" ]; then
  echo ">>> STAR index already exists, skipping..."
else
  echo ">>> Building STAR index"
  STAR --runThreadN $N_thread \
       --runMode genomeGenerate \
       --genomeDir $STAR_index_dir \
       --genomeFastaFiles ../$Genome_fasta \
       --sjdbGTFfile ../$Genome_gtf \
       --sjdbOverhang 149 \
       --genomeSAindexNbases 7
fi

# ===== [2] Build RSEM index =====
echo ">>> Checking RSEM index..."
if [ -f "${RSEM_index_name}.grp" ] && [ -f "${RSEM_index_name}.idx.fa" ]; then
  echo ">>> RSEM index already exists, skipping..."
else
  echo ">>> Building RSEM index"
  rsem-prepare-reference --gtf ../$Genome_gtf \
                         -p $N_thread \
                         ../$Genome_fasta \
                         $RSEM_index_name
fi

# ===== [3] STAR alignment + RSEM quantification =====
echo ">>> Running STAR alignment and RSEM quantification"

cd ../$Fastq_dir
SAMPLES=$(ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//' | sort)
cd - > /dev/null

for SAMPLE in $SAMPLES; do
  echo ">>> Processing sample: $SAMPLE"

  if [ "$Paired_end" == "TRUE" ]; then
    STAR --genomeDir $STAR_index_dir \
         --runThreadN $N_thread \
         --readFilesIn ../$Fastq_dir/${SAMPLE}_R1.fastq.gz ../$Fastq_dir/${SAMPLE}_R2.fastq.gz \
         --readFilesCommand zcat \
         --outFileNamePrefix ${SAMPLE}_ \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts \
         --limitBAMsortRAM 1002778363

    rsem-calculate-expression --paired-end \
         --alignments --no-bam-output -p $N_thread -q \
         ${SAMPLE}_Aligned.toTranscriptome.out.bam \
         $RSEM_index_name \
         ${SAMPLE}-RSEM

  else
    STAR --genomeDir $STAR_index_dir \
         --runThreadN $N_thread \
         --readFilesIn ../$Fastq_dir/${SAMPLE}_R1.fastq.gz \
         --readFilesCommand zcat \
         --outFileNamePrefix ${SAMPLE}_ \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts \
         --limitBAMsortRAM 1002778363

    rsem-calculate-expression \
         --alignments --no-bam-output -p $N_thread -q \
         ${SAMPLE}_Aligned.toTranscriptome.out.bam \
         $RSEM_index_name \
         ${SAMPLE}-RSEM
  fi
done

# ===== [4] Merge expression matrices and generate heatmap =====
echo ">>> Merging expression matrices and generating TPM heatmap"

Rscript --vanilla - <<EOF

required_packages <- c("data.table", "pheatmap")
installed <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(data.table)
library(pheatmap)

files <- list.files(pattern = "genes.results$")
tpm_list <- list()
count_list <- list()
fpkm_list <- list()

for (file in files) {
  dt <- fread(file)
  sample_id <- sub("-RSEM.*", "", file)
  sample_name <- paste0("sample_", sample_id)

  tpm_dt <- dt[, .(gene_id, temp_value = TPM)]
  setnames(tpm_dt, "temp_value", sample_name)
  tpm_list[[sample_name]] <- tpm_dt

  count_dt <- dt[, .(gene_id, temp_value = expected_count)]
  setnames(count_dt, "temp_value", sample_name)
  count_list[[sample_name]] <- count_dt

  fpkm_dt <- dt[, .(gene_id, temp_value = FPKM)]
  setnames(fpkm_dt, "temp_value", sample_name)
  fpkm_list[[sample_name]] <- fpkm_dt
}

merge_matrix <- function(lst) {
  merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), lst)
  gene_ids <- merged[['gene_id']]
  merged <- merged[, -1]
  merged[is.na(merged)] <- 0
  rownames(merged) <- gene_ids
  return(merged)
}

tpm_matrix <- merge_matrix(tpm_list)
count_matrix <- merge_matrix(count_list)
fpkm_matrix <- merge_matrix(fpkm_list)

fwrite(as.data.table(cbind(gene_id = rownames(count_matrix), count_matrix)),
       file = "Virus_Count_matrix.txt", sep = "\t")
fwrite(as.data.table(cbind(gene_id = rownames(tpm_matrix), tpm_matrix)),
       file = "Virus_TPM_matrix", sep = "\t")
fwrite(as.data.table(cbind(gene_id = rownames(fpkm_matrix), fpkm_matrix)),
       file = "Virus_FPKM_matrix.txt", sep = "\t")

pdf("Virus_TPM_heatmap.pdf", width=12, height=12)


tpm_matrix <- log2(tpm_matrix+1)

rownames(tpm_matrix) <- rownames(merge_matrix(tpm_list))

pheatmap(tpm_matrix,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         fontsize = 6,
         angle_col = 45,
         main = "Viral Gene Expression")
dev.off()
EOF

  echo ">>> Virus analysis complete."
  exit 0
fi

# =======================
# HOST MODULE
# =======================
if [ "$MODE" == "host" ]; then
  echo ">>> Running Host Quantification Module"
echo ">>> Merging host and virus genome and annotation"

cat $host_fasta $Genome_fasta > host_virus.fa

cat $host_gtf $Genome_gtf > host_virus.gtf

mkdir -p $Output_directory_host
cp traits.csv $Output_directory_host
cd $Output_directory_host

# ===== [1] Build STAR index =====
echo ">>> Checking STAR index..."
if [ -d "$STAR_index_dir" ] && [ -f "$STAR_index_dir/SA" ]; then
    echo ">>> STAR index already exists, skipping..."
else
    echo ">>> Building STAR index"
    STAR --runThreadN $N_thread \
         --runMode genomeGenerate \
         --genomeDir $STAR_index_dir \
         --genomeFastaFiles ../$host_virus_fasta \
         --sjdbGTFfile ../$host_virus_gtf \
         --sjdbOverhang 149
fi

# ===== [2] Build RSEM index =====
echo ">>> Checking RSEM index..."
if [ -f "${RSEM_index_name}.grp" ] && [ -f "${RSEM_index_name}.idx.fa" ]; then
    echo ">>> RSEM index already exists, skipping..."
else
    echo ">>> Building RSEM index"
    rsem-prepare-reference --gtf ../$host_virus_gtf \
                           -p $N_thread \
                           ../$host_virus_fasta \
                           $RSEM_index_name
fi

# ===== [3] STAR alignment + RSEM quantification =====
echo ">>> Running STAR alignment and RSEM quantification"

cd ../$Fastq_dir
SAMPLES=$(ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//' | sort)
cd - > /dev/null

for SAMPLE in $SAMPLES; do
  echo ">>> Processing sample: $SAMPLE"

  if [ "$Paired_end" == "TRUE" ]; then
    STAR --genomeDir $STAR_index_dir \
         --runThreadN $N_thread \
         --readFilesIn ../$Fastq_dir/${SAMPLE}_R1.fastq.gz ../$Fastq_dir/${SAMPLE}_R2.fastq.gz \
         --readFilesCommand zcat \
         --outFileNamePrefix ${SAMPLE}_ \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts

    rsem-calculate-expression --paired-end \
         --alignments --no-bam-output -p $N_thread -q \
         ${SAMPLE}_Aligned.toTranscriptome.out.bam \
         $RSEM_index_name \
         ${SAMPLE}-RSEM

  else
    STAR --genomeDir $STAR_index_dir \
         --runThreadN $N_thread \
         --readFilesIn ../$Fastq_dir/${SAMPLE}_R1.fastq.gz \
         --readFilesCommand zcat \
         --outFileNamePrefix ${SAMPLE}_ \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts

    rsem-calculate-expression \
         --alignments --no-bam-output -p $N_thread -q \
         ${SAMPLE}_Aligned.toTranscriptome.out.bam \
         $RSEM_index_name \
         ${SAMPLE}-RSEM
  fi
done

# ===== [4] Merge expression matrices and generate heatmap =====
echo ">>> Merging expression matrices"

Rscript --vanilla - <<EOF

required_packages <- c("data.table")
installed <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(data.table)

files <- list.files(pattern = "genes.results$")
tpm_list <- list()
count_list <- list()
fpkm_list <- list()

for (file in files) {
  dt <- fread(file)
  sample_id <- sub("-RSEM.*", "", file)
  sample_name <- paste0("sample_", sample_id)

  tpm_dt <- dt[, .(gene_id, temp_value = TPM)]
  setnames(tpm_dt, "temp_value", sample_name)
  tpm_list[[sample_name]] <- tpm_dt

  count_dt <- dt[, .(gene_id, temp_value = expected_count)]
  setnames(count_dt, "temp_value", sample_name)
  count_list[[sample_name]] <- count_dt

  fpkm_dt <- dt[, .(gene_id, temp_value = FPKM)]
  setnames(fpkm_dt, "temp_value", sample_name)
  fpkm_list[[sample_name]] <- fpkm_dt
}

merge_matrix <- function(lst) {
  merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), lst)
  gene_ids <- merged[['gene_id']]
  merged <- merged[, -1]
  merged[is.na(merged)] <- 0
  rownames(merged) <- gene_ids
  return(merged)
}

tpm_matrix <- merge_matrix(tpm_list)
count_matrix <- merge_matrix(count_list)
fpkm_matrix <- merge_matrix(fpkm_list)

fwrite(as.data.table(cbind(gene_id = rownames(count_matrix), count_matrix)),
       file = "Virus-Host_Count_matrix.txt", sep = "\t")
fwrite(as.data.table(cbind(gene_id = rownames(tpm_matrix), tpm_matrix)),
       file = "Virus-Host_TPM_matrix.txt", sep = "\t")
fwrite(as.data.table(cbind(gene_id = rownames(fpkm_matrix), fpkm_matrix)),
       file = "Virus-Host_FPKM_matrix.txt", sep = "\t")

## ===== Merge RSEM isoforms.results to matrices (with gene_id column) =====
library(data.table)


iso_files <- list.files(pattern = "isoforms\\.results$")  
stopifnot(length(iso_files) > 0)


iso_tpm_list   <- list()
iso_count_list <- list()
iso_fpkm_list  <- list()
iso_pct_list   <- list()
tx2gene <- NULL   


for (file in iso_files) {
  dt <- fread(file)


  sample_id   <- sub("-?RSEM.*", "", basename(file))
  sample_name <- paste0("sample_", sample_id)


  need_cols <- c("transcript_id", "gene_id", "TPM", "expected_count", "FPKM", "IsoPct")
  miss <- setdiff(need_cols, names(dt))
  if (length(miss) > 0) stop(sprintf("File %s missing columns: %s", file, paste(miss, collapse = ", ")))


  map_dt  <- unique(dt[, .(transcript_id, gene_id)])
  tx2gene <- if (is.null(tx2gene)) map_dt else unique(rbind(tx2gene, map_dt))


  tpm_dt   <- dt[, .(transcript_id, tmp = TPM)];            setnames(tpm_dt,   "tmp", sample_name)
  cnt_dt   <- dt[, .(transcript_id, tmp = expected_count)]; setnames(cnt_dt,   "tmp", sample_name)
  fpkm_dt  <- dt[, .(transcript_id, tmp = FPKM)];           setnames(fpkm_dt,  "tmp", sample_name)
  pct_dt   <- dt[, .(transcript_id, tmp = IsoPct)];         setnames(pct_dt,   "tmp", sample_name)

  iso_tpm_list[[sample_name]]   <- tpm_dt
  iso_count_list[[sample_name]] <- cnt_dt
  iso_fpkm_list[[sample_name]]  <- fpkm_dt
  iso_pct_list[[sample_name]]   <- pct_dt
}


merge_by_id <- function(lst, id = "transcript_id") {
  stopifnot(length(lst) > 0)

  m <- Reduce(function(x, y) merge(x, y, by = id, all = TRUE), lst)

  num_cols <- setdiff(names(m), id)
  for (cl in num_cols) set(m, which(is.na(m[[cl]])), cl, 0)
  setkeyv(m, id)
  m[]
}

tpm_mat   <- merge_by_id(iso_tpm_list)
count_mat <- merge_by_id(iso_count_list)
fpkm_mat  <- merge_by_id(iso_fpkm_list)
pct_mat   <- merge_by_id(iso_pct_list)


tpm_mat   <- merge(tx2gene, tpm_mat,   by = "transcript_id", all.y = TRUE)
count_mat <- merge(tx2gene, count_mat, by = "transcript_id", all.y = TRUE)
fpkm_mat  <- merge(tx2gene, fpkm_mat,  by = "transcript_id", all.y = TRUE)
pct_mat   <- merge(tx2gene, pct_mat,   by = "transcript_id", all.y = TRUE)


reorder_cols <- function(dt) {
  data.table::setcolorder(dt, c("gene_id","transcript_id",
                                setdiff(names(dt), c("gene_id","transcript_id"))))
  dt
}
tpm_mat   <- reorder_cols(tpm_mat)
count_mat <- reorder_cols(count_mat)
fpkm_mat  <- reorder_cols(fpkm_mat)
pct_mat   <- reorder_cols(pct_mat)


data.table::setorder(tpm_mat,   gene_id, transcript_id)
data.table::setorder(count_mat, gene_id, transcript_id)
data.table::setorder(fpkm_mat,  gene_id, transcript_id)
data.table::setorder(pct_mat,   gene_id, transcript_id)


fwrite(tpm_mat,   "Virus-Host_Isoform_TPM_matrix.txt",   sep = "\t")
fwrite(count_mat, "Virus-Host_Isoform_Count_matrix.txt", sep = "\t")
fwrite(fpkm_mat,  "Virus-Host_Isoform_FPKM_matrix.txt",  sep = "\t")
fwrite(pct_mat,   "Virus-Host_Isoform_IsoPct_matrix.txt",sep = "\t")

EOF
  # === WGCNA (only valid under 'host' mode) ===
  if [ "$SUBMODE" == "WGCNA" ]; then
    echo ">>> WGCNA module enabled"

Rscript TraceVirus_WGCNA.R

echo ">>> Done! WGCNA analysis completed"
  fi

  exit 0
fi

# =======================
# INVALID MODE
# =======================
echo "Error: Unknown mode '$MODE'. Use: bash TraceVirus.sh [virus|host] [WGCNA (optional)]"
exit 1
