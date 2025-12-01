#' Run the TraceVirus computational pipeline
#'
#' @description
#' This function provides a unified R interface for running the external
#' `TraceVirus.sh` pipeline, including viral quantification, host+virus
#' joint quantification, and the optional WGCNA module.  
#' It automatically generates `Parameters.txt`, copies required scripts
#' from the **SingleViral** package, and executes the pipeline inside a
#' Linux login shell using `bash -lc`, ensuring that conda environment
#' variables (PATH, LD_LIBRARY_PATH, etc.) are correctly activated.
#'
#' @param mode Character, `"virus"` or `"host"`.  
#'   Determines which module of the pipeline will be executed.
#'
#' @param run_wgcna Logical.  
#'   If `TRUE`, the WGCNA module will be triggered (only valid when `mode = "host"`).
#'
#' @param n_thread Integer. Number of threads for STAR/RSEM.
#'
#' @param genome_fasta Path to viral FASTA file.
#' @param genome_gtf   Path to viral GTF file.
#'
#' @param host_fasta Path to host FASTA file.
#' @param host_gtf   Path to host GTF file.
#'
#' @param host_virus_fasta Path to merged host+virus FASTA file (used in host mode).
#' @param host_virus_gtf   Path to merged host+virus GTF file (used in host mode).
#'
#' @param output_dir Path to output directory for viral-only analysis.
#' @param output_dir_host Path to output directory for host+virus analysis.
#'
#' @param rsem_index_name Prefix name for RSEM index.
#' @param star_index_dir Directory name for STAR index.
#'
#' @param fastq_dir Directory containing FASTQ files.
#'
#' @param paired_end Logical. Whether sequencing is paired-end (`TRUE`) or single-end (`FALSE`).
#'
#' @param traits_file Path to traits CSV file (required only when `run_wgcna = TRUE`).
#'
#' @param traits_basename Basename of the copied traits file (default `"traits.csv"`).
#'
#' @param work_dir Working directory where all pipeline files will be stored.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `work_dir` — the working directory used
#'     \item `exit_status` — exit status returned by `bash`
#'   }
#'
#' @details
#' The function wraps the external `TraceVirus.sh` pipeline into an easy-to-use
#' R interface. Because `system2("bash")` does **not** load conda environment variables,
#' the pipeline may fail to find `STAR` or `RSEM`.  
#'
#' To fix this, the function uses:
#'
#' ```
#' bash -lc "cd work_dir && bash TraceVirus.sh ..."
#' ```
#'
#' `bash -lc` launches a *login shell*, which reads:
#' - `.bashrc`
#' - `.bash_profile`
#' - Conda environment activation
#'
#' ensuring all tools are found.
#'
#' @examples
#' \dontrun{
#' RunTraceVirus(
#'   mode = "virus",
#'   fastq_dir = "./fastq/"
#' )
#' }
#'
#' @export
RunTraceVirus <- function(
    mode              = c("virus", "host"),
    run_wgcna         = FALSE,
    n_thread          = 20,
    genome_fasta      = "hhv6.fa",
    genome_gtf        = "hhv6.gtf",
    host_fasta        = "hg38.fa",
    host_gtf          = "hg38.gtf",
    host_virus_fasta  = "host_virus.fa",
    host_virus_gtf    = "host_virus.gtf",
    output_dir        = "./Virus_Track_Output",
    output_dir_host   = "./Host_Virus_Output",
    rsem_index_name   = "hhv6_RSEM",
    star_index_dir    = "hhv6_star_index",
    fastq_dir         = "",
    paired_end        = FALSE,
    traits_file       = NULL,
    traits_basename   = "traits.csv",
    work_dir          = getwd()
){
  mode <- match.arg(mode)
  
  ## -------- Create working directory --------
  if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
  
  ## -------- Write Parameters.txt --------
  param_lines <- c(
    sprintf("N_thread=%d", n_thread),
    sprintf('Genome_fasta="%s"', genome_fasta),
    sprintf('Genome_gtf="%s"', genome_gtf),
    sprintf('host_fasta="%s"', host_fasta),
    sprintf('host_gtf="%s"', host_gtf),
    sprintf('host_virus_fasta="%s"', host_virus_fasta),
    sprintf('host_virus_gtf="%s"', host_virus_gtf),
    sprintf('Output_directory="%s"', output_dir),
    sprintf('Output_directory_host="%s"', output_dir_host),
    sprintf('RSEM_index_name="%s"', rsem_index_name),
    sprintf('STAR_index_dir="%s"', star_index_dir),
    sprintf('Fastq_dir="%s"', fastq_dir),
    sprintf('Paired_end=%s', ifelse(paired_end, "TRUE", "FALSE"))
  )
  
  ## -------- Copy traits.csv if needed --------
  if (run_wgcna) {
    if (is.null(traits_file) || !file.exists(traits_file)) {
      stop("run_wgcna = TRUE but no valid 'traits_file' was provided.")
    }
    
    file.copy(traits_file, file.path(work_dir, traits_basename), overwrite = TRUE)
    param_lines <- c(param_lines, sprintf('Traits_file="%s"', traits_basename))
  }
  
  writeLines(param_lines, file.path(work_dir, "Parameters.txt"))
  
  ## -------- Copy internal scripts --------
  script_sh <- system.file("TraceVirus", "TraceVirus.sh", package = "SingleViral")
  script_R  <- system.file("TraceVirus", "TraceVirus_WGCNA.R", package = "SingleViral")
  
  file.copy(script_sh, file.path(work_dir, "TraceVirus.sh"), overwrite = TRUE)
  file.copy(script_R, file.path(work_dir, "TraceVirus_WGCNA.R"), overwrite = TRUE)
  
  ## -------- Run in a login shell (critical for conda) --------
  cmd <- sprintf(
    'cd %s && bash TraceVirus.sh %s %s',
    shQuote(normalizePath(work_dir)),
    mode,
    ifelse(run_wgcna, "WGCNA", "")
  )
  
  message(">>> Running command in login shell:\n", cmd)
  
  status <- system(
    sprintf('bash -lc %s', shQuote(cmd)),
    intern = FALSE,
    wait = TRUE
  )
  
  invisible(list(
    work_dir     = work_dir,
    exit_status  = status
  ))
}