<p align="center">
  <img width="300"  src="https://github.com/YuZhang-loki/SingleViral/blob/main/SingleViral_logo.png">
</p>

# SingleViral
## Version 0.1.0

An integrated R toolkit for viral infection quantification, classification,
and host–virus interaction analysis at the single-cell level.

## 🧬 Overview

**SingleViral** provides a complete, modular workflow to measure viral load,
define infection levels, identify infection-related genes, and explore
infection-driven functional and communication changes across single cells.

<p align="center">
  <img width="800"  src="https://github.com/YuZhang-loki/SingleViral/blob/main/SingleViral_overview.png">
</p>

### Key features
- preprocessing for multiple single-cell RNA-seq samples (`preanalysis`)
- Viral load quantification (`ViralloadScoring`)
- Viral gene dynamics analysis (`VirusGeneDynamic`)
- Infection-level classification (`TargetInfectionLevel`)
- DEG analysis and volcano plotting (`InfectionVariableGenes`)
- GSEA-based dynamic analysis (`InfectionDynamic`)
- Group-wise infection summary (`CellVariableInfection`)
- Cell-cycle scoring by infection (`CellCycleByInfection`)
- CellChat analysis by infection (`CellChatByCellInfection`)
- Compare Cell–Cell Communication (`CellChatCompare`)
- Run TraceVirus Pipeline (`RunTraceVirus`)

### Usage

A detailed usage workflow is provided in the scripts folder.

## 🚀 Installation

You can install the latest development version directly from GitHub:

```r
install.packages("devtools")
devtools::install_github("YuZhang-loki/SingleViral")
```

### Installation of other dependencies

- Install [CellChat](https://github.com/sqjin/CellChat) using `devtools::install_github("jinworks/CellChat")` if you encounter any issue.
- Install [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) using `devtools::install_github("jokergoo/ComplexHeatmap")` if you encounter any issue.

# TraceVirus
TraceVirus is an R-based computational framework built upon STAR and RSEM, designed to detect viral and host sequences from raw RNA-sequencing (scRNA-seq) data and to automate downstream analyses. The source code is available in the inst directory.

Before running TraceVirus, several dependencies must be installed :
1 . The first step is to install R software. Once this is done, several Bioconductor packages have to be installed too. To do so start a R session and type :

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi","GO.db", "preprocessCore", "impute","WGCNA"))
```

2.Spliced Transcript Alignment to A Reference (STAR) has to be installed. The full installation process is described on the STAR Github(https://github.com/alexdobin/STAR).

3.RSEM is also required. The installation is described here(https://github.com/deweylab/RSEM).

4.Place the data to be analyzed into the working directory.
For single-end sequencing data, rename the sample as ***_R1.fastq.gz (e.g., 101_R1.fastq.gz).
For paired-end sequencing data, rename the files as ***_R1.fastq.gz and ***_R2.fastq.gz (e.g., 101_R1.fastq.gz and 101_R2.fastq.gz).
After renaming, run the following function.
