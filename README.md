<p align="center">
  <img width="300"  src="https://github.com/YuZhang-loki/SingleViral/blob/main/SingleViral_logo.png">
</p>

# SingleViral
## Version 1.0.0

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
- Identification and visualization of dominant viral genes (`PlotViralTopGene`)
- Infection-level classification (`TargetInfectionLevel`)
- DEG analysis and volcano plotting (`InfectionVariableGenes`)
- GSEA-based dynamic analysis (`InfectionDynamic`)
- Group-wise infection summary (`CellVariableInfection`)
- Patient-paired pseudobulk analysis (`InfectionPseudobulk`)
- Compute correlation between gene expression and viral load (`GeneViralLoadCorrelation`)
- Infection-level progression clustering for viral and host genes (`InfectionSeries`)
- Cell-cycle scoring by infection (`CellCycleByInfection`)
- CellChat analysis by infection (`CellChatByCellInfection`)
- Compare Cell–Cell Communication (`CellChatCompare`)
- Run TraceVirus Pipeline (`RunTraceVirus`)
### Usage

A detailed usage workflow is provided in the usage folder.

## 🚀 Installation

You can install the latest development version directly from GitHub:

```r
install.packages("devtools")
devtools::install_github("YuZhang-loki/SingleViral")
```

### Installation of other dependencies

- Install [CellChat](https://github.com/sqjin/CellChat) using `devtools::install_github("jinworks/CellChat")` if you encounter any issue.
- Install [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) using `devtools::install_github("jokergoo/ComplexHeatmap")` if you encounter any issue.
For single-end sequencing data, rename the sample as ***_R1.fastq.gz (e.g., 101_R1.fastq.gz).
For paired-end sequencing data, rename the files as ***_R1.fastq.gz and ***_R2.fastq.gz (e.g., 101_R1.fastq.gz and 101_R2.fastq.gz).
After renaming, run the following function.
