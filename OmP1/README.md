# ChIP-seq, ATAC-seq, and RNA-seq Analysis Pipelines

## Description
This comprehensive suite of bioinformatics pipelines is designed for the end-to-end analysis of ChIP-seq, ATAC-seq, and RNA-seq data. It streamlines the entire workflow, from quality control and read preprocessing to alignment, peak calling, visualization, quantification, and differential expression analysis. The pipelines leverage state-of-the-art tools such as FastQC, Trimmomatic, HISAT2, MACS2, deepTools, featureCounts, and MetaSeqR2 to provide a robust and efficient analysis framework.

## Features
- Quality control with FastQC
- Read trimming and filtering with Trimmomatic
- Alignment with HISAT2
- Peak calling with MACS2 (for ChIP-seq and ATAC-seq)
- Visualization with deepTools
- Quantification with featureCounts (for RNA-seq)
- Differential expression analysis with MetaSeqR2 (for RNA-seq)
- Support for both single-end and paired-end data
- Conda environment for easy dependency management
- Nextflow workflow management for scalability and reproducibility

## Prerequisites
- [Nextflow](https://www.nextflow.io/): A workflow tool to run tasks across multiple compute infrastructures in a very portable manner.
- [Conda](https://docs.conda.io/en/latest/): Package, dependency, and environment management for any language—Python, R, Ruby, Lua, Scala, Java, JavaScript, C/ C++, FORTRAN.

## Installation
1. **Install Nextflow**:
   ```bash
   curl -s https://get.nextflow.io | bash
   ```
2. **Install Conda**:
   ```bash
   bash <conda-installer-name>-latest-Linux-x86_64.sh
   ```
3. **Install metaseqR2**:
   ```bash
   if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

   BiocManager::install("metaseqR2")
   ```
## Usage
- **For ChIP-seq**:
  ```bash
  nextflow run chip-seq-pipeline-v1.0.0.SE.nf --profile conda -c merge.config
  ```
- **For ATAC-seq**:
  ```bash
  nextflow run atac-seq-pipeline-v1.0.0.SE.nf --profile conda -c merge.config
  ```
- **For RNA-seq**:
  ```bash
  nextflow run rna-seq-pipeline-v1.2.0.SE.nf --profile conda -c merge.config
  ```
- **For Differential Expression Analysis (RNA-seq)**:
  ```bash
  Rscript metaRUN.R
  ```

Replace `.SE.nf` with `.PE.nf` for paired-end data processing.

## Scripts
- `ReadCount.R`: An R script for counting reads in BAM and FASTQ files and generating summary plots and CSV files.
- `metaRUN.R`: An R script for performing differential gene expression analysis using the MetaSeqR2 package.

## Output
The output files are organized within the `results` directory, with subdirectories corresponding to each analysis step. The differential expression analysis results can be found in the `diff_exp` directory.











