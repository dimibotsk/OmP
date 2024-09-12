# scRNA-seq Analysis Pipeline using Cell Ranger

This repository contains scripts for running Cell Ranger analysis on single-cell RNA sequencing (scRNA-seq) data. The pipeline is implemented in R and uses a YAML configuration file for easy customization.

## Overview

The pipeline automates the process of running Cell Ranger count on multiple scRNA-seq samples. It includes:

- Running Cell Ranger count for each sample
- Parallel processing of samples
- Generation of a summary HTML report

## Requirements

- R (version 3.6 or higher)
- Cell Ranger (must be installed and in your system PATH)
- R packages:
  - yaml
  - parallel

## Files

- `CRange.R`: Main R script for running the pipeline
- `cellranger_config.yaml`: Configuration file for specifying parameters
- `sample_sheet.csv`: CSV file containing sample information

## Usage

1. Clone this repository:
   ```
   git clone https://github.com/dimibotsk/scrnaseq-cellranger-pipeline.git
   cd scrnaseq-cellranger-pipeline
   ```

2. Edit the `cellranger_config.yaml` file to set your specific parameters:
   - Set the paths for your sample sheet, reference transcriptome, and output directory
   - Adjust computational resources as needed

3. Prepare your `sample_sheet.csv` file with the following columns:
   - `sample_id`: Unique identifier for each sample
   - `fastqs`: Path to the directory containing FASTQ files for each sample
   - `condition`: Experimental condition (optional, for downstream analysis)

4. Run the pipeline:
   ```
   Rscript CRange.R
   ```

## Configuration

The `cellranger_config.yaml` file contains the following parameters:

- `sample_sheet`: Path to your sample sheet CSV file
- `reference`: Path to the Cell Ranger reference transcriptome
- `outdir`: Output directory for results
- `cpus`: Number of CPUs to use for each Cell Ranger run
- `memory`: Amount of memory (in GB) to allocate for each Cell Ranger run
- `max_parallel_processes`: Maximum number of samples to process in parallel

## Output

The pipeline generates the following outputs:

- Cell Ranger count results for each sample in separate subdirectories
- An HTML report summarizing the processing of all samples

