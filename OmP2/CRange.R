#!/usr/bin/env Rscript

# Load required libraries
library(yaml)
library(parallel)

# Read configuration
config <- yaml::read_yaml("cellranger_config.yaml")

# Create output directory if it doesn't exist
dir.create(config$outdir, showWarnings = FALSE, recursive = TRUE)

# Function to run Cell Ranger count
run_cellranger_count <- function(sample) {
  sample_id <- sample$sample_id
  fastqs <- sample$fastqs
  
  # Create a subdirectory for this sample
  sample_outdir <- file.path(config$outdir, sample_id)
  dir.create(sample_outdir, showWarnings = FALSE)
  
  # Change to the sample output directory
  setwd(sample_outdir)
  
  cmd <- sprintf(
    "cellranger count --id=%s --transcriptome=%s --fastqs=%s --sample=%s --localcores=%d --localmem=%d",
    sample_id, config$reference, fastqs, sample_id, config$cpus, config$memory
  )
  
  # Run the command and capture output
  output <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
  
  # Change back to the original directory
  setwd(config$original_dir)
  
  # Return the output for logging
  return(list(sample_id = sample_id, output = output))
}

# Save the original directory
config$original_dir <- getwd()

# Read sample sheet
samples <- read.csv(config$sample_sheet, stringsAsFactors = FALSE)

# Run Cell Ranger count for each sample
results <- mclapply(
  split(samples, seq(nrow(samples))),
  run_cellranger_count,
  mc.cores = min(nrow(samples), config$max_parallel_processes)
)

# Generate HTML report
report <- c(
  "<html><body>",
  "<h1>Cell Ranger Count Report</h1>",
  sprintf("<p>Processed %d samples</p>", nrow(samples)),
  sprintf("<p>Reference: %s</p>", config$reference),
  "<ul>"
)

for (result in results) {
  report <- c(report, sprintf("<li>%s: %s</li>", result$sample_id, 
              ifelse(length(result$output) > 0, "Completed", "Failed")))
  
  # Add detailed output to report
  report <- c(report, "<pre>", result$output, "</pre>")
}

report <- c(report, "</ul></body></html>")

writeLines(report, file.path(config$outdir, "cellranger_report.html"))

cat("Cell Ranger count processes completed. Report generated at", 
    file.path(config$outdir, "cellranger_report.html"), "\n")
