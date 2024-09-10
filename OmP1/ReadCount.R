#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyr)

# Function to count reads in BAM files
count_bam_reads <- function(bam_file) {
  cmd <- paste("samtools view", bam_file, "| wc -l")
  as.numeric(system(cmd, intern = TRUE))
}

# Function to count reads in FASTQ files
count_fastq_reads <- function(fastq_file) {
  cmd <- paste("zcat", fastq_file, "| wc -l")
  as.numeric(system(cmd, intern = TRUE)) / 4  # Divide by 4 because each read in FASTQ is 4 lines
}

# Get list of filtered BAM files
bam_files <- list.files("results/filtered_bam", pattern = "\\.filtered\\.bam$", full.names = TRUE)

# Get list of FASTQ files
fastq_files <- list.files(".", pattern = "\\.fastq\\.gz$", full.names = TRUE)

# Count reads in BAM files
bam_counts <- sapply(bam_files, count_bam_reads)
bam_df <- data.frame(
  Sample = gsub("\\.filtered\\.bam$", "", basename(bam_files)),
  Reads = bam_counts,
  Type = "Filtered BAM"
)

# Count reads in FASTQ files
fastq_counts <- sapply(fastq_files, count_fastq_reads)
fastq_df <- data.frame(
  Sample = gsub("\\.fastq\\.gz$", "", basename(fastq_files)),
  Reads = fastq_counts,
  Type = "Original FASTQ"
)

# Combine data
all_counts <- rbind(bam_df, fastq_df)

# Reshape data for stacked bar plot
plot_data <- all_counts %>%
  pivot_wider(names_from = Type, values_from = Reads) %>%
  mutate(Unmapped = `Original FASTQ` - `Filtered BAM`) %>%
  pivot_longer(cols = c("Filtered BAM", "Unmapped"), names_to = "Category", values_to = "Reads")

# Create stacked bar plot
p <- ggplot(plot_data, aes(x = Sample, y = Reads, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Read Counts: Original vs Filtered",
       x = "Sample",
       y = "Number of Reads") +
  scale_fill_manual(values = c("Filtered BAM" = "blue", "Unmapped" = "red")) +
  scale_y_continuous(labels = scales::comma)

# Save the plot
ggsave("read_count_comparison.png", p, width = 12, height = 8)

# Print summary to console
print(all_counts)

# Save counts to CSV
write.csv(all_counts, "read_counts.csv", row.names = FALSE)
