# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(DoubletFinder)
library(tidyr)
library(tibble)
library(SCP)


# Read configuration file
config <- read.table("config.txt", sep = "=", row.names = 1, col.names = c("name", "path"), stringsAsFactors = FALSE)

# Create directories for pre and post QC plots
dir.create("pre_QC", showWarnings = FALSE)
dir.create("post_QC", showWarnings = FALSE)

# Function to read data and create Seurat object
create_seurat_object <- function(path, project_name) {
  data <- Read10X(data.dir = path)
  seurat_obj <- CreateSeuratObject(counts = data, project = project_name, min.cells = 3, min.features = 200)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  return(seurat_obj)
}

# Create Seurat objects
wt_rep1 <- create_seurat_object(config["WT_REP1", "path"], "WT_REP1")
wt_rep2 <- create_seurat_object(config["WT_REP2", "path"], "WT_REP2")
ko_rep1 <- create_seurat_object(config["KO_REP1", "path"], "KO_REP1")
ko_rep2 <- create_seurat_object(config["KO_REP2", "path"], "KO_REP2")

# Function to plot QC metrics
plot_qc_metrics <- function(seurat_obj, output_prefix, output_dir) {
  png(file.path(output_dir, paste0(output_prefix, "_qc_plots.png")), width = 500, height = 400)
  print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()
  
  png(file.path(output_dir, paste0(output_prefix, "_feature_scatter.png")), width = 700, height = 400)
  plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2)
  dev.off()
}

# Plot QC metrics for each sample (pre-QC)
plot_qc_metrics(wt_rep1, "WT_REP1", "pre_QC")
plot_qc_metrics(wt_rep2, "WT_REP2", "pre_QC")
plot_qc_metrics(ko_rep1, "KO_REP1", "pre_QC")
plot_qc_metrics(ko_rep2, "KO_REP2", "pre_QC")

# Function to filter cells
filter_cells <- function(seurat_obj) {
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5 & nCount_RNA < 18000)
  return(seurat_obj)
}

# Filter cells for each sample
wt_rep1_filtered <- filter_cells(wt_rep1)
wt_rep2_filtered <- filter_cells(wt_rep2)
ko_rep1_filtered <- filter_cells(ko_rep1)
ko_rep2_filtered <- filter_cells(ko_rep2)

# Function to generate summary
generate_summary <- function(original, filtered, sample_name) {
  summary <- paste0("Summary for ", sample_name, ":\n",
                    "Original cell count: ", ncol(original), "\n",
                    "Filtered cell count: ", ncol(filtered), "\n",
                    "Cells removed: ", ncol(original) - ncol(filtered), "\n\n")
  return(summary)
}

# Generate and save summary for all samples
summary_text <- c(
  generate_summary(wt_rep1, wt_rep1_filtered, "WT_REP1"),
  generate_summary(wt_rep2, wt_rep2_filtered, "WT_REP2"),
  generate_summary(ko_rep1, ko_rep1_filtered, "KO_REP1"),
  generate_summary(ko_rep2, ko_rep2_filtered, "KO_REP2")
)

writeLines(summary_text, "filtering_summary.txt")
file.rename("filtering_summary.txt", "post_QC/filtering_summary.txt")

# Save filtered Seurat objects
#saveRDS(wt_rep1_filtered, "WT_REP1_filtered.rds")
#saveRDS(wt_rep2_filtered, "WT_REP2_filtered.rds")
#saveRDS(ko_rep1_filtered, "KO_REP1_filtered.rds")
#saveRDS(ko_rep2_filtered, "KO_REP2_filtered.rds")

# Plot QC metrics for each filtered sample (post-QC)
plot_qc_metrics(wt_rep1_filtered, "WT_REP1_filtered", "post_QC")
plot_qc_metrics(wt_rep2_filtered, "WT_REP2_filtered", "post_QC")
plot_qc_metrics(ko_rep1_filtered, "KO_REP1_filtered", "post_QC")
plot_qc_metrics(ko_rep2_filtered, "KO_REP2_filtered", "post_QC")

########################################################################################################################################################

########################################################################################################################################################

########################################################################################################################################################


# Create a list of Seurat objects
seurat_list <- list(WT_REP1 = wt_rep1_filtered, 
                    WT_REP2 = wt_rep2_filtered, 
                    KO_REP1 = ko_rep1_filtered, 
                    KO_REP2 = ko_rep2_filtered)

# Normalize and identify variable features for each dataset independently
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seurat_list)

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)

# Integrate data
integrated_seurat <- IntegrateData(anchorset = anchors)

# Switch to integrated assay for downstream analysis
DefaultAssay(integrated_seurat) <- "integrated"

# Scale the data
all_genes <- rownames(integrated_seurat)
integrated_seurat <- ScaleData(integrated_seurat, features = all_genes)

# Perform linear dimensional reduction
integrated_seurat <- RunPCA(integrated_seurat, features = VariableFeatures(object = integrated_seurat))

# Determine the dimensionality of the dataset
ElbowPlot(integrated_seurat)

# Cluster the cells
integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:15)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.4)

# Run non-linear dimensional reduction (UMAP)
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:10)

# Add condition information
integrated_seurat$condition <- ifelse(grepl("^WT", integrated_seurat$orig.ident), "WT", "KO")

# Plot UMAP
umap_plot <- DimPlot(integrated_seurat, reduction = "umap", 
                     group.by = c("orig.ident", "condition", "seurat_clusters"), 
                     combine = FALSE)
ggsave("umap_plot_integrated.png", plot = wrap_plots(umap_plot, ncol = 3), width = 18, height = 6)

umap_plot <- DimPlot(integrated_seurat, reduction = "umap", label = TRUE)
ggsave("umap_plot_integrated_labeled.png")

dir.create("umaps", showWarnings = FALSE)
file.rename("umap_plot_integrated.png", "umaps/umap_plot_integrated.png")
file.rename("umap_plot_integrated_labeled.png", "umaps/umap_plot_integrated_labeled.png")


# Save the integrated and processed Seurat object
saveRDS(integrated_seurat, file = "integrated_processed_seurat.rds")

# Generate a summary of the integrated data
integrated_summary <- paste0(
  "Integrated data summary:\n",
  "Total cells: ", ncol(integrated_seurat), "\n",
  "Total features: ", nrow(integrated_seurat), "\n",
  "Cells per sample:\n",
  table(integrated_seurat$orig.ident) %>% capture.output() %>% paste(collapse = "\n"), "\n",
  "Cells per condition:\n",
  table(integrated_seurat$condition) %>% capture.output() %>% paste(collapse = "\n"), "\n",
  "Cells per cluster:\n",
  table(Idents(integrated_seurat)) %>% capture.output() %>% paste(collapse = "\n"), "\n"
)

cat(integrated_summary, file = "integrated_data_summary.txt")

# Plot QC metrics for integrated data
plot_qc_metrics(integrated_seurat, "integrated_seurat", "post_QC")


########################################################################################################################################################

########################################################################################################################################################

########################################################################################################################################################

# Normalized Clusters Percentages between conditions

# Update the get_cluster_counts function
get_cluster_counts <- function(seurat_obj) {
  cluster_counts <- table(Idents(seurat_obj), seurat_obj$condition)
  cluster_counts_df <- as.data.frame(cluster_counts)
  colnames(cluster_counts_df) <- c("Cluster", "Condition", "Count")
  
  # Calculate percentages within each condition
  cluster_counts_df <- cluster_counts_df %>%
    group_by(Condition) %>%
    mutate(Percentage = Count / sum(Count) * 100) %>%
    ungroup()
  
  # Reshape the data for easier comparison
  cluster_counts_wide <- cluster_counts_df %>%
    select(Cluster, Condition, Percentage) %>%
    pivot_wider(names_from = Condition, values_from = Percentage, names_prefix = "Percent_")
  
  # Add the raw counts
  cluster_counts_wide <- cluster_counts_df %>%
    select(Cluster, Condition, Count) %>%
    pivot_wider(names_from = Condition, values_from = Count, names_prefix = "Count_") %>%
    right_join(cluster_counts_wide, by = "Cluster")
  
  return(cluster_counts_wide)
}

cluster_stats <- get_cluster_counts(integrated_seurat)

# Create a summary string
cluster_summary <- "Cluster statistics:\n\n"
cluster_summary <- paste0(cluster_summary, "Cluster\tWT Count\tWT %\tKO Count\tKO %\n")
for (i in 1:nrow(cluster_stats)) {
  cluster_summary <- paste0(cluster_summary,
                            cluster_stats$Cluster[i], "\t",
                            cluster_stats$Count_WT[i], "\t",
                            round(cluster_stats$Percent_WT[i], 2), "%\t",
                            cluster_stats$Count_KO[i], "\t",
                            round(cluster_stats$Percent_KO[i], 2), "%\n")
}

# Write the summary to a text file
writeLines(cluster_summary, "cluster_statistics_summary.tsv")

# Create a bar plot for comparison
cluster_stats_long <- cluster_stats %>%
  tidyr::pivot_longer(cols = tidyselect::starts_with("Percent_"), 
                      names_to = "Condition", 
                      values_to = "Percentage") %>%
  mutate(Condition = gsub("Percent_", "", Condition)) %>%
  mutate(Condition = factor(Condition, levels = c("WT", "KO")))  # Setting order of conditions

comparison_plot <- ggplot(cluster_stats_long, aes(x = Cluster, y = Percentage, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  scale_fill_manual(values = c("WT" = "deepskyblue4", "KO" = "darkorange")) +  # Setting custom colors
  labs(title = "Cluster Percentages by Condition",
       y = "Percentage of Cells",
       x = "Cluster")

# Save the comparison plot
ggsave("cluster_percentages_comparison.png", comparison_plot, width = 12, height = 8, bg = "white")


# Create directories
dir.create("integrated", showWarnings = FALSE)
file.rename("cluster_statistics_summary.tsv", "integrated/cluster_statistics_summary.tsv")
file.rename("cluster_percentages_comparison.png", "integrated/cluster_percentages_comparison.png")

# split object to wt and ko
# Create WT Seurat object
wt_cells <- WhichCells(integrated_seurat, expression = condition == "WT")
WT <- subset(integrated_seurat, cells = wt_cells)

# Create KO Seurat object
ko_cells <- WhichCells(integrated_seurat, expression = condition == "KO")
KO <- subset(integrated_seurat, cells = ko_cells)

########################################################################################################################################################

########################################################################################################################################################

########################################################################################################################################################

# Differential Expression Analysis between KO and WT

# Set the default assay to "RNA" for DE analysis
DefaultAssay(integrated_seurat) <- "RNA"

# Join layers for the integrated assay if necessary
integrated_seurat <- JoinLayers(integrated_seurat)

# Create a directory for DE results if it doesn't exist
dir.create("DE_results", showWarnings = FALSE)

# Perform differential expression analysis
# 1. Genes upregulated in KO compared to WT
KO_EXPRESSED_MARKERS <- FindMarkers(integrated_seurat, 
                                    ident.1 = "KO", 
                                    ident.2 = "WT", 
                                    group.by = "condition",
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.25,
                                    only.pos = TRUE)

# Add gene names to the results
KO_EXPRESSED_MARKERS$gene <- rownames(KO_EXPRESSED_MARKERS)

# Save KO upregulated genes as TSV
write.table(KO_EXPRESSED_MARKERS, "DE_results/KO_upregulated_genes.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# 2. Genes upregulated in WT compared to KO
WT_EXPRESSED_MARKERS <- FindMarkers(integrated_seurat, 
                                    ident.1 = "WT", 
                                    ident.2 = "KO", 
                                    group.by = "condition",
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.25,
                                    only.pos = TRUE)

# Add gene names to the results
WT_EXPRESSED_MARKERS$gene <- rownames(WT_EXPRESSED_MARKERS)

# Save WT upregulated genes as TSV
write.table(WT_EXPRESSED_MARKERS, "DE_results/WT_upregulated_genes.tsv", sep="\t", quote=FALSE, row.names=FALSE)



########################################################################################################################################################

########################################################################################################################################################

########################################################################################################################################################

# Annotate clusters based on marker genes

# Function to read cell type markers from markers.txt file
read_cell_type_markers <- function(file_path = "markers.txt") {
  if (!file.exists(file_path)) {
    stop("markers.txt file not found. Please ensure it's in the current working directory.")
  }
  
  lines <- readLines(file_path)
  lines <- lines[!grepl("^#", lines)]  # Remove comment lines
  markers <- list()
  for (line in lines) {
    if (grepl("=", line)) {
      parts <- strsplit(line, "=")[[1]]
      cell_type <- trimws(parts[1])
      genes <- strsplit(trimws(parts[2]), ",")[[1]]
      markers[[cell_type]] <- trimws(genes)
    }
  }
  return(markers)
}


# Read markers
markers <- read_cell_type_markers()

# Flatten the markers list into a data frame
markers_df <- do.call(rbind, lapply(names(markers), function(cell_type) {
  data.frame(cell_type = cell_type, gene = markers[[cell_type]])
}))

# Create annotation_results directory
dir.create("annotation_results", showWarnings = FALSE)

# Function to create improved dotplots
create_dotplots <- function(integrated_seurat, condition) {
  DefaultAssay(integrated_seurat) <- "RNA"
  
  # Create dotplots with improved aesthetics
  for (cell_type in unique(markers_df$cell_type)) {
    genes <- markers_df$gene[markers_df$cell_type == cell_type]
    genes <- genes[genes %in% rownames(integrated_seurat)]
    
    if (length(genes) > 0) {
      p <- DotPlot(integrated_seurat, features = genes, group.by = "seurat_clusters", 
                   cols = c("deepskyblue", "darkorange"), dot.scale = 8) +
        ggtitle(paste("Marker Expression for", cell_type, "in", condition)) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(fill = "white", colour = NA)
        ) +
        guides(size = guide_legend(title = "Percent Expressed"),
               color = guide_colorbar(title = "Average Expression"))
      
      filename <- paste0("dotplot_", make.names(cell_type), "_", condition, ".png")
      filepath <- file.path("annotation_results", filename)
      ggsave(filepath, plot = p, bg = "white")
    }
  }
}

# Process WT condition
wt_cells <- WhichCells(integrated_seurat, expression = condition == "WT")
wt_seurat <- subset(integrated_seurat, cells = wt_cells)
create_dotplots(wt_seurat, "WT")

# Process KO condition
ko_cells <- WhichCells(integrated_seurat, expression = condition == "KO")
ko_seurat <- subset(integrated_seurat, cells = ko_cells)
create_dotplots(ko_seurat, "KO")

print("Dotplot generation completed for WT and KO separately. Improved dotplots have been saved in the 'annotation_results' directory.")
