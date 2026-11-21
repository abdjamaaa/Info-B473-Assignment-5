# INFO-B473 
# Assignment 5: Advanced R 
# Abdul Djama
# 11/21/26


# Load the packages needed for this assignment
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)

# Make a folder to save the results
dir.create("results", showWarnings = FALSE)

# Part 1 Working with the gene data

# Read in the three files
gene_data <- read_excel("Gene_Expression_Data.xlsx")
gene_info <- read.csv("Gene_Information.csv")
sample_info <- read.delim("Sample_Information.tsv",
                          skip = 1, header = FALSE)

# Give the sample information columns names
colnames(sample_info) <- c("Sample_ID", "Phenotype", "Patient")

# Make sure the sample types are written the same way
sample_info$Phenotype <- tolower(sample_info$Phenotype)

# Get the sample names from the expression data
sample_names <- colnames(gene_data)[-1]

# Match the sample information with the expression data
sample_info <- sample_info[
  match(sample_names, sample_info$Sample_ID), ]

# Rename each sample to include whether it is tumor or normal
new_names <- paste0(sample_names, "_", sample_info$Phenotype)
colnames(gene_data)[-1] <- new_names

write.csv(gene_data, "results/01_expression_renamed.csv",
          row.names = FALSE)


# Part 2: Separate tumor and normal samples

tumor_cols <- new_names[sample_info$Phenotype == "tumor"]
normal_cols <- new_names[sample_info$Phenotype == "normal"]

# Create separate datasets for both groups
tumor_data <- gene_data[, c("Probe_ID", tumor_cols)]
normal_data <- gene_data[, c("Probe_ID", normal_cols)]

write.csv(tumor_data, "results/02_tumor_data.csv",
          row.names = FALSE)

write.csv(normal_data, "results/03_normal_data.csv",
          row.names = FALSE)


# Part 3: Average gene expression

# Find the average expression for every probe
tumor_avg <- rowMeans(tumor_data[, -1], na.rm = TRUE)
normal_avg <- rowMeans(normal_data[, -1], na.rm = TRUE)

# Put the averages into one table
results <- data.frame(
  Probe_ID = gene_data$Probe_ID,
  Tumor_Average = tumor_avg,
  Normal_Average = normal_avg
)

# Calculate the log2 fold change
# The regular log2 ratio lets us see which group has higher expression
results$Log2FC <- ifelse(
  tumor_avg > 0 & normal_avg > 0,
  log2(tumor_avg / normal_avg),
  NA
)

# Also calculates the formula given in the assignment
# This formula is not defined for every gene
results$Assignment_FC <- ifelse(
  tumor_avg > normal_avg & normal_avg > 0,
  log2((tumor_avg - normal_avg) / normal_avg),
  NA
)

write.csv(results, "results/04_fold_changes.csv",
          row.names = FALSE)


# Part 4: Find the differentially expressed genes

# Add the gene information to our results
results <- left_join(results, gene_info, by = "Probe_ID")

# Keep probes with an absolute log2 fold change above 5
degs <- results %>%
  filter(is.finite(Log2FC), abs(Log2FC) > 5)

# Label which sample group has higher expression
degs$Higher_In <- ifelse(
  degs$Log2FC > 0, "Tumor", "Normal"
)

# Save the list of selected genes
write.csv(degs, "results/05_selected_degs.csv",
          row.names = FALSE)

# Print how many probes passed the cutoff
cat("Number of selected probes:", nrow(degs), "\n")


# Part 5: Explore the results with graphs

# Graph 1: Number of selected genes on each chromosome
graph1 <- ggplot(degs, aes(x = as.factor(Chromosome))) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "Differentially Expressed Genes by Chromosome",
    x = "Chromosome",
    y = "Number of Genes"
  ) +
  theme_minimal()

ggsave("results/06_chromosome_counts.png",
       graph1, width = 10, height = 6)


# Graph 2: Chromosomes separated by tumor and normal
graph2 <- ggplot(degs, aes(x = as.factor(Chromosome),
                           fill = Higher_In)) +
  geom_bar(position = "dodge") +
  labs(
    title = "Gene Expression by Chromosome and Sample Type",
    x = "Chromosome",
    y = "Number of Genes",
    fill = "Higher Expression"
  ) +
  theme_minimal()

ggsave("results/07_chromosome_groups.png",
       graph2, width = 10, height = 6)


# Graph 3: Percentage of genes higher in tumor or normal
percent_data <- degs %>%
  count(Higher_In) %>%
  mutate(Percentage = n / sum(n) * 100)

graph3 <- ggplot(percent_data,
                 aes(x = Higher_In, y = Percentage,
                     fill = Higher_In)) +
  geom_col() +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            vjust = -0.5) +
  labs(
    title = "Percentage of Differentially Expressed Genes",
    x = "Higher Expression Group",
    y = "Percentage"
  ) +
  ylim(0, 110) +
  theme_minimal()

ggsave("results/08_expression_percentages.png",
       graph3, width = 8, height = 6)


# Part 6: Heatmap and Clustermap

raw_data <- read_excel("Gene_Expression_Data.xlsx")

# Turn the expression values into a matrix
expression_matrix <- as.matrix(raw_data[, -1])
rownames(expression_matrix) <- raw_data$Probe_ID

# Log transform the expression values
log_expression <- log2(expression_matrix + 1)

# Find the genes with the most variation
gene_variation <- apply(log_expression, 1, var, na.rm = TRUE)
gene_variation[!is.finite(gene_variation)] <- -Inf

# Use the 100 genes with the greatest variation
top_genes <- order(gene_variation, decreasing = TRUE)[1:100]
heatmap_data <- log_expression[top_genes, ]

# Replace any missing values with the row average
for (i in 1:nrow(heatmap_data)) {
  missing <- is.na(heatmap_data[i, ])
  if (any(missing)) {
    heatmap_data[i, missing] <- mean(
      heatmap_data[i, !missing], na.rm = TRUE
    )
  }
}

# Scale the data to make differences easier to see
heatmap_data <- t(scale(t(heatmap_data)))
heatmap_data[!is.finite(heatmap_data)] <- 0

# Add the sample types to the heatmap
colnames(heatmap_data) <- new_names

sample_labels <- data.frame(
  Phenotype = sample_info$Phenotype,
  Patient = sample_info$Patient
)

rownames(sample_labels) <- new_names


# Heatmap without clustering
png("results/09_heatmap.png", width = 1400, height = 1100)

pheatmap(
  heatmap_data,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = sample_labels,
  show_rownames = FALSE,
  main = "Gene Expression Heatmap"
)

dev.off()


# Clustermap to compare similarities between samples
png("results/10_clustermap.png", width = 1400, height = 1100)

pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = sample_labels,
  show_rownames = FALSE,
  main = "Clustered Gene Expression Heatmap"
)

dev.off()

# Summarized Results

tumor_count <- sum(degs$Higher_In == "Tumor")
normal_count <- sum(degs$Higher_In == "Normal")

cat("Total probes:", nrow(gene_data), "\n")
cat("Differentially expressed probes:", nrow(degs), "\n")
cat("Higher in tumor:", tumor_count, "\n")
cat("Higher in normal:", normal_count, "\n")

summary_text <- c(
  paste("Total probes:", nrow(gene_data)),
  paste("Differentially expressed probes:", nrow(degs)),
  paste("Higher in tumor:", tumor_count),
  paste("Higher in normal:", normal_count),
  "The chromosome graphs show where the differentially expressed genes are found.",
  "The heatmap and clustermap compare gene expression across the samples."
)

writeLines(summary_text, "results/11_summary.txt")

