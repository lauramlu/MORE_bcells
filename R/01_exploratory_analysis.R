################################################################################
# Exploratory Data Analysis (EDA)
# Author: Laura Muñoz
# Date: 01/02/2026
# Goal:
# Perform an initial exploratory analysis of RNA expression data and sample
# metadata, including consistency checks, sample composition summaries and PCA.
#
# Input:
# - 00_data/processed/rna_norm_data.rds
# - 00_data/raw/rna_metadata.rds
#
# Output:
# - PCA plots
# - Expression density plot
# - Sample composition summaries
# - sessionInfo.txt
################################################################################

#-------------------------------------------------------------------------------
# 0. Set up
#-------------------------------------------------------------------------------
rm(list = ls())
gc()

# print R version and working directory
message("R version: ", R.version.string)
message("Working directory: ", getwd())

# Packages
library(ggplot2) ; library(dplyr) ; library(tidyr)

in_dir <- "00_data"
out_dir <- "02_results/01_exploratory"

#-------------------------------------------------------------------------------
# 1. Read data
#-------------------------------------------------------------------------------
expr <- readRDS(file.path(in_dir, "processed", "rna_norm_data.rds"))
metadata <- readRDS(file.path(in_dir, "raw", "rna_metadata.rds"))

#-------------------------------------------------------------------------------
# 2. Object type and dimensions
#-------------------------------------------------------------------------------
message("expr class:", class(expr), "\n")
message("metadata class:", class(metadata), "\n")

message("expr dim:", nrow(expr), "genes x", ncol(expr), "samples \n")
message("metadata dim:", nrow(metadata), "x", ncol(metadata), "\n")

# Rows and columns (quick peek)
message("First gene IDs:\n", head(rownames(expr)), "\n")
message("First sample IDs:\n", head(colnames(expr)), "\n")

message("Metadata fields:\n", paste(colnames(metadata), collapse = ", "), "\n")

#-------------------------------------------------------------------------------
# 3. Basic checks
#-------------------------------------------------------------------------------
message("Are there NA values in expr?:", anyNA(expr), "\n")
message("Are there NA values in metadata?:", anyNA(metadata), "\n")

message("Expression value range:", paste(range(expr, na.rm = TRUE), collapse = " .. "), "\n")

message("expr summary:\n")
print(summary(as.vector(expr)))

message("Example of expression values:\n")
print(expr[1:5, 1:5])

message("Example of metadata:\n")
print(metadata[1:5, ])

#-------------------------------------------------------------------------------
# 4. Consistency between expr and metadata
#-------------------------------------------------------------------------------
message("Samples in expr:", ncol(expr), "\n",
    "Unique samples:", length(unique(colnames(expr))), "\n")

message("Samples in metadata:", nrow(metadata), "\n",
    "Unique samples:", length(unique(metadata$Sample.Name)), "\n")

message("Any duplicates in metadata$Sample.Name?:",
    any(duplicated(metadata$Sample.Name)), "\n")

message("Any duplicated genes in expr rownames?:",
    any(duplicated(rownames(expr))), "\n")

message("Same sample set (expr vs metadata):",
    setequal(colnames(expr), metadata$Sample.Name), "\n")

message("Same sample order:",
    identical(colnames(expr), metadata$Sample.Name), "\n")

# Reorder metadata if needed
if (!identical(colnames(expr), metadata$Sample.Name)) {
  metadata <- metadata[match(colnames(expr), metadata$Sample.Name), ]
  message("Metadata reordered to match expr sample order.\n")
}

stopifnot(identical(colnames(expr), metadata$Sample.Name))

#-------------------------------------------------------------------------------
# 5. Dataset composition
#-------------------------------------------------------------------------------
n_patients <- length(unique(metadata$Donor))
message("Number of donors:", n_patients, "\n")

# Tables
samples_per_donor <- table(metadata$Donor)
celltype_table <- table(metadata$CellType)
organ_table <- table(metadata$Organ)
donor_celltype <- table(metadata$Donor, metadata$CellType)

message("\nSamples per donor:\n")
print(samples_per_donor)

message("\nCell type distribution:\n")
print(celltype_table)

message("\nOrgan distribution:\n")
print(organ_table)

message("\nDonor x CellType table:\n")
print(donor_celltype)

#-------------------------------------------------------------------------------
# 6. Plots: composition
#-------------------------------------------------------------------------------
# Samples per donor
png(file.path(out_dir, "n_samples_donor.png"), width = 1600, height = 900, res = 200)
barplot(samples_per_donor,
        main = "Number of samples per donor",
        ylab = "Number of samples",
        las = 2)
dev.off()

# CellType distribution
png(file.path(out_dir, "cell_type_distribution.png"), width = 1600, height = 900, res = 200)
barplot(celltype_table,
        main = "Cell type distribution",
        ylab = "Number of samples",
        las = 2)
dev.off()

# Organ distribution
png(file.path(out_dir, "organ_distribution.png"), width = 1600, height = 900, res = 200)
barplot(organ_table,
        main = "Organ distribution",
        ylab = "Number of samples",
        las = 2)
dev.off()

#-------------------------------------------------------------------------------
# 7. Global expression density
#-------------------------------------------------------------------------------
png(file.path(out_dir, "global_expression_density.png"), width = 1600, height = 900, res = 200)
plot(density(as.vector(expr)),
     main = "Global expression density (normalized data)",
     xlab = "Expression value",
     ylab = "Density")
dev.off()

#-------------------------------------------------------------------------------
# 8. PCA
#-------------------------------------------------------------------------------
gene_var <- apply(expr, 1, var)
message("\nGenes with zero variance:", sum(gene_var == 0), "\n")

pca <- prcomp(t(expr), center = TRUE, scale. = TRUE)

pca_res <- as.data.frame(pca$x[, 1:2])
pca_res$Sample.Name <- rownames(pca_res)

# PCA coloured by CellType
pca_celltype <- merge(
  pca_res,
  metadata[, c("Sample.Name", "CellType")],
  by = "Sample.Name"
)

p_celltype <- ggplot(pca_celltype, aes(x = PC1, y = PC2, colour = CellType)) +
  geom_point() +
  labs(title = "PCA (colour = CellType)") +
  theme_minimal()

ggsave(
  filename = file.path(out_dir, "PCA_celltype.png"),
  plot = p_celltype,
  width = 6,
  height = 5,
  dpi = 300
)

# PCA coloured by Donor
pca_donor <- merge(
  pca_res,
  metadata[, c("Sample.Name", "Donor")],
  by = "Sample.Name"
)

p_donor <- ggplot(pca_donor, aes(x = PC1, y = PC2, colour = Donor)) +
  geom_point() +
  labs(title = "PCA (colour = Donor)") +
  theme_minimal()

ggsave(
  filename = file.path(out_dir, "PCA_donor.png"),
  plot = p_donor,
  width = 6,
  height = 5,
  dpi = 300
)

#-------------------------------------------------------------------------------
# 9. Save results
#-------------------------------------------------------------------------------
eda_res <- list(
  expr_dim = dim(expr),
  metadata_dim = dim(metadata),
  samples_per_donor = samples_per_donor,
  celltype_table = celltype_table,
  organ_table = organ_table,
  donor_celltype = donor_celltype,
  gene_var_zero = sum(gene_var == 0),
  pca = pca,
  pca_scores = pca_res,
  pca_celltype = pca_celltype,
  pca_donor = pca_donor
)

saveRDS(eda_res, file.path(out_dir, "eda_results.rds"))

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)

message("\nEDA finished. Outputs saved to: ", out_dir, "\n")
