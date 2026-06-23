################################################################################
# Differential expression analysis: limma-voom
# Author: Laura Muñoz
# 11/02/2026
#
# Goal:
# Identify differentially expressed genes across B-cell developmental transitions
# and generate the DEG union used later as target features for MORE.
#
# Input:
# - 00_data/raw/rna_raw_data.rds
# - 00_data/raw/rna_metadata.rds
#
# Output:
# - deg_by_transition_adjP0.05_absFC1.5.rds
# - deg_union_adjP0.05_absFC1.5.rds
# - voom_v.rds
# - model objects and summary files
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
library(ggplot2) ; library(limma) ; library(edgeR) ; library(dplyr)

in_dir <- "00_data/raw"
out_dir <- "02_results/02_DEG"

#-------------------------------------------------------------------------------
# 1. Data
#-------------------------------------------------------------------------------
rna_counts <- readRDS(file.path(in_dir, "rna_raw_data.rds"))   # gene x sample raw counts
metadata   <- readRDS(file.path(in_dir, "rna_metadata.rds"))  # sample-level metadata

#-------------------------------------------------------------------------------
# 2. Object type, structure and dimensions 
#-------------------------------------------------------------------------------
message("rna_counts class:", class(rna_counts), "\n")
message("metadata class:", class(metadata), "\n")

message("rna_counts dim:", nrow(rna_counts), "genes x", ncol(rna_counts), "samples \n")
message("metadata dim:", nrow(metadata), "x", ncol(metadata), "\n")

message("First gene IDs:\n", head(rownames(rna_counts)), "\n")
message("First sample IDs:\n", head(colnames(rna_counts)), "\n")

# Available metadata fields
message("Information for each sample in metadata (col names):\n", 
    paste(colnames(metadata), collapse = ", "), "\n")

## -----------------------------------------------------------------------------
# 3. Basic checks
## -----------------------------------------------------------------------------
message("Are there NA values in rna_counts?:", anyNA(rna_counts), "\n")
message("Are there NA values in metadata?:", anyNA(metadata), "\n")

# Partial view of the counts matrix 
message("Example of expression values:\n")
print(rna_counts[1:5, 1:5])

# Partial view of metadata data frame
message("Example of metadata:\n")
print(metadata[1:5,])

## -----------------------------------------------------------------------------
# 4. Consistency between rna_counts and metadata
## -----------------------------------------------------------------------------
message("Samples in rna_counts:", ncol(rna_counts), "\n", 
    "Unique samples:", length(unique(colnames(rna_counts))), "\n")

message("Samples in metadata:", nrow(metadata), "\n", 
    "Unique samples:", length(unique(metadata$Sample.Name)), "\n")

message("Are there any duplicates in Sample.Name entries in metadata? \n",
    any(duplicated(metadata$Sample.Name)), "\n")

message("Are there any duplicates genes in rna_counts? \n",
    any(duplicated(row.names(rna_counts))), "\n")

## -----------------------------------------------------------------------------
# 5. Sample identity and order check
## -----------------------------------------------------------------------------
message("Same sample set (rna_counts vs metadata):",
    all(metadata$Alias %in% colnames(rna_counts)), "\n")

message("Same sample order:",
    identical(colnames(rna_counts), metadata$Alias), "\n")

# Reorder metadata to match count columns if needed
if (!identical(colnames(rna_counts), metadata$Alias)) {
  metadata <- metadata[match(colnames(rna_counts), metadata$Alias), ]
  message("Metadata reordered to match rna_counts sample order.\n")
}

stopifnot(identical(colnames(rna_counts), metadata$Alias))

## -----------------------------------------------------------------------------
# 6. Design matrix (cell types and donor block)  
## -----------------------------------------------------------------------------
# CellType and Donor must be factors for model.matrix()
metadata$CellType <- factor(metadata$CellType)
metadata$Donor <- factor(metadata$Donor)

group <- metadata$CellType
donor <- metadata$Donor

# Design without intercept:
# - One coefficient per cell state (estimated mean log2-CPM per state)
# - Donor terms adjust for repeated measures / donor-specific effects
# Note: one donor level becomes the reference (hence 10 donor columns for 11 donors)
design <- model.matrix(~ 0 + group + donor)

# Clean column names: remove group prefix
colnames(design) <- gsub("group","",colnames(design))

message("Design matrix dim: ", dim(design), "\n")
message("First 20 design columns:\n")
print(head(colnames(design), 20))

## -----------------------------------------------------------------------------
# 7. Create DGEList + filter (low expression genes)  
## -----------------------------------------------------------------------------
dge <- DGEList(counts = rna_counts, group = group)

# filterByExpr keeps genes with sufficient expression given the design/groups
keep <- filterByExpr(dge, design = design, group = group, min.count = 1)

message("Genes kept after filterByExpr:\n")
print(table(keep))

dge_f <- dge[keep, , keep.lib.sizes = FALSE]
message("Filtered DGEList dim: ", dim(dge_f), "\n")

## -----------------------------------------------------------------------------
# 8. TMM normalization  
## -----------------------------------------------------------------------------
dge_f <- calcNormFactors(dge_f, method = "TMM")

message("TMM norm factors summary:\n")
print(summary(dge_f$samples$norm.factors))

message("Library sizes and norm factors (first samples):\n")
print(head(dge_f$samples[, c("lib.size","norm.factors")]))

## -----------------------------------------------------------------------------
# 9. Voom transformation (mean-dependent variance modeling)
## -----------------------------------------------------------------------------
# voom transforms counts to log2-CPM and estimates precision weights
# to account for mean-dependent variance in RNA-seq.
v <- voom(dge_f, design = design)

message("Voom E dim: ", dim(v$E), "\n")

## -----------------------------------------------------------------------------
# 10. Fit linear model (limma) + define transition contrasts
## -----------------------------------------------------------------------------
fit <- lmFit(v, design)

message("Coefficients dim (genes x design cols): ", dim(fit$coefficients), "\n")

# Transition contrasts (state-to-state) 
contr <- makeContrasts(
  CLP_vs_HSC = CLP - HSC,
  proB_vs_CLP = proB - CLP,
  preB_vs_proB = preB - proB,
  ImmatureB_vs_preB = ImmatureB - preB,
  Transitional_B_vs_ImmatureB = Transitional_B - ImmatureB,
  Naive_CD5pos_vs_Transitional_B = Naive_CD5pos - Transitional_B,
  Naive_CD5neg_vs_Transitional_B = Naive_CD5neg - Transitional_B,
  levels = design
)

message("Contrast matrix dim:", dim(contr), "\n")
message("Contrasts:\n")
print(colnames(contr))

## -----------------------------------------------------------------------------
# 11. Apply contrasts + eBayes
## -----------------------------------------------------------------------------
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)

message("Fit2 coefficients dim (genes x contrasts): ", dim(fit2$coefficients), "\n")

## -----------------------------------------------------------------------------
# 12. DEGs per transition
## -----------------------------------------------------------------------------
# DEGs per transition following Planell et al.:
# adjusted P value (BH) <= 0.05 and absolute fold change >= 1.5
# limma reports log2 fold change, so FC >= 1.5 corresponds to:
# abs(logFC) >= log2(1.5)

padj_thr <- 0.05
fc_thr <- 1.5
logfc_thr <- log2(fc_thr)

deg_by_transition <- lapply(colnames(fit2$coefficients), function(coef_name) {
  tt <- topTable(
    fit2,
    coef = coef_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  
  tt[
    tt$adj.P.Val <= padj_thr &
      abs(tt$logFC) >= logfc_thr,
  ]
})

names(deg_by_transition) <- colnames(fit2$coefficients)

# Counts per transition
deg_counts <- sapply(deg_by_transition, nrow)

message("DEGs per transition: \n")
print(deg_counts)

# DEG union = dynamic genes, significant in at least one transition
# following Planell et al. criteria: adj.P.Val <= 0.05 and |FC| >= 1.5
deg_union <- unique(unlist(lapply(deg_by_transition, rownames)))

message(
  "Union of DEGs across transitions (unique genes), adj.P.Val <= ",
  padj_thr,
  " and |FC| >= ",
  fc_thr,
  ": ",
  length(deg_union),
  "\n"
)

## -----------------------------------------------------------------------------
# 13. Save objects and reproducibility info
## -----------------------------------------------------------------------------
saveRDS(design, file.path(out_dir, "design.rds"))
saveRDS(dge_f,  file.path(out_dir, "dge_f.rds"))
saveRDS(v,      file.path(out_dir, "voom_v.rds"))

saveRDS(fit,   file.path(out_dir, "fit_lmFit.rds"))
saveRDS(contr, file.path(out_dir, "contrasts_matrix.rds"))
saveRDS(fit2,  file.path(out_dir, "fit2_ebayes.rds"))

# DEG results
saveRDS(
  deg_by_transition,
  file.path(
    out_dir,
    sprintf("deg_by_transition_adjP%.2f_absFC%.1f.rds", padj_thr, fc_thr)
  )
)

saveRDS(
  deg_union,
  file.path(
    out_dir,
    sprintf("deg_union_adjP%.2f_absFC%.1f.rds", padj_thr, fc_thr)
  )
)

write.csv(
  deg_counts,
  file.path(
    out_dir,
    sprintf("deg_counts_by_transition_adjP%.2f_absFC%.1f.csv", padj_thr, fc_thr)
  )
)

writeLines(
  deg_union,
  file.path(
    out_dir,
    sprintf("deg_union_for_MORE_adjP%.2f_absFC%.1f.txt", padj_thr, fc_thr)
  )
)

# Summaries
write.csv(
  deg_counts,
  file.path(out_dir, sprintf("deg_counts_by_transition_adjP%.2f.csv", padj_thr))
)

writeLines(
  deg_union,
  file.path(out_dir, sprintf("deg_union_for_MORE_adjP%.2f.txt", padj_thr))
)

# Session info for reproducibility
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)


## -----------------------------------------------------------------------------
# 14. PCA: analysis of filtered data
## -----------------------------------------------------------------------------
v_expr_deg <- v$E[deg_union,]

pca <- prcomp(t(v_expr_deg), scale. = TRUE)

pca_res <- as.data.frame(pca$x[, 1:2])
pca_res$Alias <- rownames(pca_res)

# PCA coloured by CellType
pca_celltype <- merge(
  pca_res,
  metadata[, c("Sample.Name", "CellType", "Alias", "Donor")],
  by = "Alias"
)

p_celltype <- ggplot(pca_celltype, aes(x = PC1, y = PC2, colour = CellType)) +
  geom_point() +
  labs(title = "PCA - DEG filtered data") +
  theme_minimal()

p_donor <- ggplot(pca_celltype, aes(x= PC1, y = PC2, colour = Donor)) +
  geom_point() + 
  labs(title = "PCA (colour = Donor)") + 
  theme_minimal()

ggsave(
  filename = file.path(out_dir, "PCA_celltype_DEG.png"),
  plot = p_celltype,
  width = 6,
  height = 5,
  dpi = 300
)
