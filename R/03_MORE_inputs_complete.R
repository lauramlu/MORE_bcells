################################################################################
# Prepare complete multi-omic inputs for MORE
# Author: Laura Muñoz
# Date: 18/03/2026
#
# Goal:
# Prepare the complete MORE input objects by integrating RNA expression, TFBS,
# ATAC and donor information.
#
# Input:
# - 02_results/02_DEG/voom_v.rds
# - 02_results/02_DEG/deg_union_adjP0.05_absFC1.5.rds
# - 00_data/raw/rna_metadata.rds
# - 00_data/processed/atac_norm_data.rds
# - 00_data/processed/atac_dba_by_transitions.rds
# - 00_data/processed/TFBS_desviations_curated.rds
# - 00_data/processed/atac_OCR_annotation.rds
# - 00_data/processed/tf_ocr_matrix.rds
#
# Output:
# - targetData.rds
# - regulatoryData.rds
# - associations.rds
# - conditions_CellType.rds
# - clinic.rds
# - clinicType.rds
################################################################################
#-------------------------------------------------------------------------------
# 1. Set up
#-------------------------------------------------------------------------------
rm(list = ls())
gc()

# print R version and working directory
message("R version: ", R.version.string)
message("Working directory: ", getwd())

# Packages
library(dplyr) ; library(tidyr) ; library(tibble)

#-------------------------------------------------------------------------------
# 1. Paths
#-------------------------------------------------------------------------------
data_raw_dir       <- "00_data/raw"
data_processed_dir <- "00_data/processed"
results_dir        <- "02_results/02_DEG"
out_dir            <- "02_results/04_MORE_complete/00_inputs/"

#-------------------------------------------------------------------------------
# 2. Data
#-------------------------------------------------------------------------------
# RNA (voom object)
v <- readRDS(file.path(results_dir, "voom_v.rds"))
expr <- v$E  # genes x samples

# DEG union
deg <- readRDS(file.path(results_dir, "deg_union_adjP0.05_absFC1.5.rds"))

# metadata
metadata <- readRDS(file.path(data_raw_dir, "rna_metadata.rds"))

# ATAC
atac <- readRDS(file.path(data_processed_dir, "atac_norm_data.rds"))
dba_trans <- readRDS(file.path(data_processed_dir, "atac_dba_by_transitions.rds"))

# TFBS
tfbs <- readRDS(file.path(data_processed_dir, "TFBS_desviations_curated.rds"))

# prior knowledge
ocr_gene_annotation <- readRDS(file.path(data_processed_dir, "atac_OCR_annotation.rds"))

tf_ocr_annotation <- readRDS(file.path(data_processed_dir, "tf_ocr_matrix.rds"))

#-------------------------------------------------------------------------------
# 3. Consistency checks + align metadata to expr
#-------------------------------------------------------------------------------
stopifnot("Alias" %in% colnames(metadata))
stopifnot(all(c("Sample.Name", "CellType") %in% colnames(metadata)))

# Align metadata rows to expr columns
stopifnot(setequal(colnames(expr), metadata$Alias))
metadata <- metadata[match(colnames(expr), metadata$Alias), , drop = FALSE]

stopifnot(identical(colnames(expr), metadata$Alias))
stopifnot(!anyDuplicated(metadata$Sample.Name))

#-------------------------------------------------------------------------------
# 4. filtering genes (DEG)
#-------------------------------------------------------------------------------
# adj.P.Val <= 0.05 and |FC| >= 1.5

# Filter genes + rename samples
expr_deg <- expr[deg, , drop = FALSE]
colnames(expr_deg) <- metadata$Sample.Name

stopifnot(identical(colnames(expr_deg), metadata$Sample.Name))

# DENSITY PLOT
plot(density(as.vector(expr_deg)),
     main = "Global expression density (normalized and filtered DEG data)",
     xlab = "Expression value",
     ylab = "Density")

#-------------------------------------------------------------------------------
# 5. regulatoryData: TFBS and ATAC. sync samples with expr_deg
#-------------------------------------------------------------------------------
stopifnot(!anyDuplicated(colnames(expr_deg)))
stopifnot(!anyDuplicated(colnames(tfbs)))
stopifnot(!anyDuplicated(colnames(atac)))

# Common samples 
common_samples <- Reduce(
  intersect,
  list(colnames(expr_deg), colnames(tfbs), colnames(atac))
)
message("Common samples: ", length(common_samples))

stopifnot(length(common_samples) > 0)

expr_matrix <- expr_deg[, common_samples, drop = FALSE]
tf_matrix <- tfbs[ , common_samples, drop = FALSE]
atac_matrix <- atac[ , common_samples, drop = FALSE]

# Force same order
tf_matrix <- tf_matrix[ , colnames(expr_matrix), drop = FALSE]
atac_matrix <- atac_matrix[ , colnames(expr_matrix), drop = FALSE]

stopifnot(identical(colnames(expr_matrix), colnames(tf_matrix)))
stopifnot(identical(colnames(expr_matrix), colnames(atac_matrix)))

# Create regulatoryData and targetData
regulatoryData <- list(
  "TFBS" = tf_matrix,
  "ATAC" = atac_matrix)
targetData <- expr_matrix

#-------------------------------------------------------------------------------
# 6. conditions: experimental design for MORE
#-------------------------------------------------------------------------------
# Must have the same samples in the same order as targetData and regulatoryData
stopifnot(identical(rownames(metadata), metadata$Sample.Name))
conditions <- metadata[colnames(targetData), "CellType", drop = FALSE]

stopifnot(identical(rownames(conditions), colnames(targetData)))
stopifnot(identical(rownames(conditions), colnames(regulatoryData$TFBS)))
stopifnot(identical(rownames(conditions), colnames(regulatoryData$ATAC)))
conditions$CellType <- factor(conditions$CellType)
stopifnot(is.factor(conditions$CellType))

#-------------------------------------------------------------------------------
# 7. associations - TFBS
#-------------------------------------------------------------------------------
# Filter ocr_gene annotation
ocr_gene_matrix <- ocr_gene_annotation %>%
  rownames_to_column("OCR") %>%
  select(ENSEMBL,OCR) %>%
  filter(!is.na(ENSEMBL))

rownames(ocr_gene_matrix) <- NULL

# Filter tf_ocr matrix
tf_ocr_long <- tf_ocr_annotation %>%
  as.data.frame() %>%
  rownames_to_column("OCR") %>%
  pivot_longer(
    cols = -OCR,
    names_to = "TF",
    values_to = "value"
  ) %>%
  filter(value == 1) %>%
  select(OCR, TF)

tf_ocr_gene <- inner_join(tf_ocr_long, ocr_gene_matrix, by = "OCR")
tf_gene <- tf_ocr_gene %>%
  distinct(TF, ENSEMBL)

length(intersect(tf_ocr_long$OCR, ocr_gene_matrix$OCR))

# filter genes (DEG)
tf_gene_deg <- tf_gene %>%
  filter(ENSEMBL %in% deg) %>%
  select(ENSEMBL, TF)

# basic checks
sum(!tf_gene_deg$TF %in% rownames(regulatoryData$TFBS))
setdiff(unique(tf_gene_deg$TF), rownames(regulatoryData$TFBS))
grep("^NKX", rownames(regulatoryData$TFBS), value = TRUE)

tf_gene_deg_fixed <- as.data.frame(tf_gene_deg) %>%
  mutate(TF = case_when(
    TF == "NKX2.3" ~ "NKX2-3",
    TF == "NKX3.1" ~ "NKX3-1",
    TRUE ~ TF
  ))

stopifnot(all(tf_gene_deg_fixed$TF %in% rownames(regulatoryData$TFBS)))

stopifnot(all(tf_gene_deg_fixed$ENSEMBL %in% rownames(targetData)))

#-------------------------------------------------------------------------------
# 8. associations - ATAC
#-------------------------------------------------------------------------------
# Filter OCR to include only regulators, target features and area
associations_atac <- ocr_gene_annotation %>% 
  mutate(REGULATOR = rownames(ocr_gene_annotation)) %>%
  filter(!is.na(ENSEMBL), ENSEMBL %in% rownames(targetData)) %>%
  select(ENSEMBL, REGULATOR, AREA = annotation_simplified) 

rownames(associations_atac) <- NULL

# DBA: keep OCRs active in >= 1 transition
sig_cols <- grep("^sig_", colnames(dba_trans)) 

# Eliminate column "sig_HSC_Transitional.B" as in original paper
sig_cols <- sig_cols[colnames(dba_trans)[sig_cols] != "sig_HSC_Transitional.B"]

sig_dba <- dba_trans[ , sig_cols]
is_DA <- apply(sig_dba, 1, function(x) any(x != 0))
DA_OCRs <- rownames(dba_trans)[is_DA]

stopifnot(all(DA_OCRs %in% rownames(regulatoryData$ATAC)))

# keep OCRs Promoter or Distal Intergenic and significative
keep <- associations_atac$AREA == "Promoter" |
  (associations_atac$AREA == "Distal Intergenic" &
     associations_atac$REGULATOR %in% DA_OCRs)

associations_filtered <- associations_atac[keep, , drop = FALSE]

# There are OCRs in annotation not present in atac_norm_data
associations_filtered <- associations_filtered[
  associations_filtered$REGULATOR %in% rownames(regulatoryData$ATAC),
  , drop = FALSE
]

# sanity checks
stopifnot(nrow(associations_filtered) > 0)
stopifnot(all(associations_filtered$REGULATOR %in% rownames(regulatoryData$ATAC)))
stopifnot(all(associations_filtered$ENSEMBL %in% rownames(targetData)))

message("Final associations:\n")
message("Unique OCRs:", length(unique(associations_filtered$REGULATOR)), "\n")
message("Unique genes:", length(unique(associations_filtered$ENSEMBL)), "\n")

associations <- list(
  "TFBS" = tf_gene_deg_fixed,
  "ATAC" =associations_filtered
  )

stopifnot(nrow(associations_filtered) > 0)

ocrs_gen <- table(associations_filtered$ENSEMBL)
summary(as.vector(ocrs_gen))

#-------------------------------------------------------------------------------
# 9. clinic
#-------------------------------------------------------------------------------
clinic <- metadata[colnames(targetData), "Donor", drop = FALSE]

stopifnot(identical(rownames(clinic), colnames(targetData)))
stopifnot(identical(rownames(clinic), colnames(regulatoryData$TFBS)))
stopifnot(identical(rownames(clinic), colnames(regulatoryData$ATAC)))

clinic$Donor <- factor(clinic$Donor)
stopifnot(is.factor(clinic$Donor))

clinicType <- c(Donor = 1)

#-------------------------------------------------------------------------------
# 10. Save inputs
#-------------------------------------------------------------------------------
saveRDS(targetData, file.path(out_dir, "targetData.rds"))
saveRDS(conditions, file.path(out_dir, "conditions_CellType.rds"))
saveRDS(regulatoryData, file.path(out_dir, "regulatoryData.rds"))
saveRDS(associations, file.path(out_dir, "associations.rds"))
saveRDS(clinic, file.path(out_dir, "clinic.rds"))
saveRDS(clinicType, file.path(out_dir, "clinicType.rds"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("\nSaved MORE inputs to: ", out_dir)

sessionInfo()


