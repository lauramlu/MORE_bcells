################################################################################
# Analyse results of the complete global MORE model
# Author: Laura Muñoz
# Date: 01/04/2026
#
# Goal:
# Analyse the results of the complete global MORE model, including state-specific
# regulation, differential regulator analysis and functional interpretation.
#
# Input:
# - out_more_global_complete.rds
# - regpcond_complete.rds
# - regpcond_filtered_complete.rds
#
# Output:
# - regulation_by_celltype_complete.rds
# - regulator summary tables
# - functional interpretation results
# - plots
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
library(dplyr) ; library(MORE) ; library(AnnotationDbi) ; library(org.Hs.eg.db)
library(GO.db)

in_dir  <- "02_results/04_MORE_complete"
out_dir <- "02_results/04_MORE_complete/01_plots"

#-------------------------------------------------------------------------------
# 1. Data
#-------------------------------------------------------------------------------
out_more_global <- readRDS(file.path(in_dir, "out_more_global_complete.rds"))
regpcond <- readRDS(file.path(in_dir, "regpcond_complete.rds"))
regpcond_filtered <- readRDS(file.path(in_dir, "regpcond_filtered_complete.rds"))

#-------------------------------------------------------------------------------
# 2. Results
#-------------------------------------------------------------------------------
# summary
summaryPlot(out_more_global, regpcond, byTargetF = FALSE, filterR2 = 0.5)

print(summary(out_more_global))

celltypes <- c(
  "HSC", 
  "CLP", 
  "proB", 
  "preB", 
  "ImmatureB", 
  "Transitional_B", 
  "Naive_CD5pos", 
  "Naive_CD5neg"
  )

regulation_by_celltype <- vector("list", length(celltypes))
names(regulation_by_celltype) <- celltypes

for (celltype in celltypes) {
  message("Processing: ", celltype)
  
  regulation_by_celltype[[celltype]] <- 
      RegulationInCondition(regpcond_filtered, celltype)
}

saveRDS(regulation_by_celltype, file.path(
  in_dir, 
  "regulation_by_celltype_complete.rds"
  ))

#-------------------------------------------------------------------------------
# 3. preB vs proB analysis: global regulators
#-------------------------------------------------------------------------------
reg_preB <- regulation_by_celltype$preB
reg_proB <- regulation_by_celltype$proB

# regulators exclusive of each condition
preB_only <- setdiff(reg_preB$GlobalRegulators, reg_proB$GlobalRegulators)
proB_only <- setdiff(reg_proB$GlobalRegulators, reg_preB$GlobalRegulators)

#-------------------------------------------------------------------------------
# 4. differential regulators preB vs proB (activation vs repression)
#-------------------------------------------------------------------------------
differential_tf_preB_proB <- regpcond_filtered %>%
  filter(omic == "TFBS") %>%
  group_by(regulator) %>%
  summarise(
    n_targetF = n_distinct(targetF),
    
    # proB
    proB_activation = sum(Group_proB[Group_proB > 0]),
    proB_repression = sum(Group_proB[Group_proB < 0]),
    
    # preB
    preB_activation = sum(Group_preB[Group_preB > 0]),
    preB_repression = sum(Group_preB[Group_preB < 0])
    
  ) %>%
  mutate(
    delta_activation = preB_activation - proB_activation,
    delta_repression = preB_repression - proB_repression
  ) %>%
  filter(n_targetF > 10) # regulators with very few targets removed

# delta > 0 -> stronger regulation in preB
# delta < 0 -> stronger regulation in proB

top_preB_activation <- differential_tf_preB_proB %>%
  arrange(desc(delta_activation)) %>%
  dplyr::select(regulator, n_targetF, preB_activation, delta_activation)

top_proB_activation <- differential_tf_preB_proB %>%
  arrange(delta_activation) %>%
  dplyr::select(regulator, n_targetF, proB_activation, delta_activation)

top_preB_repression <- differential_tf_preB_proB %>%
  arrange(delta_repression) %>%
  dplyr::select(regulator, n_targetF, preB_repression, delta_repression)
  
top_proB_repression <- differential_tf_preB_proB %>%
  arrange(desc(delta_repression)) %>% 
  dplyr::select(regulator, n_targetF, proB_repression, delta_repression)

head(top_preB_activation,20)
head(top_preB_repression,20)
head(top_proB_activation, 20)
head(top_proB_repression, 20)

#-------------------------------------------------------------------------------
# 5. CLP vs HSC analysis: global regulators
#-------------------------------------------------------------------------------
reg_CLP <- regulation_by_celltype$CLP
reg_HSC <- regulation_by_celltype$HSC

# regulators exclusive of each condition
CLP_only <- setdiff(reg_CLP$GlobalRegulators, reg_HSC$GlobalRegulators)
HSC_only <- setdiff(reg_HSC$GlobalRegulators, reg_CLP$GlobalRegulators)

#-------------------------------------------------------------------------------
# 6. differential regulators CLP vs HSC
#-------------------------------------------------------------------------------
differential_tf_CLP_HSC <- regpcond_filtered %>%
  filter(omic == "TFBS") %>%
  group_by(regulator) %>%
  summarise(
    n_targetF = n_distinct(targetF),
    
    # HSC
    HSC_activation = sum(Group_HSC[Group_HSC > 0]),
    HSC_repression = sum(Group_HSC[Group_HSC < 0]),
    
    # CLP
    CLP_activation = sum(Group_CLP[Group_CLP > 0]),
    CLP_repression = sum(Group_CLP[Group_CLP < 0])
    
  ) %>%
  mutate(
    delta_activation = CLP_activation - HSC_activation,
    delta_repression = CLP_repression - HSC_repression
  ) %>%
  filter(n_targetF > 10) # regulators with very few targets removed

# delta > 0 -> stronger regulation in CLP
# delta < 0 -> stronger regulation in HSC

top_CLP_activation <- differential_tf_CLP_HSC %>%
  arrange(desc(delta_activation)) %>%
  dplyr::select(regulator, n_targetF, CLP_activation, delta_activation)

top_HSC_activation <- differential_tf_CLP_HSC %>%
  arrange(delta_activation) %>%
  dplyr::select(regulator, n_targetF, HSC_activation, delta_activation)

top_CLP_repression <- differential_tf_CLP_HSC %>%
  arrange(delta_repression) %>%
  dplyr::select(regulator, n_targetF, CLP_repression, delta_repression)

top_HSC_repression <- differential_tf_CLP_HSC %>%
  arrange(desc(delta_repression)) %>%
  dplyr::select(regulator, n_targetF, HSC_repression, delta_repression)

head(top_CLP_activation,20)
head(top_CLP_repression,20)
head(top_HSC_activation, 20)
head(top_HSC_repression, 20)

#-------------------------------------------------------------------------------
# 7. functional analysis (ORA)
#-------------------------------------------------------------------------------
# preB
keys <- union(
  unique(reg_preB$RegulationInCondition$targetF), 
  unique(reg_proB$RegulationInCondition$targetF)
  )

# 1. gene -> GO
ora_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys,
  columns = c("GO", "ONTOLOGY"),
  keytype = "ENSEMBL"
) %>%
  filter(!is.na(GO), ONTOLOGY == "BP") %>%
  distinct(ENSEMBL, GO)

# 2. GO -> term description
go_terms_ORA <- AnnotationDbi::select(
  GO.db,
  keys = unique(ora_annotation$GO),
  columns = "TERM",
  keytype = "GOID"
) %>%
  distinct(GOID, TERM)

# 3. final annotation matrix for MORE
ora_annotation_more <- ora_annotation %>%
  left_join(go_terms_ORA, by = c("GO" = "GOID")) %>%
  transmute(
    targetF = ENSEMBL,
    GO_term = GO,
    GO_description = TERM
  ) %>%
  dplyr::filter(
    !is.na(targetF),
    !is.na(GO_term),
    !is.na(GO_description)
  ) %>%
  dplyr::distinct()

# 4. ORA

ora_more_preB <- oraMORE(
  out_more_global, 
  reg_preB, 
  byHubs = FALSE, 
  annotation = ora_annotation_more
)

# filter
# ora_more_preB <- ora_more_preB %>%
#   dplyr::filter(adjPval < 0.05)

ora_more_proB <- oraMORE(
  out_more_global, 
  reg_proB, 
  byHubs = FALSE, 
  annotation = ora_annotation_more
)

# filter
# ora_more_proB <- ora_more_proB %>%
#   dplyr::filter(adjPval < 0.05)

# analysis
# intersect(ora_more_preB$termDescr, ora_more_proB$termDescr)
# 
# setdiff(ora_more_preB$termDescr, ora_more_proB$termDescr)
# 
# setdiff(ora_more_proB$termDescr, ora_more_preB$termDescr)

#-------------------------------------------------------------------------------
# 8. GSEA between preB and proB
#-------------------------------------------------------------------------------
gsea_more_preB_vs_proB <- gseaMORE(
  reg_preB,
  reg_proB,
  annotation = ora_annotation_more
)

# If needed, inspect significant terms
gsea_more_preB_vs_proB_sig <- gsea_more_preB_vs_proB %>%
  dplyr::filter(p.adjust < 0.05)

#-------------------------------------------------------------------------------
# 9. Plot differential TFs: preB vs proB
#-------------------------------------------------------------------------------
plot_activation_proB <- top_proB_activation %>%
  slice_head(n = 10) %>%
  transmute(
    regulator,
    value = proB_activation,
    state = "proB",
    regulation = "Activation"
  )

plot_activation_preB <- top_preB_activation %>%
  slice_head(n = 10) %>%
  transmute(
    regulator,
    value = preB_activation,
    state = "preB",
    regulation = "Activation"
  )

plot_repression_proB <- top_proB_repression %>%
  slice_head(n = 10) %>%
  transmute(
    regulator,
    value = abs(proB_repression),
    state = "proB",
    regulation = "Repression"
  )

plot_repression_preB <- top_preB_repression %>%
  slice_head(n = 10) %>%
  transmute(
    regulator,
    value = abs(preB_repression),
    state = "preB",
    regulation = "Repression"
  )

plot_tf_diff <- bind_rows(
  plot_activation_proB,
  plot_activation_preB,
  plot_repression_proB,
  plot_repression_preB
) %>%
  mutate(
    panel = paste(state, regulation, sep = " - "),
    panel = factor(
      panel,
      levels = c(
        "proB - Activation",
        "preB - Activation",
        "proB - Repression",
        "preB - Repression"
      )
    ),
    regulator_panel = paste(panel, regulator, sep = "___")
  )

# Reordenación dentro de cada panel sin forcats
plot_tf_diff$regulator_panel <- with(
  plot_tf_diff,
  reorder(regulator_panel, value)
)

p_tf_diff <- ggplot(
  plot_tf_diff,
  aes(x = regulator_panel, y = value, fill = state)
) +
  geom_col(width = 0.7, show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~ panel, scales = "free_y", ncol = 2) +
  scale_x_discrete(
    labels = function(x) sub("^.*___", "", x)
  ) +
  labs(
    title = "Top differential TF regulators in the preB-proB transition",
    subtitle = "Top 10 regulators by activation or repression score in each state",
    x = NULL,
    y = "Regulatory activity score"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

p_tf_diff

ggsave(
  filename = file.path(out_dir, "top_tf_differential_preB_proB.png"),
  plot = p_tf_diff,
  width = 12,
  height = 8,
  dpi = 300
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))