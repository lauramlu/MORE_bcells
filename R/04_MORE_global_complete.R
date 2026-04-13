################################################################################
################################################################################
# Run complete global MORE model
# Author: Laura Muñoz
# Date: 19/03/2026
#
# Goal:
# Run the final global MORE model integrating RNA expression, TFBS, ATAC and
# donor effect using the PLS1 method.
#
# Input:
# - targetData.rds
# - regulatoryData.rds
# - associations.rds
# - conditions_CellType.rds
# - clinic.rds
# - clinicType.rds
#
# Output:
# - out_more_global_complete.rds
# - regpcond_complete.rds
# - regpcond_filtered_complete.rds
################################################################################

#-------------------------------------------------------------------------------
# 0. Set up
#-------------------------------------------------------------------------------
rm(list = ls()) 
gc() 

# print R version and working directory message("R version: ", R.version.string) 
message("Working directory: ", getwd()) 

# Packages 
library(dplyr) ; library(MORE) 

in_dir <- "02_results/04_MORE_complete/00_inputs" 
out_dir <- "02_results/04_MORE_complete" 

#------------------------------------------------------------------------------- 
# 1. Data 
#------------------------------------------------------------------------------- 
targetData <- readRDS(file.path(in_dir, "targetData.rds")) 
regulatoryData <- readRDS(file.path(in_dir, "regulatoryData.rds")) 
associations <- readRDS(file.path(in_dir, "associations.rds")) 
conditions <- readRDS(file.path(in_dir, "conditions_CellType.rds"))
clinic <- readRDS(file.path(in_dir, "clinic.rds"))
clinicType <- readRDS(file.path(in_dir, "clinicType.rds"))

#-------------------------------------------------------------------------------
# 2. MORE
#-------------------------------------------------------------------------------
out <- more(
  targetData     = targetData,
  regulatoryData = regulatoryData,
  associations   = associations,
  condition      = conditions,
  clinic         = clinic, 
  clinicType     = clinicType, 
  method         = "PLS1",
  varSel         = "Jack",
  alfa           = 0.05,
  vip            = 0.8,
  minVariation   = 0,
  scaleType      = "auto",
  parallel       = FALSE
  )

regpcond <- RegulationPerCondition(out)
regpcond_filtered <- FilterRegulationPerCondition(out, regpcond, filterR2 = 0.5)

#-------------------------------------------------------------------------------
# 3. save outputs
#-------------------------------------------------------------------------------    
saveRDS(out, file = file.path(out_dir,"out_more_global_complete.rds"))
saveRDS(regpcond, file = file.path(out_dir,"regpcond_complete.rds"))
saveRDS(regpcond_filtered, file = file.path(out_dir, "regpcond_filtered_complete.rds"))

