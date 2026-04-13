# Regulatory analysis of early human B-cell differentiation using MORE

## Overview

This project performs a multi-omics regulatory analysis of early human B-cell differentiation using the MORE (Multi-Omics REgulation) framework. The aim is to identify transcriptional regulatory programs across developmental stages and to characterise key transcription factors driving cell state transitions.

This work was developed as part of an extracurricular internship in bioinformatics.

---

## Objectives

- Explore RNA-seq data across B-cell developmental stages  
- Identify differentially expressed genes (DEGs)  
- Integrate multi-omics regulatory layers (TFBS and ATAC-seq)  
- Build and evaluate a global MORE model  
- Analyse regulatory activity across cell types  
- Interpret biological results  

---

## Data

The analysis integrates multiple omics layers:

- **RNA-seq**: gene expression data (voom-normalised log2-CPM)  
- **ATAC-seq**: chromatin accessibility data  
- **TFBS deviations**: transcription factor activity estimates  

> Raw and intermediate data are not included in this repository due to size and data management constraints.  
> This repository focuses on the analysis workflow and selected outputs.

---

## Repository structure

```
.
├── MORE.Rproj
├── R
│   ├── 01_exploratory_analysis.R
│   ├── 02_DEG.R
│   ├── 03_MORE_inputs_complete.R
│   ├── 04_MORE_global_complete.R
│   └── 05_MORE_results_complete.R
├── docs
├── renv
│   ├── activate.R
│   └── settings.json
└── renv.lock
```

---

## Workflow

The analysis follows these main steps:

1. **Exploratory analysis**  
   Initial inspection and quality assessment of the dataset.

2. **Differential expression analysis (DEG)**  
   Identification of dynamically regulated genes across conditions.

3. **Input preparation for MORE**  
   Construction of:
   - targetData (RNA expression)
   - regulatoryData (TFBS + ATAC)
   - associations (TF–gene and OCR–gene relationships)
   - conditions (cell type labels)

4. **Global MORE model**  
   Integration of multi-omics data using PLS1-based modelling.

5. **Results analysis**  
   - Regulatory activity per condition  
   - Identification of key transcription factors  
   - Functional interpretation  

---

## Reproducibility

This project uses `renv` to ensure reproducibility.
A reference session snapshot is available in docs/session_Info.txt.

To reproduce the environment:

```r
renv::restore()
```

---

## Requirements

- R (≥ 4.x)

- Key packages:
  - MORE  
  - dplyr  
  - AnnotationDbi  
  - org.Hs.eg.db  
  - GO.db  

All package versions are recorded in `renv.lock`.

---

## Main reference

Planell et al. (2025)  
*Uncovering the regulatory landscape of early human B cell lymphopoiesis and its implications in the pathogenesis of B-ALL*  
Science Advances  

---

## Author

Laura Muñoz  
Bioinformatics MSc student  

---

## Notes

This repository is intended to showcase:
- multi-omics data integration  
- regulatory modelling using MORE  
- structured and reproducible analysis workflows  

It accompanies the internship report included in the `docs/` folder.
