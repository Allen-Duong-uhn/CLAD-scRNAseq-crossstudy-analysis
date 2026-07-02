# CLAD-scRNAseq-crossstudy-analysis
This repository contains the R-based analysis used to process, integrate and compare single-cell RNA sequencing (scRNAseq) datasets of lung and BAL macrophages in chronic lung allograft dysfunction (CLAD). The analysis focuses on harmonized cross-study comparison of the ISG/Super.Macro macrophages.

The R scripts used analyzed scRNAseq data from:
  1. Moshkelgosha et al. _J Heart Lung Transplant_. 2022. doi: 10.1016/j.healun.2022.05.005 
  2. Khatri et al. _JCI Insight_. 2023. doi: 10.1172/jci.insight.167082. 
  3. Yan et al. _JCI Insight_. 2025. doi: 10.1172/jci.insight.197579.

## Scripts
### 01_analysis.R
Workflow:
  1. Load datasets
  2. Perform QC and preprocessing
  3. Integrate datasets
  4. Identify myeloid/macrophage population
  6. Generate integrated objects
  7. Annotate myeloid/macrophage subsets

### 02_figures.R
Workflow:
  1. Generate UMAP for BAL samples
  2. Generate dot plot for BAL samples
  3. Generate feature plots for BAL samples
  4. Generate UMAP for lung samples, split by study
  5. Generate dot plot for lung samples
  6. Generate feature plots for lung samples, split by study

## Software

Analysis performed in R 4.4.2.
Package versions are reported in sessioninfo.txt
