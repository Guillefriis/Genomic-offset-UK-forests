# Genomic Offset Pipeline for UK Forest Trees

R workflow for genomic offset analyses under present and future climate scenarios.

This repository contains the full analytical pipeline used to:

- Process SNP allele frequencies  
- Select environmental predictors  
- Fit genotype–environment models  
- Detect candidate adaptive loci  
- Compute genomic offset metrics  

**Project:** Friis et al. 2026 (in preparation)

---

## Workflow Overview

The pipeline is modular and sequential.  
Scripts must be executed in order.

---

## 1. Climate Data Preparation

### `R_Offset1a_ClimData_Local_v3.1.R`

Preparation and formatting of climatic predictor variables.

### `R_Offset1b_SNPfreq_Cluster_v2.1.R`

Computation of population-level SNP allele frequencies.  
Designed for execution on a compute cluster.

---

## 2. Frequency Matrix Formatting

### `R_Offset2_FormatFreqs_Cluster_v2.1.R`

Reshaping and formatting allele frequency matrices for downstream modelling.

---

## 3. Environmental Variable Selection

### `R_Offset3_TestAllVar_Cluster_v2.2.R`

Iterative testing of environmental predictors.

### `R_Offset5_TestVIF_Cluster_v2.1.R`

Variance Inflation Factor (VIF) analysis to reduce multicollinearity among predictors.

---

## 4. Environmental Model Construction

### `R_Offset4_BuiltENM_Local_v3.1.R`

Construction of environmental niche models.

### `R_Offset6_DetectOutliers_Cluster_v2.2.R`

Identification of candidate adaptive loci.

---

## 5. Genomic Offset Estimation

### `R_Offset7_ModelsFreq_Local_v2.3.R`

Prediction of allele frequencies under future climate scenarios and computation of genomic offset metrics.

---

## Computational Environment

- R ≥ 4.2 recommended  
- Designed for mixed local and HPC cluster execution  

Typical R packages used across scripts include:

- tidyverse  
- vegan  
- terra  
- sf  
- AlleleShift  
- car  
- qvalue  

---

## Data Availability

Raw genomic datasets, raster layers, and intermediate objects are not included due to file size and data governance constraints.

Scripts are provided to ensure transparency and reproducibility of analytical procedures.

---
