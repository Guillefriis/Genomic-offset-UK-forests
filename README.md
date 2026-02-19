# Genomic Offset Pipeline for UK Forest Trees

R workflow for genomic offset analyses under present and future climate scenarios.

This repository contains scripts to:
- compute population SNP allele frequencies
- prepare climate predictors
- perform RDA-based environmental association analyses
- detect candidate adaptive loci
- estimate genomic offset and sapling maladaptation

Project: Friis et al. 2026 (in preparation)

## Workflow

Scripts are modular and must be run sequentially.

### Climate and Frequency Preparation
- `R_Offset1a_ClimData_Local_v3.1.R`
- `R_Offset1b_SNPfreq_Cluster_v2.1.R`
- `R_Offset2_FormatFreqs_Cluster_v2.1.R`

### Environmental Screening and Model Fitting
- `R_Offset3_TestAllVar_Cluster_v2.2.R`
- `R_Offset4_BuiltENM_Local_v3.1.R`
- `R_Offset5_TestVIF_Cluster_v2.1.R`

### Candidate Loci and Genomic Offset
- `R_Offset6_DetectOutliers_Cluster_v2.2.R`
- `R_Offset7_ModelsFreq_Local_v2.3.R`

The final stage includes:
- RDA-based candidate SNP identification
- Projection of RDA loadings onto present and future climate rasters
- Computation of genomic offset metrics
- Comparison of observed versus predicted sapling allele frequencies

## Methodological Basis

Several components of this workflow are adapted from and inspired by scripts developed by Thibaut Capblancq, particularly:

- RDA landscape genomics framework  
  https://github.com/Capblancq/RDA-landscape-genomics

- RDA genome scan approach  
  https://github.com/Capblancq/RDA-genome-scan

These repositories provide the conceptual and statistical foundations for the RDA-based candidate detection and genome–environment association procedures implemented here.

## Requirements

R ≥ 4.2 recommended.  
Common packages include tidyverse, vegan, terra, sf, car, qvalue, and related dependencies.

Cluster scripts assume high-memory parallel execution.

## Data

Large genomic datasets and raster layers are not included.

## License

MIT License.
