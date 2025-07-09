# README


## Overview

This repository contains code and data for a two-part analysis combining structural and functional brain metrics of perisylvian sub-regions (PSRs) using multivariate statistical methods. The two different types of functional metrics are as follows:

**1.	fractional amplitude of low-frequency fluctuations (fALFF)**

**2.	Functional connectivity (FC)**

## Folder Structure

```
bash

project-root/
│
├── indiv_features/
│   └── regout/
│       ├── morphologicalfeats.csv        # Structural metrics (age-regressed-out)
│       ├── Falff.csv                     # Functional metric (fALFF, age-regressed-out)
│       └── FC_XX_to_all.csv              # FC metrics (for PSRs 11 to 86, age-regressed-out)
│
├── FALFF/
│   ├── main_FALFF_PSR.R                  # Main script for Experiment 1 (fALFF) in R
│   ├── funcs_FALFF_PSR.R                 # Supporting functions in R
│   ├── feature_selection/                # Results of feature selection step (SCCA)
│   │   ├── permsinfo_com_1.csv                   # Results of selecting sparsity levels 
│   │   ├── sigcomponents_permutation.csv         # Result of selecting significant components of SCCA 
│   │   ├── allmain/                              # SCCA results on original data
│   │   ├── allboot/                              # Results of bootstrapping test for significant components of SCCA
│   │   └── allvisig/
│   │       └── selected_features.csv                  # Aggregated weights for each metric using all significant components 
│   └── CS-CF/
│       ├── FC_XX/                         # Results of CS-CF feature extraction (CCA)
│           ├── maincorr.csv                      # Canonical correlation coefficients on original data
│           ├── cca_weights_struc.csv             # Structural canonical weights on original data
│           ├── cca_weights_func.csv              # Functional canonical weights on original data
│           └── (heatmaps, histograms, regression plots, etc.)
│
├── FC/
│   ├── main_FC.R                         # Main script for Experiment 2 (FC) in R
│   ├── funcs_FC.R                        # Supporting functions in R
│   ├── feature_selection/                # Results of feature selection step (SCCA)
│   │   ├── permsinfo_seeds.csv                   # Results of selecting sparsity levels for all PSRs as seed
│   │   ├── sigcomponents_permutation_*.csv       # Result of selecting significant components of SCCA for each PSR
│   │   ├── main/                                 # SCCA results on original data
│   │   ├── boot/                                 # Results of bootstrapping test for significant components of SCCA
│   │   └── visig/
│   │       ├── FC_XX/
│   │           ├── rois_using_all_components.csv       # Selected regions/metrics
│   │           └── weights_using_all_components.csv    # Aggregated weights for each metric
│   └── CS-CF/
│       ├── FC_XX/                           # Results of CS-CF feature extraction (CCA)
│           ├── maincorr_XX.csv                  # Canonical correlation coefficients on original data
│           ├── struc_weightsXX.csv              # Structural canonical weights on original data
│           ├── func_weightsXX.csv               # Functional canonical weights on original data
│           └── (heatmaps, histograms, regression plots, etc.)
│
└── README
```

## How to Run the Code

### Dependencies

Install required R packages using: 
```
install.packages(c("foreach", "iterators", "PMA", "parallel", "doParallel",
                   "plyr", "ggplot2", "circlize", "ComplexHeatmap",
                   "yacca", "broom", "reshape2", "reporter"))
```

### Experiment 1: Using fALFF as Functional Metric

1. Navigate to the `FALFF/` folder.
2. Open the file `main_FALFF_PSR.R`.
3. On line 11, update the `maindirectory` variable to your project_root directory.

```
r
maindirectory <- "/path/to/your/project-root/"
```
### Experiment 2: Using FC values as Functional Metric

1. Navigate to the `FC/` folder.
2. Open the file `main_FC.R`.
3. On line 11, update the `maindirectory` variable to your project_root directory.

```
r
maindirectory <- "/path/to/your/project-root/"
```

## Output Description
The details of output files for each experiment in corresponding directory.

### Feature Selection (`feature_selection/`)
* **permsinfo_*.csv**: Optimal sparsity level from permutation tests (output of `Permutation` section in code).
* **sigcomponents_permutation_*.csv**: IDs of significant components (output of `checking significance of other canonical components (permutation)` section in code).
* **main/**: Sparse CCA (SCCA) results with selected sparsity level (output of `main` section in code).
* **boot/**: Bootstrapping results for significance each component of SCCA (output of `applying SCCA_bootstrapping` section in code). 
* **visig/rois_using_all_components.csv**: Final selected regions/metrics which are canonical weights with npn-zero CI (output of `visualization(feature and roi selection)` section in code). 
* **visig/weights_using_all_components.csv**: sum of the absolute canonical weights for all significant components for each metric separately (output of `visualization(feature and roi selection)` section in code).  

### CS-CF extraction feature (`CS-CF/`)
* **maincorr_*.csv**: Canonical correlation coefficients on original data.
* **struc_weight_*.csv**: Canonical weights for selected metrics using CCA for structural features.
* **func_weight_*.csv**: Canonical weights for selected metrics using CCA for functional features.
* **scatter_cca_variates.png**: Regression analysis between CS and CF features
* **hist_original_CCC.tiff**: Comparing between the correlation of individual S-F metrics (histogram) and combined S_F features (vertical red line)
* **corr_heatmap_original.tiff**: Heatmaps of correlations between individual structural anf functional metrics
