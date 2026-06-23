# Prenatal Environmental Exposure Drives Distinct MicroRNA Patterns in Maternal and Cord Serum: the NEHO Cohort

Code repository accompanying the manuscript:

**Longo V., Cosentini I., Li Vigni A., Aloi N., Sampino A., Contarino F., Colombo P., Ruggieri S., Drago G.**

*Prenatal Environmental Exposure Drives Distinct MicroRNA Patterns in Maternal and Cord Serum: the NEHO Cohort*

DOI: https://doi.org/10.21203/rs.3.rs-8107123/v1

## Study Background

Environmental exposures during pregnancy can influence maternal–foetal communication and developmental programming through epigenetic mechanisms. Among these, circulating microRNAs (miRNAs) represent promising biomarkers of both environmental exposure and biological response.

This study investigated the relationship between prenatal exposure to essential elements (EEs) and persistent organic pollutants (POPs) and circulating miRNA expression profiles in matched maternal and cord serum samples from the NEHO (Neonatal Environment and Health Outcomes) birth cohort. The analytical framework combined single-pollutant regression models, Bayesian Weighted Quantile Sum (bWQS) mixture analyses, and pathway enrichment approaches to identify exposure-associated molecular signatures in maternal and foetal compartments.

This repository provides the code used to preprocess data, perform statistical analyses, and reproduce the figures presented in the manuscript.

---

## Repository Overview

This repository contains the R scripts used to preprocess data, perform statistical analyses, and generate the figures presented in the manuscript.

Because participant-level data from the NEHO cohort cannot be publicly shared for privacy and ethical reasons, a synthetic dataset can be generated to reproduce the analytical workflow.

---

## Repository Structure

```text
.
├── scripts/
│   ├── 0_create_synthetic_data.R
│   ├── 1_0_miRNA_normalization_and_dataset_preparation.R
│   ├── 2_0_linear_regression.R
│   ├── 3_0_bwqs_analysis.R
│   └── figures/
│       └── Scripts used to generate the figures reported in the manuscript
├── supplementary_data/
│
└── README.md
```

### Scripts

| Script                                            | Description                                                                  |
| ------------------------------------------------- | ---------------------------------------------------------------------------- |
| `0_create_synthetic_data.R`                       | Generates a synthetic dataset preserving the structure of the original data. |
| `1_miRNA_normalization_and_dataset_preparation.R` | Performs miRNA normalization and prepares the analytical dataset.            |
| `2_linear_regression.R`                           | Runs single-pollutant regression analyses.                                   |
| `3_bwqs_analysis.R`                               | Runs Bayesian Weighted Quantile Sum (bWQS) mixture analyses.                 |
| `figures/`                                        | Scripts used to generate manuscript figures.                                 |

---

## Reproducing the Analysis

Run the scripts in the following order:

```r
source("scripts/0_create_synthetic_data.R")
source("scripts/1_0_miRNA_normalization_and_dataset_preparation.R")
source("scripts/2_0_linear_regression.R")
source("scripts/3_0_bwqs_analysis.R")
```

Figure-generation scripts can then be executed individually from the `figures/` directory.

---

## Data Availability

The original data used in this study are not publicly available due to ethical and privacy restrictions.

To facilitate reproducibility, this repository includes code for generating a synthetic dataset that mimics the structure and statistical properties of the original data without containing identifiable participant information.

---

## Supplementary Material

The `supplementary_data/` directory contains supplementary tables and additional material associated with the manuscript.

---

## Software Requirements

Analyses were performed using R.

Required packages include:

* tidyverse
* bwqs
* pheatmap
* ggplot2
* clusterProfiler
* multiMiR

Additional dependencies are specified within individual scripts.

---

## Citation

If you use this repository, please cite:

Longo V., Cosentini I., Li Vigni A., Aloi N., Sampino A., Contarino F., Colombo P., Ruggieri S., Drago G.

*Prenatal Environmental Exposure Drives Distinct MicroRNA Patterns in Maternal and Cord Serum: the NEHO Cohort.*

Research Square (2025).

DOI: https://doi.org/10.21203/rs.3.rs-8107123/v1
