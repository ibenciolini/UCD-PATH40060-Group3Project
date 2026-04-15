# Proteomics-Based Patient Stratification of Pediatric Brain Tumours

## Overview

This repository loosely reproduces and extends the proteomics-based patient stratification
from [Petralia et al. (Cell, 2020)](https://doi.org/10.1016/j.cell.2020.10.044) using
the CPTAC/CHOP paediatric brain cancer dataset, publicly available on
[cBioPortal](https://www.cbioportal.org/study/summary?id=brain_cptac_2020).

**Dataset at a glance**

| Feature | Value |
|---|---|
| Samples | 218 tumours from 199 patients |
| Histological subtypes | 7 (LGG, Ependymoma, HGG, Medulloblastoma, Ganglioglioma, Craniopharyngioma, ATRT) |
| Data modalities | WGS · RNA-seq · global proteomics · phosphoproteomics |
| Survival events | 40 deaths / 192 patients with OS data |

**Why proteomics for stratification?**  
Unlike sparse binary mutation profiles, proteomics data is already a dense,
continuous, pathway-integrated signal that directly captures the functional downstream
consequences of mutations, splicing, and post-translational regulation, essentially
the "smoothed" output that Network-Based Stratification (NBS) tries to recover from
mutation data alone.

**Analysis pipeline**

```
Local data files
    │  load + clean column names (janitor)
    V
Quality control
    │  remove high-missing proteins (> 30%), KNN imputation, PCA sanity check
    V
Dimensionality reduction
    │  top 50 PCs → shift to non-negative for NMF input
    V
k selection
    │  cophenetic correlation · dispersion · silhouette  (k = 2–8)
    V
Consensus NMF clustering
    │  60 NMF runs per k · consensus matrix · final cluster labels
    V
Biological validation
    ├── Kaplan–Meier survival curves (log-rank test, Cox PH model)
    ├── Histological composition per cluster
    └── Differential protein expression (limma one-vs-rest)
```

An optional NBS comparison arm (Section 9) runs mutation-based network
propagation (STRING v12, α = 0.7, igraph) and computes the Adjusted Rand Index
between proteomics- and mutation-derived subtypes to quantify how much biology
proteomics captures beyond genomics.

## Repository Structure

```
.
├── Analysis_BrainCPTAC2020.Rmd        # Main analysis notebook
├── Analysis_BrainCPTAC2020.html       # Pre-rendered HTML output
├── brain_cptac_2020/                  # Data files (downloaded from cBioPortal)
│   ├── data_clinical_patient.txt
│   ├── data_clinical_sample.txt
│   ├── data_protein_quantification.txt
│   └── data_mutations.txt
├── Cell_Petralia2020.pdf                # Petralia et al. Cell 2020
└── README.md
```

---

## Reproducing the Analysis

### 1. Get the data

Download the study from cBioPortal and place the files in `brain_cptac_2020/`:

```bash
# Option A: manual download
# Go to https://www.cbioportal.org/study/summary?id=brain_cptac_2020
# Click "Download" → download the full study zip → extract into brain_cptac_2020/

# Option B: command line (requires curl)
curl -L "https://cbioportal-datahub.s3.amazonaws.com/brain_cptac_2020.tar.gz" \
     -o brain_cptac_2020.tar.gz
tar -xzf brain_cptac_2020.tar.gz
```

Four files are required:

| File | Contents |
|---|---|
| `data_clinical_patient.txt` | OS months/status, age, sex (199 rows) |
| `data_clinical_sample.txt` | Histological subtype per sample (218 rows) |
| `data_protein_quantification.txt` | Global proteomics z-scores (6 429 proteins × 218 samples) |
| `data_mutations.txt` | Somatic mutation calls (9 951 rows) |

### 2. Install R dependencies

```r
# CRAN packages
install.packages(c(
  "tidyverse", "janitor", # data wrangling
  "NMF", "cluster", "pheatmap", # clustering & visualisation
  "RColorBrewer", # colour palettes
  "survival", "survminer", # Kaplan–Meier / Cox
  "limma" # differential protein expression
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "impute", # KNN imputation for missing proteomics values
  "limma", # also on Bioconductor
  "STRINGdb", # for NBS comparison arm
  "igraph"
))
```

### 3. Run the notebook

Open `Analysis_BrainCPTAC2020.Rmd` in RStudio and run all code chunks or knit the document.

---

## Key Results

| Finding | Detail |
|---|---|
| Optimal clusters | k = 8 (highest cophenetic correlation, matching Petralia et al.) |
| Survival separation | log-rank p < 0.0001 across 8 clusters |
| Cox concordance | C = 0.814 (cluster + age + sex model) |
| Worst-prognosis cluster | C5 — HR ≈ 10.5 vs C1; enriched for HGG, marked by IGFBP2/SHMT2/PBK |
| Cross-histology grouping | LGG, Ependymoma, and Ganglioglioma co-cluster in C1/C3/C4 |
| Therapeutic target | ERBB2 upregulated in C4 (ciliopathy/developmental programme) |

---

## Reference

Petralia F. et al. (2020). *Integrated Proteogenomic Characterization across Major
Histological Types of Pediatric Brain Cancer.*
**Cell** 183(7): 1962–1985.e31.
[doi:10.1016/j.cell.2020.10.044](https://doi.org/10.1016/j.cell.2020.10.044)

Additional methods references:

- Hofree M. et al. (2013). Network-based stratification of tumour mutations. **Nature Methods** 10, 1108–1115.
- Monti S. et al. (2003). Consensus Clustering. **Machine Learning** 52, 91–118.
- Cerami E. et al. (2012). The cBio Cancer Genomics Portal. **Cancer Discovery** 2, 401–404.
