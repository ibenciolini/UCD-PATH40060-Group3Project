# Proteomics-Based Patient Stratification and Therapeutic Target Discovery in Pediatric Brain Tumours

**Dataset:** CPTAC/CHOP Brain Cancer Cohort (`brain_cptac_2020`) [cBioPortal](https://www.cbioportal.org/study/summary?id=brain_cptac_2020)

**Reference:** Petralia F. et al. (2020) "Integrated Proteogenomic Characterization across Major Histological Types of Pediatric Brain Cancer." *Cell* 183(7), 1962–1985.e31.

---

## Overview

This repository contains a full computational analysis of the publicly available CPTAC/CHOP pediatric brain tumour dataset (218 tumour samples, 199 patients, 7 histological subtypes). Building on the unsupervised NMF proteomic subtypes from the original paper, the analysis pursues two precision-oncology questions:

1. **Therapeutic target discovery (§14).** For each proteomic subtype, identify candidate drug targets supported by: cluster-specific protein overexpression (limma), inferred kinase activity from phosphoproteomics (KSEA + OmniPath), brain-cancer cell-line dependency (DepMap CRISPR essentiality), and approved-drug coverage (OpenTargets/DGIdb). The deliverable is an integrated, evidence-scored per-cluster therapeutic panel.

2. **Molecular surrogate classification (§15).** Can the proteome predict clinically actionable molecular markers (WHO grade and canonical driver mutation status) that currently require sequencing? With ≥100 patients per outcome class, this is the only branch where a properly held-out predictive claim is statistically defensible in this cohort. A cross-validated lasso classifier (cv.glmnet) serves as the modelling backbone.

Survival analysis (§10–11) is included as a **descriptive** characterisation of the three proteomic subtypes; the cohort has only ~36 OS events, which is insufficient for a credible per-patient predictive risk model.

---

## Repository Structure

```
.
├── brain_cptac_2020/ # Raw data downloaded from cBioPortal
│   ├── data_clinical_patient.txt # Demographics, OS/DFS, WHO grade, driver-mutation flags
│   ├── data_clinical_sample.txt # Sample metadata, treatment, tissue site
│   ├── data_protein_quantification.txt # Global proteomics matrix
│   └── data_phosphoprotein_quantification.txt # Phosphosite-level abundances
│
├── data/
│   ├── Model.csv # DepMap cell-line metadata (lineage, model type)
│   ├── clinical_target_compact.csv # OpenTargets clinical-stage druggability annotations
│   └── omnipath_enzsub_human.csv # OmniPath kinase–substrate relationships (KSEA input)
│
├── figures/ # All PNG figures auto-saved by knitr (fig.path)
│
├── Analysis_BrainCPTAC2020_v3.Rmd # Main analysis notebook (source of truth)
├── Analysis_BrainCPTAC2020_v3.md # GitHub-rendered output (github_document)
└── Analysis_BrainCPTAC2020_v3.pdf # PDF output: results and figures only, no code
```

---

## Analysis Pipeline

```
Raw cBioPortal flat files
    │ load + clean column names (janitor::clean_names)
    ▼
Quality control (§3)
    │ remove proteins missing in >30% of samples; KNN imputation
    ▼
Exploratory analysis (§4)
    │ PCA, PLS-DA, sPLS-DA, UMAP
    ▼
Dimensionality reduction (§5)
    │ top 25 PCs; shift to non-negative for NMF
    ▼
k selection (§6)
    │ cophenetic / dispersion / silhouette / cluster purity across k = 2..8
    ▼
Consensus NMF clustering (§7)
    │ 20 runs at k = 3 (brunet); consensus matrix; final cluster labels
    ▼
Biological characterisation (§8–13)
    ├── Histological composition per cluster (§9)
    ├── Kaplan–Meier survival + log-rank test (§10)
    ├── Cox model: cluster + extended clinical baseline (§11)
    ├── Differential protein expression: limma one-vs-rest (§12)
    └── Variable importance: random forest (§13)
    ▼
Therapeutic target discovery (§14)
    ├── Phospho one-vs-rest contrasts per cluster
    ├── Kinase activity inference (KSEA, OmniPath substrates)
    ├── DepMap brain-line CRISPR dependency for cluster markers
    ├── Druggability annotation (OpenTargets / DGIdb)
    └── Integrated per-cluster therapeutic panel (panel_score)
    ▼
Molecular surrogate classifier (§15)
    └── Cross-validated lasso (cv.glmnet, keep=TRUE) for
        WHO grade and any-driver-mutation prediction
```

---

## Proteomic Subtypes

Three subtypes are identified at k = 3 (cophenetic, silhouette, and mean cluster purity all peak here):

| Subtype | Biological identity | Key marker proteins |
|---|---|---|
| **C1 — Proliferative/Nuclear** | High cell-cycle and nuclear activity | SEPT8, PCNA, NPM1, SON |
| **C2 — Neuronal/Synaptic** | Neuronal differentiation programme | NEFL, NEFM, CAMK2B, STX1B |
| **C3 — Mesenchymal/Microenvironment** | Extracellular matrix and immune infiltration | CSPG4, LRP1, IL1RAP, NLGN3 |

All three subtypes span multiple histological diagnoses, replicating the central finding of Petralia et al. that molecular programmes cut across classical pathological boundaries.

---

## Key Outputs

- **Consensus heatmap**: sample co-clustering stability at k = 3.
- **KM curves**: overall and disease-free survival per proteomic subtype (descriptive; log-rank p-value reported).
- **Differential expression tables**: limma one-vs-rest results per cluster (logFC, adjusted p-value).
- **KSEA results**: per-cluster kinase activity z-scores inferred from phosphoproteomics.
- **Therapeutic panel table**: integrated target list scored by `panel_score` (limma logFC + DepMap essentiality + druggability phase).
- **ROC curve**: out-of-fold held-out AUC for the WHO grade lasso classifier.

All figures are written to `./figures/` with the R chunk label as the filename.

---

## Reproducing the Analysis

### Requirements

R ≥ 4.3. The following packages are required:

**CRAN:** `tidyverse`, `janitor`, `car`, `gridExtra`, `cluster`, `pheatmap`, `RColorBrewer`, `ConsensusClusterPlus`, `umap`, `survival`, `survminer`, `limma` (via Bioconductor), `randomForest`, `randomForestSRC`, `glmnet`, `conflicted`

**Bioconductor:** `impute`, `NMF`, `mixOmics`

**Optional** (chunks degrade if absent): `OmnipathR`, `depmap`, `httr2`, `jsonlite`, `GSVA`, `msigdbr`

### Data

Download the study from cBioPortal and place the four files listed above in `./brain_cptac_2020/`. The supplementary data files in `./data/` are included in this repository.

---

## Key Configuration Parameters

All parameters are set in the `config` chunk at the top of the Rmd and documented in-line:

| Parameter | Default | Description |
|---|---|---|
| `MISSING_THRESHOLD` | 0.30 | Remove proteins missing in >30% of samples |
| `N_PCS` | 25 | PCA components fed into NMF |
| `K_MIN` / `K_MAX` | 2 / 8 | k range for the NMF sweep |
| `NMF_RUNS` | 20 | NMF initialisations per k (paper used 60) |
| `KSEA_MIN_SUBSTRATES` | 5 | Minimum substrates for a valid KSEA z-score |
| `DEPMAP_ESS_THRESH` | −0.5 | CERES/Chronos threshold for "essential" |
| `DRUGGABILITY_PHASE` | 1 | Minimum OpenTargets clinical phase |
| `CLASSIFIER_FOLDS` | 10 | Folds for cv.glmnet cross-validation |

---

## Important Caveats

- **Survival analysis is descriptive only.** ~36 OS events across 218 patients is well below the ≥10-events-per-parameter rule of thumb for a credible predictive model. KM curves and the Cox model characterise the subtypes; they do not generate a per-patient risk score.
- **The driver-mutation classifier should be interpreted cautiously.** Only 19/199 patients have a recorded driver mutation, a ~10:1 class imbalance. A CV AUC near 1.0 in this setting likely reflects the classifier learning histological subtype (which correlates with driver status) rather than driver status directly.
- **Single-cohort validation.** All results are internal. External validation in an independent cohort (PBTA, OpenPedCan, or an adult CPTAC release) is the required next step before any clinical claim.
- **The therapeutic panel is hypothesis-generating.** Targets in §14 are supported by convergent computational evidence, not experimental validation.

---

## References

- Petralia F. et al. (2020) *Cell* 183(7), 1962–1985.e31. doi:10.1016/j.cell.2020.10.044
- Monti S. et al. (2003) *Machine Learning* 52, 91–118. — Consensus clustering methodology
- Casado P. et al. (2013) *Science Signaling* 6(268), rs6. — KSEA methodology
- Türei D. et al. (2016) *Nature Methods* 13, 966–967. — OmniPath kinase–substrate annotations
- Tsherniak A. et al. (2017) *Cell* 170(3), 564–576. — DepMap CRISPR essentiality
- Ochoa D. et al. (2023) *Nucleic Acids Research* 51(D1), D1353–D1359. — Open Targets Platform
- Friedman J. et al. (2010) *Journal of Statistical Software* 33(1), 1–22. — glmnet
