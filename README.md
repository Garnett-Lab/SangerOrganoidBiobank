# 🧬 CRISPR Organoid Dependency Analysis

> Complete computational pipeline for genome-wide CRISPR screening in patient-derived organoids

---

## 🎯 Overview

Comprehensive R pipeline for analyzing CRISPR dependency screens in organoid models, identifying essential genes, biomarkers, and therapeutic vulnerabilities.

**Publication:** A tumour-derived organoid biobank defines a functional map of gene dependencies | *Nature* (2025)

---

## 🚀 Quick Start

### Installation

```r
# Install CRAN packages
install.packages(c("data.table", "tidyverse", "scales", "lmtest", "digest"))

# Install Bioconductor packages
BiocManager::install(c("CoRe", "msigdbr", "clusterProfiler", "CRISPRcleanR", "sva"))
```

### Run Pipeline

#### 1. Process raw CRISPR data
source("scripts/CRISPR_pipeline.R")
source(run_BAGEL2_bf.sh)
source(run_BAGEL2_pr.sh)

        CRISPR_pipeline.R (R)
            ↓
        Normalized fold changes (TSV files)
            ↓
        run_BAGEL2_bf.sh (Bash) ← YOU ARE HERE
            ↓
        BAGEL2 Bayes Factors (.bf files)
            ↓
        BAGEL2 Precision-Recall (.pr files)
            ↓
        CRISPR_pipeline.R (R) - Import results
            ↓
        Binary essentiality calls

#### 2. Identify core genes
source("scripts/core_fitness_genes_analysis.R")

#### 3. Find differential dependencies
source("scripts/CRISPR_differential_dependency_analysis.R")

#### 4. Prepare features
source("scripts/input_biomarker_analysis.R")

#### 5. Discover biomarkers
source("scripts/CRISPR_biomarker_analysis.R")

####6. Reproduce main figure of the paper
surce(Main_figures.Rmd)


---

## 🔬 Key Methods

### CRISPR Processing
- **Normalization**: CRISPRcleanR → Spline → ComBat
- **QC**: Replicate distance, ROC curves, Cohen's D
- **Calling**: BAGEL2 (FDR < 0.05)

### Biomarker Discovery
- **Model**: Linear regression with likelihood ratio test
- **Features**: Mutations, CNVs, expression, clinical
- **Classification**: A (strong), B (moderate), C (significant)
- **Thresholds**: FDR < 0.05, |Δmean| > 0.5

---

## 📥 Input Data

Required data types:
- ✅ Raw CRISPR counts (sgRNA level) - Not provided for this project
- ✅ Somatic mutations (VCF-derived)
- ✅ Copy number alterations
- ✅ RNA-seq expression (TPM)
- ✅ Clinical metadata
- ✅ Control gene lists

- ** All imput data provided in:** figshare

---

## 📤 Output Files

### Supplementary Tables

| Table | Description |
|-------|-------------|
| S4.2 | Core fitness genes in organoids |
| S5 | Binary essentiality calls |
| S6 | Normalized LFCs (fitness scores) |
| S7.1-3 | Differential dependencies by cancer type |
| S8.1-4 | Biomarkers and priority scores |

---

## 🔄 Reproducibility

### Session Information
- session_info_R.txt


**Key Methods:**
- ADAM: Behan et al., Nature 2019
- CRISPRcleanR: Iorio et al., BMC Genomics 2018
- BAGEL2: Kim & Hart, Genome Medicine 2021
- ComBat: Johnson et al., Biostatistics 2007

---

## 📧 Contact

**Questions?** Open an [issue](https://github.com/Garnett-Lab/SangerOrganoidBiobank) or email: ch24@sanger.ac.uk

---

**Last Updated:** November 2025 | **Status:** Reviewed
