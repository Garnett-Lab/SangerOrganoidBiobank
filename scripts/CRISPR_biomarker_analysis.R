################################################################################
# Script: CRISPR_biomarker_analysis.R
# Purpose: Discover genomic and clinical biomarkers of gene dependencies
# Description: Performs comprehensive biomarker analysis using linear modeling
#              to identify genomic features, clinical variables, and expression
#              patterns that predict gene essentiality in organoids
#
# Workflow:
# 1. Load multi-omic features (mutations, CNVs, expression, clinical)
# 2. Load CRISPR dependency data (LFCs from pipeline)
# 3. Filter features by sample size and variance
# 4. Perform linear regression (LFC ~ feature + covariates)
# 5. Apply likelihood ratio test and FDR correction
# 6. Classify biomarkers by effect size (Classes A, B, C)
# 7. Calculate priority scores for each gene-biomarker pair
# 8. Generate final biomarker tables for publication
#
# Statistical Framework:
# - Linear Model: LFC ~ feature + ploidy + library [+ MSI for colorectal]
# - Test: Likelihood Ratio Test (LRT) comparing models with/without feature
# - Correction: Benjamini-Hochberg FDR within each feature type
# - Effect Size: Cohen's coefficient, mean difference, correlation
################################################################################

# Load required libraries
library(data.table)
library(tidyverse)
library(digest)      # For duplicate feature detection
library(scales)      # For feature scaling
library(lmtest)      # For likelihood ratio tests
library(dplyr)

################################################################################
# 1. LOAD CLINICAL AND MOLECULAR FEATURES
################################################################################

## Clinical features

# Tumor stage (individual stages I-IV)
stages1 <- as.data.frame(
  fread('/path_to_folder/data_binary-stages_ALL-organoids.csv')
)

# Tumor stage (early vs advanced)
stages2 <- as.data.frame(
  fread('/path_to_folder/data_binary-stages2_ALL-organoids.csv')
)

# Patient age (young vs old)
age <- as.data.frame(
  fread('/path_to_folder/data_binary-age_ALL-organoids.csv')
)

# Patient sex
gender <- as.data.frame(
  fread('/path_to_folder/data_binary-sex_ALL-organoids.csv')
)

# Colorectal tumor location (proximal, distal, rectum)
localization_COADREAD <- as.data.frame(
  fread('/path_to_folder/data_binary-sides_COLO-organoids.csv')
)

# Metastasis status (colorectal only)
metastasis_COADREAD <- as.data.frame(
  fread('/path_to_folder/data_binary-metastasis_COLO-organoids.csv')
)

# Barrett's esophagus diagnosis (esophageal only)
barrett_ESCA <- as.data.frame(
  fread('/path_to_folder/data_binary-barrett_OESO-organoids.csv')
)

## Molecular subtypes (colorectal only)

# CMS (Consensus Molecular Subtypes)
CMS_COADREAD <- as.data.frame(
  fread('/path_to_folder/data_CMS-subtypes_organoids.csv')
)

# CRIS (Cancer-related intrinsic subtypes)
CRIS_COADREAD <- as.data.frame(
  fread('/path_to_folder/data_CRIS-subtypes_organoids.csv')
)

## Merge all clinical features
clinical <- merge(stages1, stages2, by = "sample_ID")
clinical <- merge(clinical, age, by = "sample_ID", all.x = TRUE)
clinical <- merge(clinical, gender, by = "sample_ID", all.x = TRUE)
clinical <- merge(clinical, localization_COADREAD, by = "sample_ID")
clinical <- merge(clinical, metastasis_COADREAD, by = "sample_ID", all.x = TRUE)
clinical <- merge(clinical, barrett_ESCA, by = "sample_ID", all.x = TRUE)
clinical <- merge(clinical, CMS_COADREAD, by = "sample_ID", all.x = TRUE)
clinical <- merge(clinical, CRIS_COADREAD, by = "sample_ID", all.x = TRUE)

# Remove redundant stage and age columns 
clinical$stage_II <- NULL
clinical$stage_III <- NULL
clinical$old <- NULL
clinical$advanced_stage <- NULL



################################################################################
# 2. LOAD GENOMIC FEATURES
################################################################################

## Driver gene mutations (binary: 0/1)
gene_mutations <- as.data.frame(
  fread('/path_to_folder/data_driver-mutations_ALL-organoids.csv')
)

## Specific driver variants (e.g., KRAS G12D, BRAF V600E)
allele_mutations <- as.data.frame(
  fread('/path_to_folder/data_driver-variants_ALL-organoids.csv')
)

# Merge mutation data
mutations <- merge(gene_mutations, allele_mutations, by = "sample_ID", all = TRUE)

## Somatic copy number alterations (SCNAs) + MSI status
SCNA <- as.data.frame(
  fread('/path_to_folder/data_CNV-genes_ALL-organoids.csv')
)

## RNA-seq expression (log2 TPM)
RNASEQ <- as.data.frame(
  fread('/path_to_folder/data_RNAseq_ALL-organoids.csv')
)

## Mutational signatures (binary: signature present >5% of mutations)
signatures <- as.data.frame(
  fread('/path_to_folder/data_mutational-signatures_ALL-organoids.csv')
)

## Loss-of-Function (LoF) and Gain-of-Function (GoF) events
LOF <- as.data.frame(
  fread('/path_to_folder/data_LOF_ALL-organoids.csv')
)

GOF <- as.data.frame(
  fread('/path_to_folder/data_GOF_ALL-organoids.csv')
)

LoF_GoF <- merge(LOF, GOF, by = "sample_ID", all = TRUE)



################################################################################
# 3. LOAD COVARIATES
################################################################################

## Covariates: Variables to control for in all models
# These are NOT tested as biomarkers, but adjust for technical/biological confounders

covariates <- as.data.frame(
  fread('/path_to_folder/data_covariates_ALL-organoids.csv')
)

# Extract MSI status separately (will be used as covariate for colorectal)
MSI <- covariates[, c("sample_ID", "msStatus")]



################################################################################
# 4. MERGE ALL BINARY FEATURES
################################################################################

## Combine all binary/categorical features into one matrix
# This will be filtered later based on sample size requirements
MERGE_binary <- merge(mutations, SCNA, by = "sample_ID", all = TRUE)
MERGE_binary <- merge(MERGE_binary, LoF_GoF, by = "sample_ID", all = TRUE)
MERGE_binary <- merge(MERGE_binary, signatures, by = "sample_ID", all = TRUE)
MERGE_binary <- merge(MERGE_binary, clinical, by = "sample_ID", all = TRUE)

# Remove MSI from features (will be added as covariate)
MERGE_binary$msStatus <- NULL



################################################################################
# 5. LOAD CRISPR DEPENDENCY DATA (RESPONSES)
################################################################################

## Load gene-level LFCs from CRISPR pipeline (supplementary_table_6)
LFC <- as.data.frame(
  fread('/path_to_folder/supplementary_table_6.csv')
)

## Remove duplicate samples and outliers
# These samples were either screened with multiple libraries or failed QC
LFC_filtered <- LFC[!(LFC$sample_ID == "HCM-SANG-0267-D12" & LFC$library == "minLib"), ]
LFC_filtered <- LFC_filtered[!(LFC_filtered$sample_ID == "HCM-SANG-0268-C18" & LFC_filtered$library == "v1_1"), ]
LFC_filtered <- LFC_filtered[!(LFC_filtered$sample_ID == "HCM-SANG-0270-C20" & LFC_filtered$library == "minLib"), ]
LFC_filtered <- LFC_filtered[!(LFC_filtered$sample_ID == "HCM-SANG-0272-C20" & LFC_filtered$library == "v1_1"), ]
LFC_filtered <- LFC_filtered[!(LFC_filtered$sample_ID == "HCM-SANG-0273-C18" & LFC_filtered$library == "minLib"), ]
LFC_filtered <- LFC_filtered[!(LFC_filtered$sample_ID == "HCM-SANG-0276-C18" & LFC_filtered$library == "minLib"), ]
LFC_filtered <- LFC_filtered[!(LFC_filtered$sample_ID == "HCM-SANG-0282-C18" & LFC_filtered$library == "minLib"), ]
LFC_filtered <- LFC_filtered[!(LFC_filtered$sample_ID == "HCM-SANG-0529-C18" & LFC_filtered$library == "v1_1"), ]
LFC_filtered <- LFC_filtered[!(LFC_filtered$sample_ID == "HCM-SANG-0284-C18" & LFC_filtered$library == "minLib"), ]
LFC_filtered <- LFC_filtered[!(LFC_filtered$sample_ID == "HCM-SANG-0288-C18" & LFC_filtered$library == "v1_1"), ]
LFC_filtered <- LFC_filtered[!(LFC_filtered$sample_ID == "WTSI-COLO_278" & LFC_filtered$library == "minLib"), ]
LFC_filtered <- LFC_filtered[!(LFC_filtered$sample_ID == "WTSI-COLO_376" & LFC_filtered$library == "minLib"), ]
LFC_filtered <- LFC_filtered[!(LFC_filtered$sample_ID == "HCM-SANG-0291-C15" & LFC_filtered$library == "v1_1"), ]
LFC_filtered <- LFC_filtered[!(LFC_filtered$sample_ID == "HCM-SANG-0299-C15-A" & LFC_filtered$library == "minLib"), ]
LFC_filtered <- LFC_filtered[!(LFC_filtered$sample_ID == "HCM-SANG-0299-C15-B" & LFC_filtered$library == "minLib"), ]
LFC_filtered <- LFC_filtered[!(LFC_filtered$sample_ID == "HCM-SANG-0306-C15" & LFC_filtered$library == "minLib"), ]

# Extract library information
library <- unique(LFC_filtered[, c("sample_ID", "library")])
covariates <- merge(covariates, library, by = "sample_ID", all.x = TRUE)

# Keep only essential columns
LFC_filtered <- LFC_filtered[, c("sample_ID", "gene", "LFC")]



################################################################################
# 6. LOAD DIFFERENTIAL DEPENDENCY ANALYSIS RESULTS
################################################################################

## These results identify genes with variable dependencies
# Used to focus biomarker analysis on genes with heterogeneous effects

# Gastrointestinal cancers (all)
differential_dep_gastrointestinal <- as.data.frame(
  fread("/path_to_folder/supplementary_table_7.1.csv")
)
differential_dep_gastrointestinal$ct <- "Gastrointestinal"

# Colorectal cancers
differential_dep_COADREAD <- as.data.frame(
  fread("/path_to_folder/supplementary_table_7.2.csv")
)
differential_dep_COADREAD$ct <- "Colorectal"

# Esophageal cancers
differential_dep_ESCA <- as.data.frame(
  fread("/path_to_folder/supplementary_table_7.3.csv")
)
differential_dep_ESCA$ct <- "Oesophageal"

# Combine all differential dependency results
differential_dep <- rbind(
  differential_dep_gastrointestinal, 
  differential_dep_COADREAD, 
  differential_dep_ESCA
)



################################################################################
# 7. LOAD CANCER TYPE ANNOTATIONS
################################################################################

cancer_type <- as.data.frame(
  fread('/path_to_folder/supplementary_table_2.2.csv')
)
cancer_type <- cancer_type[, c(1, 3)]

# Create gastrointestinal cancer category
cancer_type2 <- cancer_type[
  cancer_type$primary_tumour_type %in% c("Colorectal", "Oesophageal", "Gastric", "Pancreatic"), 
]
cancer_type2$primary_tumour_type <- "Gastrointestinal"

cancer_type <- rbind(cancer_type, cancer_type2)



################################################################################
# 8. LINEAR REGRESSION BIOMARKER ANALYSIS
################################################################################

## Analysis will be performed separately for each cancer type
# Focus on cancer types with sufficient sample sizes
TYPES <- c("Colorectal", "Oesophageal", "Gastrointestinal")

# Initialize storage for results
linear_model_results_RESPONSES <- NULL

## Main analysis loop: Iterate through cancer types
for (ct in TYPES) {
  
  print(paste("Processing cancer type:", ct))
  
  ##########################################################
  ## 8a. Filter data for current cancer type
  ##########################################################
  LFC_filtered_ct <- LFC_filtered[
    LFC_filtered$sample_ID %in% cancer_type[cancer_type$primary_tumour_type == ct, ]$sample_ID, 
  ]
  
  # Convert to wide format (genes as columns, samples as rows)
  LFC_filtered_ct <- spread(LFC_filtered_ct, "gene", "LFC")
  
  # Filter features and covariates for current cancer type
  MERGE_binary_ct <- MERGE_binary[MERGE_binary$sample_ID %in% LFC_filtered_ct$sample_ID, ]
  RNASEQ_ct <- RNASEQ[RNASEQ$sample_ID %in% LFC_filtered_ct$sample_ID, ]
  differential_dep_ct <- differential_dep[differential_dep$ct == ct, ]
  covariates_ct <- covariates[covariates$sample_ID %in% LFC_filtered_ct$sample_ID, ]
  
  # For esophageal, remove MSI covariate (insufficient variability)
  if (ct == "Oesophageal") {
    covariates_ct$msStatus <- NULL
  }
  
  ############################################################################
  ## 8b. Filter genes: Keep only differentially dependent genes
  ############################################################################
  # Requirement: Gene must be depleted in >=3 samples AND not depleted in >=3 samples
  # Rationale: Need sufficient variability to detect biomarker associations
  LFC_filtered_ct <- LFC_filtered_ct[
    , colnames(LFC_filtered_ct) %in% c(
      "sample_ID", 
      differential_dep_ct[
        differential_dep_ct$N_not_depleted > 2 & differential_dep_ct$N_depleted > 2, 
      ]$g
    )
  ]
  
  ############################################################################
  # 8c. FILTER BINARY FEATURES
  ############################################################################
  
  ## Goal: Remove features with insufficient sample size or redundancy
  
  if (!is.null(MERGE_binary_ct)) {
    
    print("Filtering binary features...")
    
    ## Filter 1: Remove features with insufficient samples in either group
    # Requirement: >= 3 samples with feature present AND >= 3 samples without
    # Rationale: Minimum sample size for statistical power and stable estimates
    binary_features <- c()
    
    for (c in colnames(MERGE_binary_ct)[2:ncol(MERGE_binary_ct)]) {
      
      vector_binary <- MERGE_binary_ct[, c]
      vector_binary <- vector_binary[!is.na(vector_binary)]
      vector_SUM <- sum(vector_binary)
      
      # Check if feature has >= 3 in each group
      if ((vector_SUM >= 3 & vector_SUM <= (length(vector_binary) - 3)) & 
          vector_SUM != length(vector_binary)) {
        binary_features <- c(binary_features, c)
      }
    }
    
    MERGE_binary_ct_temp <- MERGE_binary_ct[
      , colnames(MERGE_binary_ct) %in% c("sample_ID", binary_features)
    ]
    
    ## Filter 2: Remove features identical to covariates
    # Why: Covariates are already in the model; identical features cause multicollinearity
    remove_columns <- c()
    
    for (c3 in colnames(MERGE_binary_ct_temp)[2:ncol(MERGE_binary_ct_temp)]) {
      
      df_temp1 <- MERGE_binary_ct_temp[
        , colnames(MERGE_binary_ct_temp) %in% c("sample_ID", c3)
      ]
      
      for (c4 in colnames(covariates)[2:ncol(covariates)]) {
        
        df_temp2 <- covariates[, colnames(covariates) %in% c("sample_ID", c4)]
        df_temp3 <- merge(df_temp1, df_temp2, by = "sample_ID")
        df_temp3 <- df_temp3[complete.cases(df_temp3), ]
        
        # Check if vectors are identical
        if (identical(df_temp3[[2]], df_temp3[[3]])) {
          remove_columns <- c(remove_columns, c3)
        }
      }
    }
    
    MERGE_binary_ct_temp <- MERGE_binary_ct_temp[
      , !colnames(MERGE_binary_ct_temp) %in% remove_columns
    ]
    
    ## Filter 3: Remove duplicate features (perfect correlation)
    # Uses digest() to compare entire vectors efficiently
    nondups <- MERGE_binary_ct_temp[!duplicated(lapply(MERGE_binary_ct_temp, digest))]
    dups <- MERGE_binary_ct_temp[duplicated(lapply(MERGE_binary_ct_temp, digest))]
    
    # If duplicates exist, identify and record them
    if (dim(dups)[1] > 0) {
      
      duplicated_names <- NULL
      for (i in 1:ncol(nondups)) {
        for (j in 1:ncol(dups)) {
          if (FALSE %in% paste0(nondups[, i] == dups[, j])) {
            NULL
          } else {
            duplicated_names <- rbind(
              duplicated_names, 
              data.frame(names(nondups[i]), names(dups[j]))
            )
          }
        }
      }
      
      colnames(duplicated_names) <- c("Feature", "Feature_duplicated")
      MERGE_binary_ct_temp <- MERGE_binary_ct_temp[
        , !colnames(MERGE_binary_ct_temp) %in% colnames(dups)
      ]
    }
    
    ## Filter 4: Identify highly correlated features using chi-square test
    # Goal: Flag features that may represent redundant information
    # Threshold: FDR-adjusted p-value < 0.5
    
    CHI_SQ <- NULL
    couples <- c()
    
    for (c1 in colnames(MERGE_binary_ct_temp)[2:ncol(MERGE_binary_ct_temp)]) {
      for (c2 in colnames(MERGE_binary_ct_temp)[2:ncol(MERGE_binary_ct_temp)]) {
        
        # Create unique pair identifier
        couple_paste <- ifelse(
          c1 < c2, 
          paste(c1, c2, sep = "-"), 
          paste(c2, c1, sep = "-")
        )
        couple_paste2 <- unlist(strsplit(couple_paste, "-"))
        
        # Test each pair only once
        if (c1 != c2 & !(couple_paste %in% couples)) {
          
          df_temp4 <- MERGE_binary_ct_temp[, c(c1, c2)]
          df_temp4 <- df_temp4[complete.cases(df_temp4), ]
          df_temp4 <- df_temp4 %>% mutate_all(as.numeric)
          
          # Perform chi-square test if both features have variability
          if ((sum(df_temp4[, c1]) < nrow(df_temp4) & sum(df_temp4[, c1]) > 0) | 
              (sum(df_temp4[, c2]) < nrow(df_temp4) & sum(df_temp4[, c2]) > 0)) {
            
            chi_sq_table <- table(df_temp4[, c1], df_temp4[, c2])
            chi_qr_test <- chisq.test(chi_sq_table)
            chi_qr_result <- data.frame(
              chi_qr_test$p.value, 
              couple_paste2[1], 
              couple_paste2[2]
            )
            CHI_SQ <- rbind(CHI_SQ, chi_qr_result)
            couples <- c(couples, paste(c1, c2, sep = "-"), paste(c2, c1, sep = "-"))
            
          } else {
            # No variability - set p-value to 0
            chi_qr_result_0 <- data.frame(0, couple_paste2[1], couple_paste2[2])
            colnames(chi_qr_result_0)[1] <- "chi_qr_test.p.value"
            CHI_SQ <- rbind(CHI_SQ, chi_qr_result_0)
            couples <- c(couples, paste(c1, c2, sep = "-"), paste(c2, c1, sep = "-"))
          }
        }
      }
    }
    
    # Apply FDR correction
    CHI_SQ$p.value.adj <- p.adjust(CHI_SQ$chi_qr_test.p.value, method = "BH")
    CHI_SQ$group <- paste(CHI_SQ$couple_paste2.1., CHI_SQ$couple_paste2.2., sep = "-")
    
    # Identify significantly correlated feature pairs (FDR < 0.5)
    # Note: Using 0.5 threshold to be conservative about flagging correlations
    SIG_group <- CHI_SQ[CHI_SQ$p.value.adj < 0.5, ]$group
  }
  
  ############################################################################
  # 8d. FILTER EXPRESSION FEATURES
  ############################################################################
  
  ## Goal: Remove lowly expressed or invariant genes
  
  if (!is.null(RNASEQ_ct)) {
    
    print("Filtering expression features...")
    
    ## Filter 1: Remove genes not expressed in >= 5% of samples
    # Rationale: Genes expressed in <5% of samples have insufficient data
    no_expression <- as.data.frame(colSums(RNASEQ_ct < 1))
    no_expression$gene <- rownames(no_expression)
    no_expression <- no_expression[
      no_expression$`colSums(RNASEQ_ct < 1)` < nrow(RNASEQ_ct) * 0.95, 
    ]
    
    RNASEQ_ct_temp <- RNASEQ_ct[
      , colnames(RNASEQ_ct) %in% no_expression$gene | colnames(RNASEQ_ct) == "sample_ID"
    ]
    
    ## Filter 2: Require >= 3 samples in different expression categories
    # Categories: Not expressed (<1), Expressed (1-5), Overexpressed (>5)
    # Rationale: Need variability across expression levels for biomarker detection
    categories <- NULL
    
    for (col in colnames(RNASEQ_ct_temp)[2:ncol(RNASEQ_ct_temp)]) {
      
      categories_temp <- RNASEQ_ct_temp[
        , colnames(RNASEQ_ct_temp) %in% c("sample_ID", col)
      ]
      colnames(categories_temp)[2] <- "gene"
      
      # Count samples in each expression category
      not_expressed <- nrow(categories_temp[categories_temp$gene < 1, ])
      expressed <- nrow(categories_temp[categories_temp$gene >= 1 & categories_temp$gene <= 5, ])
      over_expressed <- nrow(categories_temp[categories_temp$gene > 5, ])
      
      # Keep gene if it has >= 3 samples in at least two different categories
      if ((not_expressed >= 3 & expressed >= 3) | 
          (not_expressed >= 3 & over_expressed >= 3) | 
          (expressed >= 3 & over_expressed >= 3)) {
        categories <- c(categories, col)
      }
    }
    
    RNASEQ_ct_temp <- RNASEQ_ct_temp[
      , colnames(RNASEQ_ct_temp) %in% categories | colnames(RNASEQ_ct_temp) == "sample_ID"
    ]
    
    ## Filter 3: Remove genes with low variance using Z-score test
    # Goal: Keep only genes with high variability across samples
    # Method: Calculate Z-score of standard deviation for each gene
    # Threshold: Keep genes with Z-score >= 2 (top ~2.5% most variable)
    SD_values <- NULL
    
    for (col in colnames(RNASEQ_ct_temp)[2:ncol(RNASEQ_ct_temp)]) {
      
      vector <- RNASEQ_ct_temp[, col]
      SD <- sd(vector)
      
      SD_table <- data.frame(col, SD)
      SD_values <- rbind(SD_values, SD_table)
    }
    
    # Calculate Z-scores for standard deviations
    MEAN_all_SD <- mean(SD_values$SD)
    SD_all_SD <- sd(SD_values$SD)
    SD_values$Z.SCORE <- (SD_values$SD - MEAN_all_SD) / SD_all_SD
    
    # Keep high-variance genes OR genes that are response variables
    MERGE_numeric_ct_temp <- RNASEQ_ct_temp[
      , colnames(RNASEQ_ct_temp) %in% SD_values[SD_values$Z.SCORE >= 2, ]$col | 
        colnames(RNASEQ_ct_temp) %in% paste(colnames(LFC_filtered_ct), "exp", sep = "_") | 
        colnames(RNASEQ_ct_temp) == "sample_ID"
    ]
  }
  
  ############################################################################
  # 8e. SCALE ALL FEATURES TO [0,1] FOR COMPARABLE COEFFICIENTS
  ############################################################################
  
  ## Why scale? Allows direct comparison of effect sizes across different feature types
  # All features (binary, expression) are scaled to 0-1 range
  DATA_scale <- merge(MERGE_binary_ct_temp, MERGE_numeric_ct_temp, by = "sample_ID", all = TRUE)
  DATA_scale[2:ncol(DATA_scale)] <- data.frame(
    lapply(DATA_scale[2:ncol(DATA_scale)], function(x) rescale(x, to = c(0, 1)))
  )
  
  ############################################################################
  # 8f. CLASSIFY FEATURE TYPES
  ############################################################################
  
  ## Assign each feature to a type for organized analysis and reporting
  print("Classifying feature types...")
  
  FEATURES_TYPE <- as.data.frame(colnames(DATA_scale)[2:ncol(DATA_scale)])
  colnames(FEATURES_TYPE)[1] <- "FEATURES"
  FEATURES_TYPE$FEATURES <- as.character(FEATURES_TYPE$FEATURES)
  
  # Default: Mutation
  FEATURES_TYPE$type <- "Mut"
  
  # Specific variant (e.g., KRAS G12D)
  if (!is.null(allele_mutations) && 
      nrow(FEATURES_TYPE[FEATURES_TYPE$FEATURES %in% colnames(allele_mutations)[2:ncol(allele_mutations)], ]) > 0) {
    FEATURES_TYPE[FEATURES_TYPE$FEATURES %in% colnames(allele_mutations)[2:ncol(allele_mutations)], ]$type <- "Var"
  }
  
  # Copy number alteration
  if (!is.null(SCNA) && 
      nrow(FEATURES_TYPE[FEATURES_TYPE$FEATURES %in% colnames(SCNA)[2:ncol(SCNA)], ]) > 0) {
    FEATURES_TYPE[FEATURES_TYPE$FEATURES %in% colnames(SCNA)[2:ncol(SCNA)], ]$type <- "CN"
  }
  
  # Loss-of-Function or Gain-of-Function event
  if (!is.null(LoF_GoF) && 
      nrow(FEATURES_TYPE[FEATURES_TYPE$FEATURES %in% colnames(LoF_GoF)[2:ncol(LoF_GoF)], ]) > 0) {
    FEATURES_TYPE[FEATURES_TYPE$FEATURES %in% colnames(LoF_GoF)[2:ncol(LoF_GoF)], ]$type <- "LoF_GoF"
  }
  
  # Mutational signature
  if (!is.null(signatures) && 
      nrow(FEATURES_TYPE[FEATURES_TYPE$FEATURES %in% colnames(signatures)[2:ncol(signatures)], ]) > 0) {
    FEATURES_TYPE[FEATURES_TYPE$FEATURES %in% colnames(signatures)[2:ncol(signatures)], ]$type <- "Sig"
  }
  
  # Clinical variable
  if (!is.null(clinical) && 
      nrow(FEATURES_TYPE[FEATURES_TYPE$FEATURES %in% colnames(clinical)[2:ncol(clinical)], ]) > 0) {
    FEATURES_TYPE[FEATURES_TYPE$FEATURES %in% colnames(clinical)[2:ncol(clinical)], ]$type <- "Clin"
  }
  
  # Expression (RNA-seq)
  if (!is.null(MERGE_numeric_ct_temp) && 
      nrow(FEATURES_TYPE[FEATURES_TYPE$FEATURES %in% colnames(MERGE_numeric_ct_temp)[2:ncol(MERGE_numeric_ct_temp)], ]) > 0) {
    FEATURES_TYPE[FEATURES_TYPE$FEATURES %in% colnames(MERGE_numeric_ct_temp)[2:ncol(MERGE_numeric_ct_temp)], ]$type <- "Expr"
  }
  
  # List of response variables (gene LFCs to test)
  RESPONSES <- colnames(LFC_filtered_ct)[2:ncol(LFC_filtered_ct)]
  
  ############################################################################
  # 8g. LINEAR REGRESSION: MAIN ANALYSIS LOOP
  ############################################################################
  
  ## Nested loop structure:
  ## FOR each gene (response)
  ##   FOR each feature type
  ##     FOR each feature
  ##       - Fit two models: LM0 (without feature) and LM1 (with feature)
  ##       - Perform Likelihood Ratio Test (LRT)
  ##       - Calculate effect sizes and statistics
  
  print("Running linear regression analysis...")
  
  for (response in RESPONSES) {
    
    if (match(response, RESPONSES) %% 100 == 0) {
      print(paste("Processing gene", match(response, RESPONSES), "of", length(RESPONSES)))
    }
    
    # Extract LFC for current gene
    LFC_temp <- LFC_filtered_ct[, colnames(LFC_filtered_ct) %in% c("sample_ID", response)]
    LFC_temp <- LFC_temp[complete.cases(LFC_temp), ]
    
    # Storage for current gene
    DF_TYPE <- NULL
    
    ## Loop through feature types
    for (type in unique(FEATURES_TYPE$type)) {
      
      # Get features of current type
      DATA_scale_temp <- DATA_scale[
        , colnames(DATA_scale) %in% c("sample_ID", FEATURES_TYPE[FEATURES_TYPE$type == type, ]$FEATURES)
      ]
      
      FEATURES <- colnames(DATA_scale_temp)[2:ncol(DATA_scale_temp)]
      
      # Storage for current feature type
      linear_model_results_FEATURES <- NULL
      
      ## Loop through individual features
      for (feature in FEATURES) {
        
        # Extract current feature
        DATA_scale_temp2 <- DATA_scale_temp[
          , colnames(DATA_scale_temp) %in% c("sample_ID", feature)
        ]
        
        # Merge LFC + feature + covariates
        TABLE <- merge(LFC_temp, DATA_scale_temp2, by = "sample_ID")
        colnames(TABLE)[2:3] <- c("response", "feature")
        
        # Add covariates
        if (sum(covariates$msStatus, na.rm = TRUE) == 0 | 
            sum(covariates$msStatus, na.rm = TRUE) == nrow(covariates)) {
          # MSI has no variability - don't include
          TABLE <- merge(TABLE, covariates, by = "sample_ID", all.x = TRUE)
        } else {
          # Include MSI as covariate
          TABLE <- merge(TABLE, covariates, by = "sample_ID", all.x = TRUE)
        }
        
        TABLE$sample_ID <- NULL
        TABLE <- TABLE[!is.na(TABLE$feature), ]
        
        ###################################################################
        # FIT LINEAR MODELS
        ###################################################################
        
        ## Model 0 (null): LFC ~ covariates only
        # Does not include the feature being tested
        LM0 <- lm(response ~ . - feature, data = TABLE, na.action = "na.exclude")
        
        ## Model 1 (full): LFC ~ covariates + feature
        # Includes the feature being tested
        LM1 <- lm(response ~ ., data = TABLE, na.action = "na.exclude")
        
        ## Likelihood Ratio Test: Is Model 1 significantly better than Model 0?
        # H0: Feature does not improve model fit
        # H1: Feature significantly improves model fit
        LRT <- lrtest(LM1, LM0)
        
        ###################################################################
        # EXTRACT MODEL COEFFICIENTS
        ###################################################################
        
        # Extract coefficient for the feature (beta, p-value, R-squared)
        temp_coef <- data.frame(
          response, 
          feature, 
          type, 
          summary(LM1)$coefficients[, 1],  # Coefficient estimate
          summary(LM1)$coefficients[, 4],  # P-value
          summary(LM1)$r.squared           # R-squared
        )
        
        temp_coef$variable <- rownames(temp_coef)
        temp_coef <- temp_coef[temp_coef$variable == "feature", ]
        temp_coef$variable <- NULL
        
        # Handle cases where model failed to fit
        if (dim(temp_coef)[1] == 0) {
          temp_coef <- data.frame(response, feature, type, NA, NA, NA)
        }
        
        # Add LRT p-value
        temp_coef <- data.frame(temp_coef, LRT$`Pr(>Chisq)`[2])
        colnames(temp_coef) <- c(
          'NAME',           # Gene name (response)
          'Feature',        # Feature name
          'Feature_type',   # Feature type (Mut, CN, Expr, etc.)
          'b_coef',         # Beta coefficient
          'P_val_b',        # P-value for coefficient
          'R_squared_LR',   # R-squared
          "P_val_LRT"       # P-value from Likelihood Ratio Test
        )
        
        ###################################################################
        # CALCULATE EFFECT SIZE STATISTICS
        ###################################################################
        
        ## For binary features: Calculate group statistics
        if (feature %in% colnames(MERGE_binary_ct_temp)) {
          
          # Group statistics (mean, median, min, max, SD) for each group
          stats <- TABLE %>% 
            dplyr::group_by(feature) %>% 
            dplyr::summarise(
              MEAN = mean(response),
              MEDIAN = median(response),
              MIN = min(response),
              MAX = max(response),
              SD = sd(response)
            )
          
          # Group 0 (feature absent)
          stats_group1 <- stats[1, ]
          colnames(stats_group1) <- paste(colnames(stats_group1), "group1", sep = '_')
          stats_group1$feature_group1 <- NULL
          
          # Group 1 (feature present)
          stats_group2 <- stats[2, ]
          colnames(stats_group2) <- paste(colnames(stats_group2), "group2", sep = '_')
          stats_group2$feature_group2 <- NULL
          
          # Combine and calculate differences
          STATS <- cbind(stats_group1, stats_group2)
          STATS$MEAN_difference <- STATS$MEAN_group2 - STATS$MEAN_group1
          STATS$MEDIAN_difference <- STATS$MEDIAN_group2 - STATS$MEDIAN_group1
          
          # Cohen's coefficient (effect size)
          # Measures strength of association
          CORR <- cor(TABLE[, 1], TABLE[, 2])
          STATS$cohens_coef <- (CORR^2) / (1 - (CORR)^2)
          
        } else {
          ## For continuous features: Use quantile comparison
          # Compare top 25% vs bottom 25% of expression
          
          TABLE2 <- TABLE
          TABLE2$upper_quant <- TABLE2$feature >= quantile(TABLE2$feature, probs = 0.75)
          TABLE2$lower_quant <- TABLE2$feature <= quantile(TABLE2$feature, probs = 0.25)
          TABLE2$quants <- "X"
          TABLE2[TABLE2$upper_quant == "TRUE", ]$quants <- "upper"
          TABLE2[TABLE2$lower_quant == "TRUE", ]$quants <- "lower"
          
          # Calculate statistics for upper and lower quartiles
          stats <- TABLE2 %>% 
            dplyr::group_by(quants) %>% 
            dplyr::summarise(
              MEAN = mean(response),
              MEDIAN = median(response),
              MIN = min(response),
              MAX = max(response),
              SD = sd(response)
            )
          
          stats <- stats[!stats$quants == "X", ]
          
          # Lower quartile (group 1)
          stats_group1 <- stats[1, ]
          colnames(stats_group1) <- paste(colnames(stats_group1), "group1", sep = '_')
          stats_group1$quants_group1 <- NULL
          
          # Upper quartile (group 2)
          stats_group2 <- stats[2, ]
          colnames(stats_group2) <- paste(colnames(stats_group2), "group2", sep = '_')
          stats_group2$quants_group2 <- NULL
          
          # Combine and calculate differences
          STATS <- cbind(stats_group1, stats_group2)
          STATS$MEAN_difference <- STATS$MEAN_group2 - STATS$MEAN_group1
          STATS$MEDIAN_difference <- STATS$MEDIAN_group2 - STATS$MEDIAN_group1
          
          # Cohen's coefficient
          CORR <- cor(TABLE[, 1], TABLE[, 2])
          STATS$cohens_coef <- (CORR^2) / (1 - (CORR)^2)
        }
        
        # Combine model results with effect size statistics
        DF <- data.frame(temp_coef, STATS)
        linear_model_results_FEATURES <- rbind(linear_model_results_FEATURES, DF)
      }
      
      ## Apply FDR correction within feature type
      # Why: Different feature types may have different multiple testing burdens
      linear_model_results_FEATURES$P_val_LRT_BH <- p.adjust(
        linear_model_results_FEATURES$P_val_LRT, 
        method = "BH"
      )
      
      DF_TYPE <- rbind(DF_TYPE, linear_model_results_FEATURES)
    }
    
    ###################################################################
    # IDENTIFY HIGHLY CORRELATED FEATURES
    ###################################################################
    
    ## For binary features, flag those that showed significant correlation
    # in chi-square test (from earlier filtering step)
    print("Flagging correlated features...")
    
    DF_TYPE$similar_variable <- NA
    
    # Get binary features sorted by p-value (worst to best)
    binary_variables <- DF_TYPE[DF_TYPE$Feature %in% colnames(MERGE_binary_ct_temp), ]
    binary_variables <- binary_variables[order(binary_variables$P_val_LRT, decreasing = TRUE), ]
    
    for (feature1 in binary_variables$Feature) {
      
      # Remove current feature from list
      binary_variables <- binary_variables[!binary_variables$Feature == feature1, ]
      
      # Check correlation with remaining features
      for (feature2 in binary_variables$Feature) {
        
        # Create pair identifier (alphabetically sorted)
        B <- ifelse(
          feature1 < feature2, 
          paste(feature1, feature2, sep = "-"), 
          paste(feature2, feature1, sep = "-")
        )
        
        # If pair was significantly correlated in chi-square test
        if (B %in% SIG_group) {
          
          # Flag as similar variable
          DF_TYPE[DF_TYPE$Feature == feature1, ]$similar_variable <- paste(
            DF_TYPE[DF_TYPE$Feature == feature1, ]$similar_variable, 
            feature2, 
            sep = "/"
          )
          
          print(paste(response, B, sep = ","))
        }
      }
    }
    
    # Clean up NA/ prefix
    DF_TYPE$similar_variable <- gsub("^NA/", "", DF_TYPE$similar_variable)
    DF_TYPE$Feature <- as.character(DF_TYPE$Feature)
    
    # Accumulate results for all genes
    linear_model_results_RESPONSES <- rbind(linear_model_results_RESPONSES, DF_TYPE)
  }
  
  ## Add information about duplicate features (if any were found)
  if (dim(dups)[1] > 0) {
    linear_model_results_RESPONSES <- merge(
      linear_model_results_RESPONSES, 
      duplicated_names, 
      by = "Feature", 
      all.x = TRUE
    )
  }
  
  # Add cancer type label
  linear_model_results_RESPONSES$cancer_type <- ct
}



################################################################################
# 9. BIOMARKER CLASSIFICATION
################################################################################

## Classify associations as biomarkers based on statistical and biological criteria

## Criteria for biomarker classification:
# 1. FDR < 0.05 (statistically significant)
# 2. Mean LFC in one group < -0.5 (biologically meaningful dependency)
# 3. Groups show opposite directions (one depleted, one not)

linear_model_results_RESPONSES$biomarker_class <- "No"

# Biomarker: Group 2 shows strong dependency (< -0.5) AND is more dependent than Group 1
linear_model_results_RESPONSES[
  linear_model_results_RESPONSES$P_val_LRT_BH < 0.05 & 
  linear_model_results_RESPONSES$MEAN_group2 < -0.5 & 
  linear_model_results_RESPONSES$MEAN_group1 > linear_model_results_RESPONSES$MEAN_group2, 
]$biomarker_class <- "Yes"

# Biomarker: Group 1 shows strong dependency (< -0.5) AND is more dependent than Group 2
linear_model_results_RESPONSES[
  linear_model_results_RESPONSES$P_val_LRT_BH < 0.05 & 
  linear_model_results_RESPONSES$MEAN_group1 < -0.5 & 
  linear_model_results_RESPONSES$MEAN_group2 > linear_model_results_RESPONSES$MEAN_group1, 
]$biomarker_class <- "Yes"

# Create unique association identifier
linear_model_results_RESPONSES$association <- paste(
  linear_model_results_RESPONSES$NAME, 
  linear_model_results_RESPONSES$Feature, 
  sep = "-"
)



################################################################################
# 10. REFINE BIOMARKER CLASSIFICATION BY EFFECT SIZE
################################################################################

## Extract significant associations (those classified as biomarkers)
significant <- linear_model_results_RESPONSES[
  linear_model_results_RESPONSES$biomarker_class == "Yes", 
]

## Calculate effect size quantiles
# These will be used to classify biomarkers into tiers (A, B, C)
significant <- significant %>% 
  mutate(cohens_coef_quantile_90 = abs(cohens_coef) > quantile(abs(cohens_coef), probs = 0.90, na.rm = TRUE))

significant <- significant %>% 
  mutate(cohens_coef_quantile_95 = abs(cohens_coef) > quantile(abs(cohens_coef), probs = 0.95, na.rm = TRUE))

## Determine direction of association
# "Increased Dep." = Feature associated with stronger dependency (more negative LFC)
# "Decreased Dep." = Feature associated with weaker dependency (less negative LFC)
significant$ASSOCIATION_EFFECT <- "Decreased Dep."
significant[significant$b_coef < 0, ]$ASSOCIATION_EFFECT <- "Increased Dep."

## Filter for high-effect associations
# Requirement: Cohen's coefficient in top 10% (quantile >= 0.90)
# Rationale: Focus on strongest, most robust associations
significant <- significant[significant$cohens_coef_quantile_90 == "TRUE", ]



################################################################################
# 11. CLASSIFY BIOMARKERS INTO PRIORITY TIERS
################################################################################

## Tier C (Default): Significant but modest effect
# Already assigned "Yes" from previous step

## Tier B: Moderate effect size
# Criteria: |Mean difference| > 0.5
# Interpretation: Half-unit difference in scaled LFC between groups
significant[abs(significant$MEAN_difference) > 0.5, ]$biomarker_class <- "B"

## Tier A: Strong effect size
# Criteria: Cohen's coefficient in top 5% AND |Mean difference| > 0.75
# Interpretation: Strongest, most clinically relevant associations
significant[
  significant$cohens_coef_quantile_95 == "TRUE" & 
    abs(significant$MEAN_difference) > 0.75, 
]$biomarker_class <- "A"

## Remaining associations stay as Tier C
significant[significant$biomarker_class == "Yes", ]$biomarker_class <- "C"

## Summary of classification:
# Class A: Top-tier biomarkers (largest effects, highest confidence)
# Class B: Medium-tier biomarkers (moderate effects)
# Class C: Lower-tier biomarkers (significant but smaller effects)

# NOTE: Save significant as supplementary_table_8.1



################################################################################
# 12. PRIORITY SCORING SYSTEM
################################################################################

## Goal: Select single "best" biomarker for each gene
# Challenge: Genes may have multiple significant biomarkers
# Solution: Prioritize by class, then feature type, then p-value

## Load cancer driver gene annotations
# Used to flag biomarkers for known cancer genes
cancer_driver_genes <- as.data.frame(
  fread('/path_to_folder/cancer_driver_genes.tsv')
)

## Test association with MSI status
# Why: MSI is a major confounding factor in colorectal cancer
# Method: Fisher's exact test for each binary feature
# Threshold: FDR < 0.25 (exploratory, will flag rather than exclude)

# Get binary features tested in this cancer type
binary_features_MSI <- MERGE_binary[
  , colnames(MERGE_binary) %in% 
    unique(linear_model_results_RESPONSES[
      linear_model_results_RESPONSES$cancer_type == "Gastrointestinal", 
    ]$feature) | colnames(MERGE_binary) == "sample_ID"
]

binary_features_MSI <- binary_features_MSI[
  binary_features_MSI$sample_ID %in% LFC_filtered$sample_ID, 
]

binary_features_MSI <- merge(MSI, binary_features_MSI, by = "sample_ID", all.y = TRUE)

# Test each feature for association with MSI
MSI_features <- NULL

for (f in colnames(binary_features_MSI)[3:ncol(binary_features_MSI)]) {
  
  df_temp <- binary_features_MSI[
    , colnames(binary_features_MSI) %in% c("sample_ID", "msStatus", f)
  ]
  
  # Chi-square test
  CHI.TEST <- chisq.test(df_temp[[2]], df_temp[[3]])
  
  # Fisher's exact test (more appropriate for small counts)
  FISHER.TEST <- fisher.test(df_temp[[2]], df_temp[[3]])
  
  MSI_features <- rbind(MSI_features, data.frame(f, FISHER.TEST$p.value))
}

# Apply FDR correction
MSI_features$p.adj <- p.adjust(MSI_features$FISHER.TEST.p.value, method = "BH")



################################################################################
# 13. PRIORITY SCORING - GASTROINTESTINAL CANCERS
################################################################################

## Subset to gastrointestinal cancers
significant_GAST <- significant[significant$ct == "gastrointestinal", ]

## Priority scoring algorithm for each gene
# Decision tree:
# 1. Select best biomarker class (A > B > C)
# 2. If multiple biomarkers in same class:
#    a. Prefer genomic features (Mut, CN, Clin, Var) over expression
#    b. If still tied, select by lowest p-value
# 3. Flag if associated with MSI (FDR < 0.25)

significant_GAST2 <- NULL

for (x in unique(significant_GAST$NAME)) {
  
  print(x)
  
  # Get all biomarkers for current gene
  df_temp <- significant_GAST[significant_GAST$NAME == x, ]
  
  # Sort by biomarker class (A < B < C alphabetically, so A is best)
  df_temp <- df_temp[order(df_temp$biomarker_class), ]
  
  # Select best class
  MAX_biomarker <- min(df_temp$biomarker_class)
  df_temp2 <- df_temp[df_temp$biomarker_class == MAX_biomarker, ]
  
  # Count number of feature types in best class
  TYPES <- length(unique(df_temp2$Feature_type))
  
  df_temp3 <- NULL
  
  if (TYPES > 1) {
    ## Multiple feature types: Prioritize genomic over expression
    # Rationale: Genomic features are more stable and actionable
    
    TEMP2 <- df_temp2[df_temp2$Feature_type %in% c("Mut", "CN", "Clin", "Var"), ]
    
    if (nrow(TEMP2) > 0) {
      # Keep genomic features
      df_temp3 <- rbind(df_temp3, TEMP2)
    } else {
      # No genomic features, keep all
      df_temp3 <- rbind(df_temp3, df_temp2)
    }
    
  } else {
    ## Single feature type: Randomly select one if multiple features
    # Note: In practice, could select by lowest p-value instead
    df_temp2 <- df_temp2 %>% sample_n(1)
    df_temp3 <- rbind(df_temp3, df_temp2)
  }
  
  # Record number of biomarkers selected
  N_MAX <- nrow(df_temp3)
  
  # Identify feature with lowest p-value (most significant)
  MAX_pvalue <- min(df_temp3$P_val_LRT_BH)
  MAX_feature <- df_temp3[df_temp3$P_val_LRT_BH == MAX_pvalue, ]$Feature
  MAX_feature <- paste(MAX_feature, collapse = "-")
  MAX_feature_type <- df_temp3[df_temp3$P_val_LRT_BH == MAX_pvalue, ]$Feature_type
  MAX_feature_type <- paste(MAX_feature_type, collapse = "-")
  
  # Combine selected biomarkers with metadata
  significant_GAST2 <- rbind(
    significant_GAST2, 
    data.frame(df_temp3, MAX_feature, MAX_feature_type, N_MAX)
  )
}

## Flag MSI-associated features
# Threshold: FDR < 0.25 (exploratory threshold)
significant_GAST2$MSI_feature <- "No"
significant_GAST2[
  significant_GAST2$Feature %in% MSI_features[MSI_features$p.adj < 0.25, ]$f, 
]$MSI_feature <- "Yes"

## Merge with differential dependency results
# Adds effect size information from differential dependency analysis
priority_score_GAST <- merge(
  differential_dep_gastrointestinal, 
  significant_GAST2, 
  by = "NAME", 
  all.x = TRUE
)

## Format priority score table
# Select key columns for publication
priority_score_GAST2 <- priority_score_GAST[, c(
  "NAME",                          # Gene name
  "MAX_feature",                   # Best biomarker feature
  "MAX_feature_type",              # Feature type
  "Feature",                       # All selected features
  "Feature_type",                  # All feature types
  "delta_MEAN",                    # Effect size from diff. dep.
  "MEAN_LFC_depleted_group",       # Mean LFC in depleted group
  "MEAN_LFC_not_depleted_group",   # Mean LFC in non-depleted group
  "N_depleted",                    # N samples depleted
  "N_not_depleted",                # N samples not depleted
  "biomarker_class",               # Class (A, B, or C)
  "association",                   # Gene-feature pair
  "MSI_feature"                    # MSI association flag
)]

## Collapse multiple entries per gene (if gene has multiple selected biomarkers)
# Use "/" as delimiter between multiple values
priority_score_GAST3 <- priority_score_GAST2 %>% 
  group_by(NAME) %>% 
  summarise_all(funs(paste(na.omit(.), collapse = "/")))

# Remove duplicates within collapsed fields
priority_score_GAST3$biomarker_class <- sapply(
  priority_score_GAST3$biomarker_class, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)

priority_score_GAST3$delta_MEAN <- sapply(
  priority_score_GAST3$delta_MEAN, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)

priority_score_GAST3$MEAN_LFC_depleted_group <- sapply(
  priority_score_GAST3$MEAN_LFC_depleted_group, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)

priority_score_GAST3$MEAN_LFC_not_depleted_group <- sapply(
  priority_score_GAST3$MEAN_LFC_not_depleted_group, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)

priority_score_GAST3$N_depleted <- sapply(
  priority_score_GAST3$N_depleted, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)

priority_score_GAST3$N_not_depleted <- sapply(
  priority_score_GAST3$N_not_depleted, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)

priority_score_GAST3$MAX_feature <- sapply(
  priority_score_GAST3$MAX_feature, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)

priority_score_GAST3$MAX_feature_type <- sapply(
  priority_score_GAST3$MAX_feature_type, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)

# Convert back to numeric where appropriate
priority_score_GAST3$delta_MEAN <- as.numeric(priority_score_GAST3$delta_MEAN)
priority_score_GAST3$MEAN_LFC_depleted_group <- as.numeric(priority_score_GAST3$MEAN_LFC_depleted_group)
priority_score_GAST3$MEAN_LFC_not_depleted_group <- as.numeric(priority_score_GAST3$MEAN_LFC_not_depleted_group)
priority_score_GAST3$N_depleted <- as.numeric(priority_score_GAST3$N_depleted)
priority_score_GAST3$N_not_depleted <- as.numeric(priority_score_GAST3$N_not_depleted)

## Flag cancer driver genes
priority_score_GAST3$cancer_driver <- "No"
priority_score_GAST3[
  priority_score_GAST3$NAME %in% cancer_driver_genes$gene, 
]$cancer_driver <- "Yes"

# Handle genes with no biomarkers (all values empty string)
priority_score_GAST3[priority_score_GAST3$biomarker_class == "", ]$biomarker_class <- "No"

# NOTE: Save priority_score_GAST3 as supplementary_table_8.2



################################################################################
# 14. PRIORITY SCORING - COLORECTAL CANCERS
################################################################################

## Repeat priority scoring for colorectal-specific analysis
# Same algorithm as gastrointestinal, but cancer type-specific

significant_COLO <- significant[significant$ct == "Colorectal", ]
significant_COLO$cancer_subtype_target <- paste(significant_COLO$cohort, significant_COLO$NAME, sep = "-")

significant_COLO2 <- NULL

for (x in unique(significant_COLO$NAME)) {
  
  print(x)
  
  df_temp <- significant_COLO[significant_COLO$NAME == x, ]
  df_temp <- df_temp[order(df_temp$biomarker_class), ]
  
  MAX_biomarker <- min(df_temp$biomarker_class)
  df_temp2 <- df_temp[df_temp$biomarker_class == MAX_biomarker, ]
  
  TYPES <- length(unique(df_temp2$Feature_type))
  
  df_temp3 <- NULL
  
  if (TYPES > 1) {
    TEMP2 <- df_temp2[df_temp2$Feature_type %in% c("Mut", "CN", "Clin", "Var"), ]
    
    if (nrow(TEMP2) > 0) {
      df_temp3 <- rbind(df_temp3, TEMP2)
    } else {
      df_temp3 <- rbind(df_temp3, df_temp2)
    }
  } else {
    df_temp2 <- df_temp2 %>% sample_n(1)
    df_temp3 <- rbind(df_temp3, df_temp2)
  }
  
  N_MAX <- nrow(df_temp3)
  MAX_pvalue <- min(df_temp3$P_val_LRT_BH)
  
  MAX_feature <- df_temp3[df_temp3$P_val_LRT_BH == MAX_pvalue, ]$Feature
  MAX_feature <- paste(MAX_feature, collapse = "-")
  MAX_feature_type <- df_temp3[df_temp3$P_val_LRT_BH == MAX_pvalue, ]$Feature_type
  MAX_feature_type <- paste(MAX_feature_type, collapse = "-")
  
  significant_COLO2 <- rbind(
    significant_COLO2, 
    data.frame(df_temp3, MAX_feature, MAX_feature_type, N_MAX)
  )
}

## Flag MSI-associated features
significant_COLO2$MSI_feature <- "No"
significant_COLO2[
  significant_COLO2$Feature %in% MSI_features[MSI_features$p.adj < 0.25, ]$f, 
]$MSI_feature <- "Yes"

## Merge and format
priority_score_COLO <- merge(
  differential_dep_COADREAD, 
  significant_COLO2, 
  by = "NAME", 
  all.x = TRUE
)

priority_score_COLO2 <- priority_score_COLO[, c(
  "NAME", "MAX_feature", "MAX_feature_type", "Feature", "Feature_type", 
  "delta_MEAN", "MEAN_LFC_depleted_group", "MEAN_LFC_not_depleted_group", 
  "N_depleted", "N_not_depleted", "biomarker_class", "association", "MSI_feature"
)]

# Collapse multiple entries
priority_score_COLO3 <- priority_score_COLO2 %>% 
  group_by(NAME) %>% 
  summarise_all(funs(paste(na.omit(.), collapse = "/")))

# Remove duplicates within collapsed fields
priority_score_COLO3$biomarker_class <- sapply(
  priority_score_COLO3$biomarker_class, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)
priority_score_COLO3$delta_MEAN <- sapply(
  priority_score_COLO3$delta_MEAN, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)
priority_score_COLO3$MEAN_LFC_depleted_group <- sapply(
  priority_score_COLO3$MEAN_LFC_depleted_group, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)
priority_score_COLO3$MEAN_LFC_not_depleted_group <- sapply(
  priority_score_COLO3$MEAN_LFC_not_depleted_group, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)
priority_score_COLO3$N_depleted <- sapply(
  priority_score_COLO3$N_depleted, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)
priority_score_COLO3$N_not_depleted <- sapply(
  priority_score_COLO3$N_not_depleted, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)
priority_score_COLO3$MAX_feature <- sapply(
  priority_score_COLO3$MAX_feature, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)
priority_score_COLO3$MAX_feature_type <- sapply(
  priority_score_COLO3$MAX_feature_type, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)

# Convert to numeric
priority_score_COLO3$delta_MEAN <- as.numeric(priority_score_COLO3$delta_MEAN)
priority_score_COLO3$MEAN_LFC_depleted_group <- as.numeric(priority_score_COLO3$MEAN_LFC_depleted_group)
priority_score_COLO3$MEAN_LFC_not_depleted_group <- as.numeric(priority_score_COLO3$MEAN_LFC_not_depleted_group)
priority_score_COLO3$N_depleted <- as.numeric(priority_score_COLO3$N_depleted)
priority_score_COLO3$N_not_depleted <- as.numeric(priority_score_COLO3$N_not_depleted)

## Flag cancer driver genes
priority_score_COLO3$cancer_driver <- "No"
priority_score_COLO3[
  priority_score_COLO3$NAME %in% cancer_driver_genes$gene, 
]$cancer_driver <- "Yes"

priority_score_COLO3[priority_score_COLO3$biomarker_class == "", ]$biomarker_class <- "No"

# NOTE: Save priority_score_COLO3 as supplementary_table_8.3



################################################################################
# 15. PRIORITY SCORING - OESOPHAGEAL CANCERS
################################################################################

## Repeat priority scoring for esophageal-specific analysis

significant_OESO <- significant[significant$ct == "Oesophageal", ]
significant_OESO$cancer_subtype_target <- paste(significant_OESO$cohort, significant_OESO$NAME, sep = "-")

significant_OESO2 <- NULL

for (x in unique(significant_OESO$NAME)) {
  
  print(x)
  
  df_temp <- significant_OESO[significant_OESO$NAME == x, ]
  df_temp <- df_temp[order(df_temp$biomarker_class), ]
  
  MAX_biomarker <- min(df_temp$biomarker_class)
  df_temp2 <- df_temp[df_temp$biomarker_class == MAX_biomarker, ]
  
  TYPES <- length(unique(df_temp2$Feature_type))
  
  df_temp3 <- NULL
  
  if (TYPES > 1) {
    TEMP2 <- df_temp2[df_temp2$Feature_type %in% c("Mut", "CN", "Clin", "Var"), ]
    
    if (nrow(TEMP2) > 0) {
      df_temp3 <- rbind(df_temp3, TEMP2)
    } else {
      df_temp3 <- rbind(df_temp3, df_temp2)
    }
  } else {
    df_temp2 <- df_temp2 %>% sample_n(1)
    df_temp3 <- rbind(df_temp3, df_temp2)
  }
  
  N_MAX <- nrow(df_temp3)
  MAX_pvalue <- min(df_temp3$P_val_LRT_BH)
  
  MAX_feature <- df_temp3[df_temp3$P_val_LRT_BH == MAX_pvalue, ]$Feature
  MAX_feature <- paste(MAX_feature, collapse = "-")
  MAX_feature_type <- df_temp3[df_temp3$P_val_LRT_BH == MAX_pvalue, ]$Feature_type
  MAX_feature_type <- paste(MAX_feature_type, collapse = "-")
  
  significant_OESO2 <- rbind(
    significant_OESO2, 
    data.frame(df_temp3, MAX_feature, MAX_feature_type, N_MAX)
  )
}

## Merge and format
priority_score_OESO <- merge(
  differential_dep_ESCA, 
  significant_OESO2, 
  by = "NAME", 
  all.x = TRUE
)

priority_score_OESO2 <- priority_score_OESO[, c(
  "NAME", "MAX_feature", "MAX_feature_type", "Feature", "Feature_type", 
  "delta_MEAN", "MEAN_LFC_depleted_group", "MEAN_LFC_not_depleted_group", 
  "N_depleted", "N_not_depleted", "biomarker_class", "association"
)]

# Collapse multiple entries
priority_score_OESO3 <- priority_score_OESO2 %>% 
  group_by(NAME) %>% 
  summarise_all(funs(paste(na.omit(.), collapse = "/")))

# Remove duplicates within collapsed fields
priority_score_OESO3$biomarker_class <- sapply(
  priority_score_OESO3$biomarker_class, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)
priority_score_OESO3$delta_MEAN <- sapply(
  priority_score_OESO3$delta_MEAN, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)
priority_score_OESO3$MEAN_LFC_depleted_group <- sapply(
  priority_score_OESO3$MEAN_LFC_depleted_group, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)
priority_score_OESO3$MEAN_LFC_not_depleted_group <- sapply(
  priority_score_OESO3$MEAN_LFC_not_depleted_group, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)
priority_score_OESO3$N_depleted <- sapply(
  priority_score_OESO3$N_depleted, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)
priority_score_OESO3$N_not_depleted <- sapply(
  priority_score_OESO3$N_not_depleted, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)
priority_score_OESO3$MAX_feature <- sapply(
  priority_score_OESO3$MAX_feature, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)
priority_score_OESO3$MAX_feature_type <- sapply(
  priority_score_OESO3$MAX_feature_type, 
  function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
)

# Convert to numeric
priority_score_OESO3$delta_MEAN <- as.numeric(priority_score_OESO3$delta_MEAN)
priority_score_OESO3$MEAN_LFC_depleted_group <- as.numeric(priority_score_OESO3$MEAN_LFC_depleted_group)
priority_score_OESO3$MEAN_LFC_not_depleted_group <- as.numeric(priority_score_OESO3$MEAN_LFC_not_depleted_group)
priority_score_OESO3$N_depleted <- as.numeric(priority_score_OESO3$N_depleted)
priority_score_OESO3$N_not_depleted <- as.numeric(priority_score_OESO3$N_not_depleted)

## Flag cancer driver genes
priority_score_OESO3$cancer_driver <- "No"
priority_score_OESO3[
  priority_score_OESO3$NAME %in% cancer_driver_genes$gene, 
]$cancer_driver <- "Yes"

# NOTE: Save priority_score_OESO3 as supplementary_table_8.4



################################################################################
# SUMMARY OF OUTPUT
################################################################################

# 1. linear_model_results_RESPONSES (all_biomarker-analysis_results.zip)
#    - Complete results for all gene-feature associations tested
#    - Includes non-significant associations
#    - Used for exploratory analysis and QC

# 2. significant (supplementary_table_8.1)
#    - All significant biomarker associations (Classes A, B, C)
#    - Filtered for FDR < 0.05 and biological effect
#    - Multiple biomarkers per gene possible

# 3. priority_score_GAST3 (supplementary_table_8.2)
#    - Priority-scored biomarkers for gastrointestinal cancers
#    - One "best" biomarker per gene
#    - Includes MSI association flags

# 4. priority_score_COLO3 (supplementary_table_8.3)
#    - Priority-scored biomarkers for colorectal cancers
#    - Includes MSI association flags

# 5. priority_score_OESO3 (supplementary_table_8.4)
#    - Priority-scored biomarkers for esophageal cancers
