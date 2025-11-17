################################################################################
# Script: CRISPR_differential_dependency_analysis.R
# Purpose: Identify genes with differential dependencies across organoid models
# Description: Performs differential dependency analysis to identify genes that
#              show significant differences in essentiality between subsets of
#              organoid models, stratified by cancer type
################################################################################

# Load required libraries
library(data.table)
library(tidyverse)

################################################################################
# 1. LOAD INPUT DATA
################################################################################

## Load binary dependency matrix (gene depletion calls from BAGEL2)
BINARY <- as.data.frame(
  fread('/path_to_folder/supplementary_table_5.csv')
)

## Load Log2 Fold Change (LFC) matrix
LFC <- as.data.frame(
  fread('/path_to_folder/supplementary_table_6.csv')
)

## Merge binary calls with LFC values
MERGE_BINARY_LFC <- merge(BINARY, LFC)

## Load cancer type annotations
cancer_type <- as.data.frame(
  fread('/path_to_folder/supplementary_table_2.2.csv')
)
cancer_type <- cancer_type[, c(1, 3)]

## Load gene classification lists

# Essential genes (positive controls)
essential <- as.data.frame(
  fread('/path_to_folder/essential_genes_organoids_downsampling.csv')
)

# Non-essential genes (negative controls)
non_essential <- as.data.frame(
  fread('/path_to_folder/non-essential_genes_organoids_downsampling.csv')
)

# Pan-cancer core fitness genes from cell lines
load('/path_to_folder/10_PANCANCER_coreFitness_genes.RData')

# Pan-cancer core fitness genes from organoids
ORG_PanCancerCoreFitnessGenes <- as.data.frame(
  fread('/path_to_folder/supplementary_table_4.2.csv')
)

## Load RNA-seq data to filter lowly expressed genes
RNAseq <- as.data.frame(
  fread('/path_to_folder/data_RNAseq_ALL-organoids.csv')
)

# Remove expression suffix from column names
colnames(RNAseq) <- gsub("_exp", "", colnames(RNAseq))



################################################################################
# 2. CALCULATE MEAN EXPRESSION LEVELS
################################################################################

## Calculate mean expression for gastrointestinal cancers (all)
RNAseq2 <- RNAseq[RNAseq$sample_ID %in% cancer_type[cancer_type$primary_tumour_type %in% c("Colorectal", "Oesophageal", "Gastric", "Pancreatic"), ]$sample_ID,]
RNAseq2 <- gather(RNAseq2, "gene", "log2TPM", 2:ncol(RNAseq2))

RNAseq_mean <- RNAseq2 %>% 
  group_by(gene) %>% 
  dplyr::summarise(mean_log2TPM = mean(log2TPM))

# Identify genes with low/absent expression (mean log2TPM < 0.1)
RNAseq_not_expressed <- RNAseq_mean[RNAseq_mean$mean_log2TPM < 0.1, ]

## Calculate mean expression for colorectal cancers
RNAseq_COLO <- RNAseq[RNAseq$sample_ID %in% cancer_type[cancer_type$primary_tumour_type == "Colorectal", ]$sample_ID,]
RNAseq_COLO <- gather(RNAseq_COLO, "gene", "log2TPM", 2:ncol(RNAseq_COLO))

RNAseq_mean_COLO <- RNAseq_COLO %>% 
  group_by(gene) %>% 
  dplyr::summarise(mean_log2TPM = mean(log2TPM))

RNAseq_not_expressed_COLO <- RNAseq_mean_COLO[RNAseq_mean_COLO$mean_log2TPM < 0.1, ]

## Calculate mean expression for oesophageal cancers
RNAseq_OESO <- RNAseq[RNAseq$sample_ID %in% cancer_type[cancer_type$primary_tumour_type == "Oesophageal", ]$sample_ID,]
RNAseq_OESO <- gather(RNAseq_OESO, "gene", "log2TPM", 2:ncol(RNAseq_OESO))

RNAseq_mean_OESO <- RNAseq_OESO %>% 
  group_by(gene) %>% 
  dplyr::summarise(mean_log2TPM = mean(log2TPM))

RNAseq_not_expressed_OESO <- RNAseq_mean_OESO[RNAseq_mean_OESO$mean_log2TPM < 0.1, ]



################################################################################
# 3. FILTER GENES FOR ANALYSIS
################################################################################

## Exclude core fitness and control genes
# These genes are constitutively essential and not of interest for biomarker discovery
genes_to_remove <- unique(c(
  ORG_PanCancerCoreFitnessGenes$gene, 
  essential$GENE, 
  non_essential$GENE, 
  PanCancerCoreFitnessGenes
))

MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC[!MERGE_BINARY_LFC$gene %in% genes_to_remove,]

################################################################################
# 4. REMOVE DUPLICATE SAMPLES
################################################################################

## Remove technical replicates and low-quality samples
# These samples were either screened with multiple libraries or failed QC

MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[!(MERGE_BINARY_LFC2$sample_ID == "HCM-SANG-0267-D12" & MERGE_BINARY_LFC2$library == "minLib"), ]
MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[!(MERGE_BINARY_LFC2$sample_ID == "HCM-SANG-0268-C18" & MERGE_BINARY_LFC2$library == "v1_1"), ]
MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[!(MERGE_BINARY_LFC2$sample_ID == "HCM-SANG-0270-C20" & MERGE_BINARY_LFC2$library == "minLib"), ]
MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[!(MERGE_BINARY_LFC2$sample_ID == "HCM-SANG-0272-C20" & MERGE_BINARY_LFC2$library == "v1_1"), ]
MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[!(MERGE_BINARY_LFC2$sample_ID == "HCM-SANG-0273-C18" & MERGE_BINARY_LFC2$library == "minLib"), ]
MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[!(MERGE_BINARY_LFC2$sample_ID == "HCM-SANG-0276-C18" & MERGE_BINARY_LFC2$library == "minLib"), ]
MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[!(MERGE_BINARY_LFC2$sample_ID == "HCM-SANG-0282-C18" & MERGE_BINARY_LFC2$library == "minLib"), ]
MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[!(MERGE_BINARY_LFC2$sample_ID == "HCM-SANG-0529-C18" & MERGE_BINARY_LFC2$library == "v1_1"), ]
MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[!(MERGE_BINARY_LFC2$sample_ID == "HCM-SANG-0284-C18" & MERGE_BINARY_LFC2$library == "minLib"), ]
MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[!(MERGE_BINARY_LFC2$sample_ID == "HCM-SANG-0288-C18" & MERGE_BINARY_LFC2$library == "v1_1"), ]
MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[!(MERGE_BINARY_LFC2$sample_ID == "WTSI-COLO_278" & MERGE_BINARY_LFC2$library == "minLib"), ]
MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[!(MERGE_BINARY_LFC2$sample_ID == "WTSI-COLO_376" & MERGE_BINARY_LFC2$library == "minLib"), ]
MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[!(MERGE_BINARY_LFC2$sample_ID == "HCM-SANG-0291-C15" & MERGE_BINARY_LFC2$library == "v1_1"), ]
MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[!(MERGE_BINARY_LFC2$sample_ID == "HCM-SANG-0299-C15-A" & MERGE_BINARY_LFC2$library == "minLib"), ]
MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[!(MERGE_BINARY_LFC2$sample_ID == "HCM-SANG-0299-C15-B" & MERGE_BINARY_LFC2$library == "minLib"), ]
MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[!(MERGE_BINARY_LFC2$sample_ID == "HCM-SANG-0306-C15" & MERGE_BINARY_LFC2$library == "minLib"), ]

## Add cancer type annotations
MERGE_BINARY_LFC2 <- merge(MERGE_BINARY_LFC2, cancer_type, by = "sample_ID")



################################################################################
# 5. STRATIFY BY CANCER TYPE
################################################################################

## Focus on gastrointestinal cancers (sufficient sample sizes)
# Remove ovarian models
MERGE_BINARY_LFC2 <- MERGE_BINARY_LFC2[MERGE_BINARY_LFC2$primary_tumour_type %in% c("Colorectal", "Oesophageal", "Gastric", "Pancreatic"),]

## Create cancer type-specific subsets
MERGE_BINARY_LFC2_COLO <- MERGE_BINARY_LFC2[MERGE_BINARY_LFC2$primary_tumour_type == "Colorectal",]

MERGE_BINARY_LFC2_OESO <- MERGE_BINARY_LFC2[MERGE_BINARY_LFC2$primary_tumour_type == "Oesophageal",]



################################################################################
# 6. DIFFERENTIAL DEPENDENCY ANALYSIS - GASTROINTESTINAL
################################################################################

## Perform differential dependency analysis across all GI cancers
SUMMARY <- NULL

for (g in unique(MERGE_BINARY_LFC2$gene)) {
  
  print(g)
  
  # Subset data for current gene
  df_temp <- MERGE_BINARY_LFC2[MERGE_BINARY_LFC2$gene == g, ]
  id <- unique(df_temp[, "ensembl_id"])
  
  # Calculate overall mean and median LFC
  MEAN_LFC <- mean(df_temp$LFC)
  MEDIAN_LFC <- median(df_temp$LFC)
  
  # Check if gene has both depleted and non-depleted samples
  N <- length(unique(df_temp$is_depleted))
  
  if (N == 2) {  # Proceed only if variability exists
    
    # Calculate statistics by depletion status
    df_temp2 <- df_temp %>% 
      group_by(is_depleted) %>% 
      dplyr::summarise(
        MEAN_LFC = mean(LFC),
        MEDIAN_LFC = median(LFC)
      )
    
    # Reshape to wide format
    df_temp2 <- pivot_wider(
      data = df_temp2, 
      names_from = is_depleted,
      id_cols = NULL,
      values_from = c("MEAN_LFC", "MEDIAN_LFC")
    )
    
    # Calculate effect sizes (delta between groups)
    df_temp2$delta_MEAN <- df_temp2[[2]] - df_temp2[[1]]
    df_temp2$delta_MEDIAN <- df_temp2[[4]] - df_temp2[[3]]
    
    # Count samples in each group
    N_not_depleted <- nrow(df_temp[df_temp$is_depleted == "0", ])
    N_depleted <- nrow(df_temp[df_temp$is_depleted == "1", ])
    
    # Perform Wilcoxon test if both groups have >1 sample
    if (N_depleted > 1 & N_not_depleted > 1) {
      
      W.TEST <- wilcox.test(df_temp$LFC ~ as.factor(df_temp$is_depleted))
      P.VALUE <- W.TEST$p.value
      
    } else {
      
      P.VALUE <- NA
      
    }
    
    # Save results
    SUMMARY <- rbind(
      SUMMARY, 
      data.frame(g, id, MEAN_LFC, MEDIAN_LFC, df_temp2, N_not_depleted, N_depleted, P.VALUE)
    )
  }
}

# Rename columns for clarity
colnames(SUMMARY) <- c(
  "gene", "ensembl_id", "MEAN_LFC", "MEDIAN_LFC", 
  "MEAN_LFC_not_depleted_group", "MEAN_LFC_depleted_group", 
  "MEDIAN_LFC_not_depleted_group", "MEDIAN_LFC_depleted_group", 
  "delta_MEAN", "delta_MEDIAN", 
  "N_not_depleted", "N_depleted", "p.value"
)

## Apply quality filters
# Remove genes with low expression and require minimum sample counts
SUMMARY_filter <- SUMMARY[!SUMMARY$gene %in% RNAseq_not_expressed$gene & SUMMARY$N_not_depleted > 1 & SUMMARY$N_depleted > 1,]

# Apply FDR correction
SUMMARY_filter$p.value.adj <- p.adjust(SUMMARY_filter$p.value, method = "fdr")

# Filter for significant genes (FDR < 0.05)
SUMMARY_filter2 <- SUMMARY_filter[SUMMARY_filter$p.value.adj < 0.05, ]

# NOTE: Save SUMMARY_filter as supplementary_table_7.1



################################################################################
# 7. DIFFERENTIAL DEPENDENCY ANALYSIS - COLORECTAL
################################################################################

## Repeat analysis for colorectal cancers specifically
SUMMARY_COLO <- NULL

for (g in unique(MERGE_BINARY_LFC2_COLO$gene)) {
  
  print(g)
  
  df_temp <- MERGE_BINARY_LFC2_COLO[MERGE_BINARY_LFC2_COLO$gene == g, ]
  id <- unique(df_temp[, "ensembl_id"])
  
  MEAN_LFC <- mean(df_temp$LFC)
  MEDIAN_LFC <- median(df_temp$LFC)
  
  N <- length(unique(df_temp$is_depleted))
  
  if (N == 2) {
    
    df_temp2 <- df_temp %>% 
      group_by(is_depleted) %>% 
      dplyr::summarise(
        MEAN_LFC = mean(LFC),
        MEDIAN_LFC = median(LFC)
      )
    
    df_temp2 <- pivot_wider(
      data = df_temp2, 
      names_from = is_depleted,
      id_cols = NULL,
      values_from = c("MEAN_LFC", "MEDIAN_LFC")
    )
    
    df_temp2$delta_MEAN <- df_temp2[[2]] - df_temp2[[1]]
    df_temp2$delta_MEDIAN <- df_temp2[[4]] - df_temp2[[3]]
    
    N_not_depleted <- nrow(df_temp[df_temp$is_depleted == "0", ])
    N_depleted <- nrow(df_temp[df_temp$is_depleted == "1", ])
    
    if (N_depleted > 1 & N_not_depleted > 1) {
      
      W.TEST <- wilcox.test(df_temp$LFC ~ as.factor(df_temp$is_depleted))
      P.VALUE <- W.TEST$p.value
      
    } else {
      
      P.VALUE <- NA
      
    }
    
    SUMMARY_COLO <- rbind(
      SUMMARY_COLO, 
      data.frame(g, id, MEAN_LFC, MEDIAN_LFC, df_temp2, N_not_depleted, N_depleted, P.VALUE)
    )
  }
}

colnames(SUMMARY_COLO) <- c(
  "gene", "ensembl_id", "MEAN_LFC", "MEDIAN_LFC", 
  "MEAN_LFC_not_depleted_group", "MEAN_LFC_depleted_group", 
  "MEDIAN_LFC_not_depleted_group", "MEDIAN_LFC_depleted_group", 
  "delta_MEAN", "delta_MEDIAN", 
  "N_not_depleted", "N_depleted", "p.value"
)

SUMMARY_COLO_filter <- SUMMARY_COLO[!SUMMARY_COLO$gene %in% RNAseq_not_expressed_COLO$gene & SUMMARY_COLO$N_not_depleted > 1 & SUMMARY_COLO$N_depleted > 1,]

SUMMARY_COLO_filter$p.value.adj <- p.adjust(SUMMARY_COLO_filter$p.value, method = "fdr")

SUMMARY_COLO_filter2 <- SUMMARY_COLO_filter[SUMMARY_COLO_filter$p.value.adj < 0.05, ]

# NOTE: Save SUMMARY_COLO_filter as supplementary_table_7.2



################################################################################
# 8. DIFFERENTIAL DEPENDENCY ANALYSIS - OESOPHAGEAL
################################################################################

## Repeat analysis for oesophageal cancers specifically
SUMMARY_OESO <- NULL

for (g in unique(MERGE_BINARY_LFC2_OESO$gene)) {
  
  print(g)
  
  df_temp <- MERGE_BINARY_LFC2_OESO[MERGE_BINARY_LFC2_OESO$gene == g, ]
  id <- unique(df_temp[, "ensembl_id"])
  
  MEAN_LFC <- mean(df_temp$LFC)
  MEDIAN_LFC <- median(df_temp$LFC)
  
  N <- length(unique(df_temp$is_depleted))
  
  if (N == 2) {
    
    df_temp2 <- df_temp %>% 
      group_by(is_depleted) %>% 
      dplyr::summarise(
        MEAN_LFC = mean(LFC),
        MEDIAN_LFC = median(LFC)
      )
    
    df_temp2 <- pivot_wider(
      data = df_temp2, 
      names_from = is_depleted,
      id_cols = NULL,
      values_from = c("MEAN_LFC", "MEDIAN_LFC")
    )
    
    df_temp2$delta_MEAN <- df_temp2[[2]] - df_temp2[[1]]
    df_temp2$delta_MEDIAN <- df_temp2[[4]] - df_temp2[[3]]
    
    N_not_depleted <- nrow(df_temp[df_temp$is_depleted == "0", ])
    N_depleted <- nrow(df_temp[df_temp$is_depleted == "1", ])
    
    if (N_depleted > 1 & N_not_depleted > 1) {
      
      W.TEST <- wilcox.test(df_temp$LFC ~ as.factor(df_temp$is_depleted))
      P.VALUE <- W.TEST$p.value
      
    } else {
      
      P.VALUE <- NA
      
    }
    
    SUMMARY_OESO <- rbind(
      SUMMARY_OESO, 
      data.frame(g, id, MEAN_LFC, MEDIAN_LFC, df_temp2, N_not_depleted, N_depleted, P.VALUE)
    )
  }
}

colnames(SUMMARY_OESO) <- c(
  "gene", "ensembl_id", "MEAN_LFC", "MEDIAN_LFC", 
  "MEAN_LFC_not_depleted_group", "MEAN_LFC_depleted_group", 
  "MEDIAN_LFC_not_depleted_group", "MEDIAN_LFC_depleted_group", 
  "delta_MEAN", "delta_MEDIAN", 
  "N_not_depleted", "N_depleted", "p.value"
)

SUMMARY_OESO_filter <- SUMMARY_OESO[!SUMMARY_OESO$gene %in% RNAseq_not_expressed_OESO$gene & SUMMARY_OESO$N_not_depleted > 1 & SUMMARY_OESO$N_depleted > 1,]

SUMMARY_OESO_filter$p.value.adj <- p.adjust(SUMMARY_OESO_filter$p.value, method = "fdr")

SUMMARY_OESO_filter2 <- SUMMARY_OESO_filter[SUMMARY_OESO_filter$p.value.adj < 0.05, ]

# NOTE: Save SUMMARY_OESO_filter as supplementary_table_7.3
