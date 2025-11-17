################################################################################
# Script: core_fitness_genes_analysis.R
# Purpose: Identify core fitness genes in organoid models
# Description: Uses ADAM (Adaptive Daisy Model) method to identify pan-cancer
#              core fitness genes specific to organoids, followed by pathway
#              enrichment analysis
################################################################################

# Load required libraries
library(data.table)
library(tidyverse)
library(CoRe)          # For ADAM analysis
library(msigdbr)       # For pathway database access
library(clusterProfiler)  # For gene ID conversion

################################################################################
# 1. LOAD INPUT DATA
################################################################################

## Load binary dependency matrix (gene depletion calls)
BINARY <- as.data.frame(
  fread('/path_to_folder/supplementary_table_5.csv')
)

## Load curated essential gene lists for validation
# BAGEL essential genes (positive controls)
curated_BAGEL_essential <- as.data.frame(
  fread('/path_to_folder/curated_BAGEL_essential_genes_organoids_downsampling.csv')
)

# BAGEL non-essential genes (negative controls)
curated_BAGEL_nonEssential <- as.data.frame(
  fread('/path_to_folder/curated_BAGEL_non-essential_genes_organoids_downsampling.csv')
)

# Additional control genes
control_genes <- as.data.frame(
  fread('/path_to_folder/control-genes_organoids.csv')
)
control_genes <- control_genes$gene

# Convert gene symbols to Ensembl IDs
control_genes <- bitr(
  control_genes, 
  fromType = "SYMBOL", 
  toType = "ENSEMBL", 
  OrgDb = "org.Hs.eg.db"
)

## Load pan-cancer core fitness genes from cell lines (ADAM method)
load('/path_to_folder/10_PANCANCER_coreFitness_genes.RData')

# Convert to Ensembl IDs
PanCancerCoreFitnessGenes <- bitr(
  PanCancerCoreFitnessGenes, 
  fromType = "SYMBOL", 
  toType = "ENSEMBL", 
  OrgDb = "org.Hs.eg.db"
)

# Keep only genes present in organoid analysis
PanCancerCoreFitnessGenes2 <- PanCancerCoreFitnessGenes[PanCancerCoreFitnessGenes$ENSEMBL %in% BINARY$ensembl_id,]

## Load pan-cancer common essential genes from cell lines (AUC method)
CommonEssentialsAUC <- as.data.frame(
  fread('/path_to_folder/CRISPRInferredCommonEssentials.csv')
)
CommonEssentialsAUC <- c(CommonEssentialsAUC$V1)

# Convert to Ensembl IDs
CommonEssentialsAUC <- bitr(
  CommonEssentialsAUC, 
  fromType = "SYMBOL", 
  toType = "ENSEMBL", 
  OrgDb = "org.Hs.eg.db"
)

# Keep only genes present in organoid analysis
CommonEssentialsAUC2 <- CommonEssentialsAUC[CommonEssentialsAUC$ENSEMBL %in% BINARY$ensembl_id,]

## Load cancer type annotations
cancer_type <- as.data.frame(
  fread('/path_to_folder/supplementary_table_2.2.csv')
)
cancer_type <- cancer_type[, c(1, 4)]
colnames(cancer_type) <- c("model_name", "tissue")  # Required names for CoRe



################################################################################
# 2. FILTER SAMPLES
################################################################################

## Remove duplicate samples and outliers
# Note: These samples were screened with multiple libraries or failed QC

# Remove duplicates
BINARY2 <- BINARY[!(BINARY$sample_ID == "HCM-SANG-0267-D12" & BINARY$library == "minLib"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "HCM-SANG-0268-C18" & BINARY2$library == "v1_1"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "HCM-SANG-0270-C20" & BINARY2$library == "minLib"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "HCM-SANG-0272-C20" & BINARY2$library == "v1_1"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "HCM-SANG-0273-C18" & BINARY2$library == "minLib"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "HCM-SANG-0276-C18" & BINARY2$library == "minLib"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "HCM-SANG-0282-C18" & BINARY2$library == "minLib"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "HCM-SANG-0529-C18" & BINARY2$library == "v1_1"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "HCM-SANG-0284-C18" & BINARY2$library == "minLib"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "HCM-SANG-0288-C18" & BINARY2$library == "v1_1"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "WTSI-COLO_278" & BINARY2$library == "minLib"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "WTSI-COLO_376" & BINARY2$library == "minLib"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "HCM-SANG-0291-C15" & BINARY2$library == "v1_1"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "HCM-SANG-0299-C15-A" & BINARY2$library == "minLib"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "HCM-SANG-0299-C15-B" & BINARY2$library == "minLib"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "HCM-SANG-0306-C15" & BINARY2$library == "minLib"), ]

# Remove samples with insufficient replicates per cancer type
# Gastric (only 3 samples)
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "WTSI-GAST_183" & BINARY2$library == "minLib"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "WTSI-GAST_187" & BINARY2$library == "minLib"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "WTSI-GAST_186" & BINARY2$library == "minLib"), ]

# Pancreatic (only 4 samples)
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "HCM-SANG-1095-C25" & BINARY2$library == "minLib"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "WTSI-PANC_067" & BINARY2$library == "minLib"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "HCM-SANG-1306-C25" & BINARY2$library == "minLib"), ]
BINARY2 <- BINARY2[!(BINARY2$sample_ID == "HCM-SANG-1320-C25" & BINARY2$library == "minLib"), ]

# Remove library column (no longer needed)
BINARY2$library <- NULL



################################################################################
# 3. PREPARE DATA FOR ADAM ANALYSIS
################################################################################

## Convert to wide format (genes as rows, samples as columns)
BINARY3 <- spread(BINARY2, sample_ID, is_depleted)

# Save gene annotation information
gene_names <- as.data.frame(
  unique(BINARY3[, c("gene", "ensembl_id", "HGNC_ID")])
)

# Set gene names as row names and remove annotation columns
rownames(BINARY3) <- BINARY3$gene
BINARY3$gene <- NULL
BINARY3$ensembl_id <- NULL
BINARY3$HGNC_ID <- NULL

## Prepare cancer type annotations for ADAM
# Filter for samples present in the binary matrix
cancer_type2 <- cancer_type[cancer_type$model_name %in% colnames(BINARY3), ]

# Map tissue types to cancer type names required by CoRe
cancer_type2$cancer_type <- "Colorectal Carcinoma"
cancer_type2[cancer_type2$tissue == "Oesophageal", ]$cancer_type <- "Esophageal Carcinoma"
cancer_type2[cancer_type2$tissue == "Ovarian", ]$cancer_type <- "Ovarian Carcinoma"



################################################################################
# 4. RUN ADAM ANALYSIS
################################################################################

## Identify pan-cancer core fitness genes using ADAM method
# ADAM = Adaptive Daisy Model
# This method identifies genes that are essential across multiple cancer types
PanCancerCoreFitnessGenes_ORG <- CoRe.PanCancer_ADaM(
  pancan_depMat = BINARY3,                           # Binary dependency matrix
  tissues_ctypes = c("Colorectal", "Oesophageal", "Ovarian"),  # Cancer types
  clannotation = cancer_type2,                       # Sample annotations
  TruePositives = curated_BAGEL_essential$gene,     # Positive control genes
  display = TRUE                                     # Show progress
)

# Convert to data frame
PanCancerCoreFitnessGenes_ORG <- as.data.frame(PanCancerCoreFitnessGenes_ORG)
colnames(PanCancerCoreFitnessGenes_ORG) <- "gene"

# Add gene identifiers
PanCancerCoreFitnessGenes_ORG <- merge(
  PanCancerCoreFitnessGenes_ORG, 
  gene_names
)



################################################################################
# 5. CLASSIFY ORGANOID-SPECIFIC GENES
################################################################################

## Label genes based on presence in cell line core fitness lists

# Check if gene is exclusive to organoids (not in cell line ADAM list)
PanCancerCoreFitnessGenes_ORG$exclusive_organoids <- "Yes"
PanCancerCoreFitnessGenes_ORG[PanCancerCoreFitnessGenes_ORG$ensembl_id %in% PanCancerCoreFitnessGenes2$ENSEMBL,]$exclusive_organoids <- "No"

# Check if gene is exclusive when excluding AUC common essentials
PanCancerCoreFitnessGenes_ORG$exclusive_organoids_exclude_CommonEssentials <- PanCancerCoreFitnessGenes_ORG$exclusive_organoids

PanCancerCoreFitnessGenes_ORG[PanCancerCoreFitnessGenes_ORG$ensembl_id %in% CommonEssentialsAUC2$ENSEMBL,]$exclusive_organoids_exclude_CommonEssentials <- "No"

# NOTE: Save PanCancerCoreFitnessGenes_ORG as supplementary_table_4.2



################################################################################
# 6. PATHWAY ENRICHMENT ANALYSIS
################################################################################

## Load pathway databases from MSigDB

# Biological Process (GO-BP) gene sets
BP_db <- msigdbr(
  species = "Homo sapiens", 
  category = "C5", 
  subcategory = "GO:BP"
)
BP_conv <- unique(BP_db[, c("gene_symbol", "gs_name", "ensembl_gene")])
BP_conv$type <- "BP"

# KEGG pathways
KEGG_db <- msigdbr(
  species = "Homo sapiens", 
  category = "C2", 
  subcategory = "CP:KEGG_LEGACY"
)
KEGG_conv <- unique(KEGG_db[, c("gene_symbol", "gs_name", "ensembl_gene")])
KEGG_conv$type <- "KEGG"

# Hallmark gene sets
H_db <- msigdbr(
  species = "Homo sapiens", 
  category = "H"
)
H_conv <- unique(H_db[, c("gene_symbol", "gs_name", "ensembl_gene")])
H_conv$type <- "H"

# Combine all pathway databases
all_pathways <- rbind(BP_conv, KEGG_conv, H_conv)

# Keep only pathways with genes in the CRISPR screen
all_pathways <- all_pathways[all_pathways$ensembl_gene %in% BINARY2$ensembl_id, ]
all_genes <- unique(BINARY2$ensembl_id)


################################################################################
# 7. ENRICHMENT ANALYSIS FOR ORGANOID-SPECIFIC GENES
################################################################################

## Test enrichment of organoid-exclusive genes in each pathway

# Select organoid-exclusive core fitness genes
exclusive_organoids <- PanCancerCoreFitnessGenes_ORG[PanCancerCoreFitnessGenes_ORG$exclusive_organoids_exclude_CommonEssentials == "Yes",]$ensembl_id

# Perform Fisher's exact test for each pathway
results_PANCAN_ORG_spec <- NULL
for (i in unique(all_pathways$gs_name)) {
  
  # Define background genes (all genes except organoid-specific)
  all_genes_rmv <- all_genes[!all_genes %in% exclusive_organoids]
  
  # Get genes in current pathway
  set <- all_pathways[all_pathways$gs_name %in% i,]
  TYPE <- unique(set$type)
  
  # Build 2x2 contingency table:
  # N1: Organoid-specific genes in pathway
  N1 <- length(set[set$ensembl_gene %in% exclusive_organoids,]$ensembl_gene)
  
  # Get gene names for reporting
  PATHWAY_GENES <- paste0(
    set[set$ensembl_gene %in% exclusive_organoids,]$gene_symbol, 
    collapse = ","
  )
  
  # N2: Non-organoid-specific genes in pathway
  N2 <- length(set[set$ensembl_gene %in% all_genes_rmv, ]$ensembl_gene)
  
  # Total genes in pathway
  N_PATHWAY <- N1 + N2
  
  # N3: Organoid-specific genes NOT in pathway
  N3 <- length(exclusive_organoids[!exclusive_organoids %in% set$ensembl_gene])
  
  # N4: Non-organoid-specific genes NOT in pathway
  N4 <- length(all_genes_rmv[!all_genes_rmv %in% set$ensembl_gene])
  
  # Create contingency table
  DF <- data.frame(
    "ORG_YES" = c(N1, N3),
    "ORG_NO" = c(N2, N4),
    row.names = c("PATHWAY_YES", "PATHWAY_NO"),
    stringsAsFactors = FALSE
  )
  
  # Perform Fisher's exact test
  FT <- fisher.test(DF)
  
  # Save results
  results_PANCAN_ORG_spec <- rbind(
    results_PANCAN_ORG_spec, 
    data.frame(i, TYPE, N1, PATHWAY_GENES, N_PATHWAY, FT$estimate, FT$p.value)
  )
}

## Apply FDR correction separately for each pathway database
results_PANCAN_ORG_spec_FDR <- NULL
for (j in unique(results_PANCAN_ORG_spec$TYPE)) {
  
  # Subset by pathway type
  results_PANCAN_ORG_spec_temp <- results_PANCAN_ORG_spec[results_PANCAN_ORG_spec$TYPE == j,]
  
  # Apply Benjamini-Hochberg FDR correction
  results_PANCAN_ORG_spec_temp$FDR <- p.adjust(
    results_PANCAN_ORG_spec_temp$FT.p.value, 
    method = "fdr"
  )
  
  results_PANCAN_ORG_spec_FDR <- rbind(
    results_PANCAN_ORG_spec_FDR, 
    results_PANCAN_ORG_spec_temp
  )
}

## Filter for significant enrichments
# Odds ratio > 2 and FDR < 0.05
results_PANCAN_ORG_spec_FDR_sig <- results_PANCAN_ORG_spec_FDR[results_PANCAN_ORG_spec_FDR$FT.estimate > 2 & results_PANCAN_ORG_spec_FDR$FDR < 0.05,]

# NOTE: Save results_PANCAN_ORG_spec_FDR_sig as supplementary_table_4.3
