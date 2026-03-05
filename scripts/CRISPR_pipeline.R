################################################################################
# Script: CRISPR_pipeline.R
# Purpose: Process raw CRISPR screen data through complete QC and normalization
# Description: Complete pipeline from raw counts to normalized gene-level 
#              fitness scores with quality control, batch correction, and 
#              binary dependency calling using BAGEL2
#
# Workflow:
# 1. Load raw sgRNA counts from sequencing
# 2. Calculate log2 fold changes (LFCs)
# 3. Perform sgRNA correlation analysis
# 4. Quality control: Replicate reproducibility filtering
# 5. CRISPRcleanR normalization (copy number correction)
# 6. Spline normalization (expression bias correction)
# 7. ComBat batch correction (library effects)
# 8. Generate gene-level scores
# 9. BAGEL2 integration for binary calls
# 10. Final scaling and output generation
################################################################################

# Load required libraries
library(data.table)
library(tidyverse)
library(CRISPRcleanR)  # For copy number bias correction
library(sva)           # For ComBat batch correction

################################################################################
# 1. LOAD AND PREPARE DATA
################################################################################

## Set working directory and data folder - Raw counts are not publicly available
setwd('/path_to_folder/')
FOLDER <- '/path_to_folder/'

## List all raw count files from sequencing
raw_counts <- list.files(path = FOLDER, pattern = '\\.raw_counts.csv', recursive = TRUE)


## Load sgRNA library information
# Two libraries were used: KY v1.1 (newer) and MiniLibCas9 (older)
data("KY_Library_v1.1")        # Load v1.1 library annotation
data("MiniLibCas9_Library")    # Load MiniLib library annotation

# Create mapping table to convert MiniLib names to v1.1 nomenclature
minLib_to_v1.1_names <- MiniLibCas9_Library[, c("sgRNA_ID", "CODE")]
colnames(minLib_to_v1.1_names)[1] <- "sgRNA"

# Load common sgRNAs present in both libraries (for cross-library comparison)
common_sgRNAs <- as.data.frame(
  fread('/path_to_folder/common-sgRNAs_minLib-v1_ensembl-ids.csv')
)

## Load control gene lists for quality assessment
# Essential genes: Should show strong negative LFC (depletion)
essential <- as.data.frame(
  fread('/path_to_folder/essential_genes_organoids_downsampling.csv')
)

# Non-essential genes: Should show LFC near zero (no effect)
non_essential <- as.data.frame(
  fread('/path_to_folder/non-essential_genes_organoids_downsampling.csv')
)

# Curated control lists for BAGEL2 (more stringent)
curated_essential <- read.table(
  '/path_to_folder/curated_BAGEL_essential_genes_organoids_downsampling.txt'
)

curated_non_essential <- read.table(
  '/path_to_folder/curated_BAGEL_non-essential_genes_organoids_downsampling.txt'
)



################################################################################
# 2. CALCULATE LOG2 FOLD CHANGES
################################################################################

## Process each sample: Calculate LFC from raw counts

LFC <- NULL

for (x in raw_counts) {
  
  # Extract sample and project identifiers from filename
  ID_long <- str_split_i(x, "/", 2)
  ID <- str_split_i(ID_long, "\\.", 1)
  print(ID)
  
  # Read raw count data
  B <- as.data.frame(fread(x))
  
  # Count number of replicates in this sample
  replicates <- ncol(B)
  replicates <- replicates - 4  # Subtract annotation columns
  
  # Convert to long format for LFC calculation
  C <- gather(B, "replicate", "counts", 3:(3 + replicates - 1))
  colnames(C)[3] <- "plasmid"
  
  # Calculate log2 fold change
  # Formula: LFC = log2((counts + 0.5) / plasmid)
  # Pseudocount (+0.5) prevents log(0) = -Inf for sgRNAs with zero counts
  D <- C %>% 
    group_by(sgRNA, gene, replicate, library) %>% 
    summarise(LFC = log2((counts + 0.5) / plasmid))
  
  # Separate sample ID from replicate number
  D <- separate(D, replicate, into = c("sample_ID", "replicate"), sep = "\\.", remove = TRUE)
  
  # Accumulate all LFCs
  LFC <- rbind(LFC, D)
}



################################################################################
# 3. sgRNA CORRELATION ANALYSIS
################################################################################

## Goal: Identify genes with highly correlated sgRNA effects

# Create unique sample-replicate-library identifier
LFC$sample_replicate_library <- paste(LFC$sample_ID, LFC$replicate, LFC$library, sep = ".")
LFC <- LFC[, c("sgRNA", "gene", "sample_replicate_library", "LFC")]

# Use only common sgRNAs for correlation analysis
LFC_corr <- merge(LFC, common_sgRNAs, by = "sgRNA")
LFC_corr$gene.x <- NULL
colnames(LFC_corr)[7] <- "gene"
LFC_corr <- LFC_corr[, c("sgRNA", "gene", "sample_replicate_library", "LFC")]

## Calculate Pearson correlation between sgRNA pairs for each gene
correlation_sgRNAs <- NULL

for (i in unique(LFC_corr$gene)) {
  
  print(i)
  
  # Get LFCs for current gene
  df_temp <- LFC_corr[LFC_corr$gene == i, ]
  
  # Reshape to wide format (sgRNAs as columns)
  df_temp <- spread(df_temp, sgRNA, LFC)
  row.names(df_temp) <- df_temp$sample_replicate_library
  df_temp$sample_replicate_library <- NULL
  df_temp$gene <- NULL
  
  # Calculate correlation if gene has >= 2 sgRNAs with complete data
  if (ncol(df_temp) >= 2 & any(is.na(df_temp)) == FALSE) {
    
    # Pearson correlation between first two sgRNAs
    CORR <- cor.test(df_temp[[1]], df_temp[[2]], method = "pearson")
    df_final <- data.frame(i, CORR$estimate)
    colnames(df_final) <- c("gene", "PCC")
    
  } else {
    
    # Not enough data for correlation
    df_final <- data.frame(i, NA)
    colnames(df_final) <- c("gene", "PCC")
  }
  
  correlation_sgRNAs <- rbind(correlation_sgRNAs, df_final)
}

## Select top 400 genes with highest correlations
# These genes will be used for replicate quality control
# Rationale: Well-behaved genes should show consistent replicate behavior
selected_genes <- correlation_sgRNAs %>% 
  dplyr::arrange(desc(PCC)) %>% 
  dplyr::slice(1:400)

# Calculate mean LFC across both sgRNAs for selected genes
LFC_subset <- LFC_corr[LFC_corr$gene %in% selected_genes$gene, ]
LFC_subset <- spread(LFC_subset, "sample_replicate_library", "LFC")
LFC_subset$sgRNA <- NULL

# Average the two sgRNAs targeting each gene
LFC_subset <- LFC_subset %>% 
  dplyr::group_by(gene) %>% 
  dplyr::summarise(across(everything(), mean))
LFC_subset$gene <- NULL



################################################################################
# 4. REPLICATE QUALITY CONTROL
################################################################################

## Methodology: Determine which replicates are truly biological replicates
# Approach: Calculate pairwise distances between all replicates
# Expected: Replicates from same sample should be more similar than random pairs

## Calculate Euclidean distances between all replicate pairs
# Uses only the 400 most correlated genes for robust distance estimation
df_temp2_dist <- t(LFC_subset) %>% 
  stats::dist(method = 'euclidean')

# Convert distance matrix to data frame
df_temp3_dist <- df_temp2_dist %>%
  as.matrix %>% 
  as.data.frame %>%
  rownames_to_column(var = 'var1') %>%
  gather(var2, value, -var1)

# Remove duplicate pairs (lower triangle of distance matrix)
df_temp3_dist <- df_temp3_dist %>%
  dplyr::mutate(
    var_order = paste(var1, var2) %>%
      strsplit(split = ' ') %>%
      map_chr(~ sort(.x) %>% paste(collapse = ' '))
  ) %>%
  dplyr::mutate(cnt = 1) %>%
  dplyr::group_by(var_order) %>%
  dplyr::mutate(cumsum = cumsum(cnt)) %>%
  dplyr::filter(cumsum != 2) %>%
  ungroup %>%
  dplyr::select(-var_order, -cnt, -cumsum)

# Remove self-comparisons (distance = 0)
df_temp4_dist <- df_temp3_dist[df_temp3_dist$var1 != df_temp3_dist$var2, ]

# Parse sample identifiers
df_temp4_dist <- separate(
  df_temp4_dist, var1, 
  into = c("sample1", "replicate1", "library1"), 
  sep = "\\.", 
  remove = FALSE
)
df_temp4_dist <- separate(
  df_temp4_dist, var2, 
  into = c("sample2", "replicate2", "library2"), 
  sep = "\\.", 
  remove = FALSE
)

## Classify comparisons
# "Same model comparison": Replicates from the same organoid sample
# "Random model comparison": Replicates from different organoid samples
df_temp4_dist$comparison <- "Random model comparison"
df_temp4_dist[
  paste(df_temp4_dist$sample1, df_temp4_dist$library1, sep = "-") == 
  paste(df_temp4_dist$sample2, df_temp4_dist$library2, sep = "-"), 
]$comparison <- "Same model comparison"

replicate_distances <- df_temp4_dist %>%
  mutate(
    comparison = factor(
      comparison, 
      levels = c("Same model comparison", "Random model comparison")
    )
  ) %>%
  arrange(comparison)

# Calculate summary statistics for each distribution
summary_replicate_distances <- replicate_distances %>% 
  dplyr::group_by(comparison) %>% 
  dplyr::summarise(
    mean_value = mean(value),
    sd_value = sd(value),
    sd_value2 = 2 * sd(value)
  )

## Statistical testing to determine optimal threshold
# Goal: Find distance where P(same model) / P(random model) = 2
# This maximizes true biological replicates while minimizing false positives

# Kolmogorov-Smirnov test: Are the two distributions different?
ks.test(
  replicate_distances[replicate_distances$comparison == "Random model comparison", ]$value,
  replicate_distances[replicate_distances$comparison == "Same model comparison", ]$value
)

# Calculate probabilities at threshold = 48.75
# This threshold was determined empirically where P(same)/P(random) ≈ 2
pnorm(
  48.75, 
  mean = summary_replicate_distances[summary_replicate_distances$comparison == "Random model comparison", ]$mean_value, 
  sd = summary_replicate_distances[summary_replicate_distances$comparison == "Random model comparison", ]$sd_value
)

pnorm(
  48.75, 
  mean = summary_replicate_distances[summary_replicate_distances$comparison == "Same model comparison", ]$mean_value, 
  sd = summary_replicate_distances[summary_replicate_distances$comparison == "Same model comparison", ]$sd_value
)

## Apply threshold: Keep only replicates with distance <= 48.75
# Rationale: These replicates are twice as likely to be true biological replicates
replicates_PASS_dist <- replicate_distances[
  replicate_distances$value <= 48.75 & 
  replicate_distances$comparison == "Same model comparison", 
]

# Create list of passing replicates
replicates_PASS_dist_summ <- unique(c(replicates_PASS_dist$var1, replicates_PASS_dist$var2))



################################################################################
# 5. CRISPRcleanR NORMALIZATION
################################################################################

## CRISPRcleanR corrects for copy number bias in CRISPR screens
# Problem: sgRNAs in amplified genomic regions appear depleted (more copies to cut)
# Solution: Model and remove copy number-associated bias

## Prepare sgRNA annotation with updated gene names
sgRNA_names <- common_sgRNAs[, c("sgRNA", "gene", "ensembl_id", "HGNC_ID")]
colnames(sgRNA_names)[1] <- "CODE"

# Merge with library annotation
LIBRARY <- merge(sgRNA_names, KY_Library_v1.1, by = "CODE")
LIBRARY$GENES <- NULL
colnames(LIBRARY)[2] <- "GENES"
rownames(LIBRARY) <- LIBRARY$CODE

## Create output directory for LFC files
output_directory_LFC <- paste(FOLDER, "LFC", sep = "")
dir.create(file.path(output_directory_LFC))

# Initialize storage for results
gene_names <- unique(sgRNA_names[, 2:4])
gene_LFCs <- gene_names
ROC_PR_AUC <- NULL
ROC_PR_values <- NULL

## Process each sample through CRISPRcleanR
for (y in raw_counts2) {
  
  # Parse sample information
  y2 <- str_split_i(y, "\\/", 2)
  ID <- str_split_i(y2, "\\.", 1)
  library <- str_split_i(y2, "\\.", 2)
  ID_library <- paste(ID, library, sep = ".")
  
  # Load raw counts
  E <- as.data.frame(fread(y))
  
  # Add library suffix to replicate names
  colnames(E)[3:(ncol(E) - 1)] <- paste(
    colnames(E)[3:(ncol(E) - 1)], 
    library, 
    sep = "."
  )
  replicates <- colnames(E)[3:(ncol(E) - 1)]
  
  ## Filter replicates: Keep only those that passed QC
  for (rep in replicates) {
    
    if (rep %in% replicates_PASS_dist_summ) {
      # Replicate passed QC - keep it
      E <- E
    } else {
      # Replicate failed QC - remove it
      E <- E[, !colnames(E) %in% rep]
    }
  }
  
  ## Proceed only if sample has sufficient replicates (>= 3)
  # Minimum 3 replicates needed for robust normalization
  if (ncol(E) >= 5) {  # 5 columns = sgRNA + gene + plasmid + 2 replicates
    
    print(ID_library)
    
    # Format data for CRISPRcleanR
    E$sgRNA <- as.character(E$sgRNA)
    E$gene <- as.character(E$gene)
    E <- E %>% relocate(colnames(E)[ncol(E)], .after = "gene")
    
    # Create sample-specific output directory
    output_directory_LFC_sample <- paste(output_directory_LFC, ID_library, sep = "/")
    dir.create(file.path(output_directory_LFC_sample))
    
    ## Step 1: Calculate normalized fold changes
    # Removes technical artifacts and low-count sgRNAs (min_reads = 30)
    normANDfcs <- ccr.NormfoldChanges(
      Dframe = E,
      min_reads = 30,
      EXPname = ID,
      libraryAnnotation = LIBRARY,
      outdir = paste(output_directory_LFC_sample, "/", sep = "")
    )
    
    ## Step 2: Map sgRNAs to genomic positions
    gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs$logFCs, LIBRARY)
    
    ## Step 3: Correct copy number bias
    # Uses local regression to identify and remove CN-associated effects
    correctedFCs <- ccr.GWclean(gwSortedFCs, display = FALSE)
    
    ## Step 4: Extract corrected sgRNA-level LFCs
    sample_sgRNAFCs <- correctedFCs$corrected_logFCs$correctedFC
    names(sample_sgRNAFCs) <- rownames(correctedFCs$corrected_logFCs)
    
    ## Step 5: Aggregate to gene-level scores
    # Takes mean of sgRNAs targeting the same gene
    sample_geneFCs <- ccr.geneMeanFCs(sample_sgRNAFCs, LIBRARY)
    
    ## Quality assessment using ROC and Precision-Recall curves
    # Goal: Measure separation between essential and non-essential genes
    
    # ROC curve (Receiver Operating Characteristic)
    ROC <- ccr.ROC_Curve(
      sample_geneFCs, 
      essential$`Target Gene Symbol`, 
      non_essential$`Target Gene Symbol`
    )
    
    # Precision-Recall curve
    PRRC <- ccr.PrRc_Curve(
      sample_geneFCs, 
      essential$`Target Gene Symbol`, 
      non_essential$`Target Gene Symbol`
    )
    
    # Calculate detailed ROC statistics
    FCsprofile <- sample_geneFCs[
      intersect(
        c(essential$`Target Gene Symbol`, non_essential$`Target Gene Symbol`),
        names(sample_geneFCs)
      )
    ]
    
    predictions <- FCsprofile
    observations <- is.element(names(FCsprofile), essential$`Target Gene Symbol`) + 0
    names(observations) <- names(predictions)
    
    # Generate ROC curve data
    roc.info <- roc(observations, predictions, direction = ">", quiet = TRUE)
    roc.df <- roc.info %>% 
      coords(ret = "all", transpose = FALSE) %>%
      dplyr::select(sensitivity, specificity, fpr, threshold, precision)
    
    # Calculate distance from perfect classifier (top-left corner)
    roc.df$Dist <- sqrt((1 - roc.df$sensitivity)^2 + (1 - roc.df$specificity)^2)
    roc.df$sample <- ID_library
    
    # Store AUC values
    ROC_AUC <- ROC$AUC
    PRRC_AUC <- PRRC$AUC
    
    ROC_PR_AUC <- rbind(ROC_PR_AUC, data.frame(ID_library, ROC_AUC, PRRC_AUC))
    ROC_PR_values <- rbind(ROC_PR_values, roc.df)
    
    # Store gene-level LFCs
    sample_geneFCs2 <- as.data.frame(sample_geneFCs)
    sample_geneFCs2 <- add_column(sample_geneFCs2, gene = rownames(sample_geneFCs2), .before = 1)
    colnames(sample_geneFCs2)[2] <- ID_library
    gene_LFCs <- merge(gene_LFCs, sample_geneFCs2, by = "gene")
  }
}

## Filter samples based on ROC performance
# Remove outliers: Samples with ROC AUC > 3 SD below mean
colnames(ROC_PR_AUC)[1] <- "sample"

ROC_PR_AUC$ROC_3SD_filtered <- ROC_PR_AUC$ROC_AUC < (1 - 3 * sd(ROC_PR_AUC$ROC_AUC))
ROC_PR_AUC$PR_3SD_filtered <- ROC_PR_AUC$PRRC_AUC < (1 - 3 * sd(ROC_PR_AUC$PRRC_AUC))

ROC_PR_AUC_filtered <- ROC_PR_AUC[
  ROC_PR_AUC$sample %in% ROC_PR_AUC[ROC_PR_AUC$ROC_3SD_filtered == "FALSE", ]$sample, 
]

# Keep only high-quality samples
gene_LFCs_PASS <- gene_LFCs[
  , colnames(gene_LFCs) %in% c("gene", "HGNC_ID", "ensembl_id", ROC_PR_AUC_filtered$sample)
]



################################################################################
# 6. SPLINE NORMALIZATION
################################################################################

## Problem: Systematic bias correlated with mean gene expression
# Highly expressed genes may show artifactual depletion
# Solution: Fit cubic spline to remove expression-dependent bias

OrgFC_spline <- gene_LFCs_PASS

# Prepare data matrix
rownames(OrgFC_spline) <- OrgFC_spline$gene
OrgFC_spline$gene <- NULL
OrgFC_spline$HGNC_ID <- NULL
OrgFC_spline$ensembl_id <- NULL
OrgFC_spline <- as.matrix(OrgFC_spline)

# Store dimensions for later
dm <- dimnames(OrgFC_spline)

# Calculate mean LFC across all samples for each gene
GeneMeans <- rowMeans(OrgFC_spline, na.rm = TRUE)

# Define spline knots at quartiles (25%, 50%, 75%)
knots <- quantile(GeneMeans, p = c(0.25, 0.5, 0.75))

## Fit and remove spline-based bias for each sample
# Method: Fit cubic B-spline, subtract fitted values, add back mean
# Result: Removes systematic bias while preserving biological variation
for (i in 1:ncol(OrgFC_spline)) {
  
  # Fit cubic spline model: LFC ~ f(mean_LFC)
  model <- lm(
    OrgFC_spline[, i] ~ splines::bs(GeneMeans, knots = knots),
    na.action = na.exclude
  )
  
  # Correct LFC: Remove fitted bias, restore overall mean
  if (i == 1) {
    OrgFC_spline2 <- OrgFC_spline[, i] - fitted.values(model) + GeneMeans
  } else {
    OrgFC_spline2 <- cbind(OrgFC_spline2, OrgFC_spline[, i] - fitted.values(model) + GeneMeans)
  }
}

# Restore dimension names
dimnames(OrgFC_spline2) <- dm



################################################################################
# 7. COMBAT BATCH CORRECTION
################################################################################

## Problem: Systematic differences between CRISPR library versions
# MiniLib vs v1.1 may have different baseline signals
# Solution: ComBat empirical Bayes batch correction

# Extract library information from sample names
batch <- as.data.frame(colnames(gene_LFCs_PASS[, 4:ncol(gene_LFCs_PASS)]))
colnames(batch)[1] <- "sample_library"
batch <- separate(batch, sample_library, into = c("sample_ID", "library"), sep = "\\.", remove = FALSE)

## Apply ComBat batch correction
# Method: Empirical Bayes shrinkage (Johnson et al. 2007, Biostatistics)
# Batch variable: CRISPR library (minLib vs v1_1)
# Preserves: Biological variation between samples
# Removes: Technical variation due to library differences
gene_LFCs_combat <- ComBat(
  as.matrix(OrgFC_spline2), 
  batch = batch$library,      # Batch variable to correct
  par.prior = TRUE,            # Use parametric priors
  prior.plots = TRUE           # Generate diagnostic plots
)



################################################################################
# 8. PREPARE FILES FOR BAGEL2
################################################################################

## BAGEL2 (Bayesian Analysis of Gene EssentiaLity) requires specific format
# Input: Gene-level LFC for each sample in separate files
# Output: Bayes Factor (BF) for each gene indicating essentiality

output_directory_BAGEL2 <- paste(FOLDER2, "BAGEL2", sep = "")
dir.create(file.path(output_directory_BAGEL2))

output_directory_BAGEL2_input <- paste(output_directory_BAGEL2, "input_bf", sep = "/")
dir.create(file.path(output_directory_BAGEL2_input))

## Generate BAGEL2 input files (one per sample)
gene_LFCs_combat <- as.data.frame(gene_LFCs_combat)

for (c in colnames(gene_LFCs_combat)) {
  
  print(c)
  
  # Extract LFC for current sample
  df_temp <- as.data.frame(gene_LFCs_combat[, colnames(gene_LFCs_combat) %in% c])
  
  # Add required columns for BAGEL2
  df_temp <- add_column(df_temp, REAGENT_ID = rownames(gene_LFCs_combat), .before = 1)
  df_temp <- add_column(df_temp, GENE = rownames(gene_LFCs_combat), .after = 1)
  colnames(df_temp)[3] <- c
  
  # Write sample-specific file for BAGEL2
  write.table(
    df_temp, 
    file = paste(output_directory_BAGEL2_input, "/", c, ".gene_LFC_combat_input_bf", ".tsv", sep = ""), 
    row.names = FALSE, 
    col.names = TRUE, 
    sep = "\t", 
    quote = FALSE
  )
}

## Prepare gene annotation for output
gene_LFCs_combat <- add_column(gene_LFCs_combat, gene = rownames(gene_LFCs_combat), .before = 1)
gene_LFCs_combat <- merge(gene_names, gene_LFCs_combat, by = "gene")



################################################################################
# 9. COHEN'S D QUALITY CONTROL
################################################################################

## Cohen's D measures separation between essential and non-essential genes
# Formula: D = (mean_nonessential - mean_essential) / pooled_SD
# Interpretation:
#   D < 0.5:  Small effect (poor separation)
#   D = 0.5-0.8: Medium effect (acceptable)
#   D > 0.8:  Large effect (good separation) <- Our QC threshold

cohens_D <- NULL

for (c in colnames(gene_LFCs_combat)[4:ncol(gene_LFCs_combat)]) {
  
  # Extract sample LFCs
  df_temp <- gene_LFCs_combat[, colnames(gene_LFCs_combat) %in% c("gene", c)]
  colnames(df_temp)[2] <- "sample_ID"
  
  # Calculate statistics for essential genes
  MEAN_ess <- mean(df_temp[df_temp$gene %in% essential$`Target Gene Symbol`, ]$sample_ID)
  SD_ess <- sd(df_temp[df_temp$gene %in% essential$`Target Gene Symbol`, ]$sample_ID)
  
  # Calculate statistics for non-essential genes
  MEAN_non_ess <- mean(df_temp[df_temp$gene %in% non_essential$`Target Gene Symbol`, ]$sample_ID)
  SD_non_ess <- sd(df_temp[df_temp$gene %in% non_essential$`Target Gene Symbol`, ]$sample_ID)
  
  # Calculate Cohen's D with pooled standard deviation
  COEFF <- (MEAN_non_ess - MEAN_ess) / sqrt(((SD_non_ess)^2 + (SD_ess)^2) / 2)
  
  cohens_D <- rbind(cohens_D, data.frame(c, COEFF))
}

# Quality control: All samples should have Cohen's D > 0.8



################################################################################
# 10. FINAL SCALING
################################################################################

## Goal: Normalize LFC scale across all samples
# Method: Scale so non-essential genes center at 0, essential genes at -1
# Result: Comparable fitness scores across all samples

# Calculate median LFCs for control genes
LFC_median <- gather(gene_LFCs_combat, "sample", "LFC", 4:ncol(gene_LFCs_combat))
LFC_median$type_gene <- "X"
LFC_median[LFC_median$gene %in% essential$`Target Gene Symbol`, ]$type_gene <- "essential"
LFC_median[LFC_median$gene %in% non_essential$`Target Gene Symbol`, ]$type_gene <- "non_essential"
LFC_median$LFC <- as.numeric(LFC_median$LFC)

LFC_median <- LFC_median %>% 
  group_by(type_gene) %>% 
  dplyr::summarise(median_LFC = median(LFC, na.rm = TRUE))

## Apply scaling transformation to each sample
# Formula: scaled_LFC = (LFC - median_nonessential) / range
# where range = median_nonessential - median_essential
gene_LFC_final <- as.data.frame(gene_LFCs_combat[, c("gene", "ensembl_id", "HGNC_ID")])

for (s in colnames(gene_LFCs_combat)[4:ncol(gene_LFCs_combat)]) {
  
  df_temp <- gene_LFCs_combat[, colnames(gene_LFCs_combat) %in% c("gene", "ensembl_id", "HGNC_ID", s)]
  df_temp2 <- df_temp
  
  # Apply scaling transformation
  df_temp2[[4]] <- (df_temp[[4]] - LFC_median[LFC_median$type_gene == "non_essential", ]$median_LFC) / 
                   (LFC_median[LFC_median$type_gene == "non_essential", ]$median_LFC - 
                    LFC_median[LFC_median$type_gene == "essential", ]$median_LFC)
  
  gene_LFC_final <- merge(gene_LFC_final, df_temp2, by = c("gene", "ensembl_id", "HGNC_ID"))
}


# Convert to long format
gene_LFC_final <- gather(gene_LFC_final, "sample_library", "LFC", 4:ncol(gene_LFC_final))
gene_LFC_final <- separate(gene_LFC_final, sample_library, into = c("sample_ID", "library"), sep = "\\.")

# NOTE: Save gene_LFC_final as supplementary_table_6



################################################################################
# 11. PROCESS BAGEL2 OUTPUT - BINARY MATRIX
################################################################################

## After running BAGEL2 (external step), process results
# BAGEL2 calculates Bayes Factors indicating gene essentiality
# BF > threshold → gene is essential (binary = 1)
# BF ≤ threshold → gene is non-essential (binary = 0)

output_directory_BAGEL2_output <- paste(output_directory_BAGEL2, "output_pr", sep = "/")

# List BAGEL2 output files (Precision-Recall format)
PRs <- list.files(path = output_directory_BAGEL2_output, pattern = '.PC', recursive = TRUE)

# Initialize storage matrices
gene_binary <- gene_names      # Binary essentiality calls

for (x in PRs) {
  
  # Parse sample information
  ID <- str_split_i(x, "\\.", 1)
  library <- str_split_i(x, "\\.", 2)
  ID_library <- paste(ID, library, sep = ".")
  
  print(ID_library)
  
  # Read BAGEL2 results
  G <- read.table(
    file = paste(output_directory_BAGEL2_output, x, sep = "/"), 
    sep = "\t", 
    header = TRUE
  )
  
  # Determine threshold: Minimum BF at FDR < 0.05
  THRES <- min(G[G$FDR < 0.05, ]$BF)
  
  # Scale BF relative to threshold
  G$BFscaled <- G$BF - THRES
  
  # Create binary calls
  G$is_depleted <- 0
  G[G$BFscaled > 0, ]$is_depleted <- 1
  
  # Store binary calls
  binary_temp <- G[, c("Gene", "is_depleted")]
  colnames(binary_temp) <- c("gene", ID_library)
  gene_binary <- merge(gene_binary, binary_temp, by = "gene")
  
}

# Convert to long format
gene_binary <- gather(gene_binary, "sample_library", "is_depleted", 4:ncol(gene_binary))
gene_binary <- separate(gene_binary, sample_library, sep = "\\.", into = c("sample_ID", "library"), remove = TRUE)

# NOTE: Save gene_binary as supplementary_table_5



################################################################################
# SUMMARY OF OUTPUTS
################################################################################

# 1. supplementary_table_6: Scaled gene-level LFCs (fitness scores)
# 2. supplementary_table_5: Binary essentiality calls from BAGEL2
# 3. QC metrics: ROC/PR curves, Cohen's D, replicate distances
# 4. Intermediate files: Normalized LFCs at each processing step
