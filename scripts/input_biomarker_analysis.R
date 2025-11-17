################################################################################
# Script: input_biomarker_analysis.R
# Purpose: Prepare biomarker input data for CRISPR biomarker analysis
# Description: This script processes genomic, transcriptomic, and clinical data
#              from organoid models to create binary and continuous feature 
#              matrices for downstream biomarker analysis
################################################################################

# Load required libraries
library(data.table)
library(tidyverse)

################################################################################
# 1. LOAD REFERENCE DATA
################################################################################

# Load cancer driver genes reference
cancer_driver_genes <- as.data.frame(
  fread('/path_to_folder/cancer_driver_genes.tsv')
)



################################################################################
# 2. TIDY DRIVER GENOMIC ALTERATIONS
################################################################################

# Load genomic alterations (mutations, CNVs, and structural variants)
mut_CNV_SV <- as.data.frame(
  fread('/path_to_folder/supplementary_table_3.4.csv')
)

# Filter for organoid models only
mut_CNV_SV <- mut_CNV_SV[mut_CNV_SV$model == "organoid", ]

# Create sample ID reference list
SAMPLES <- as.data.frame(unique(mut_CNV_SV[, "sample_ID"]))
colnames(SAMPLES) <- "sample_ID"

# Standardize genomic alteration nomenclature
mut_CNV_SV2 <- mut_CNV_SV %>% 
  dplyr::mutate(Genomic_alterations = case_when(
    Genomic_alterations == "AMP" ~ "Amplification",
    Genomic_alterations == "DEL" ~ "Large deletion",
    Genomic_alterations == "HOM_DISRUPTION" ~ "SV",
    Genomic_alterations == "PARTIAL_AMP" ~ "Amplification",
    .default = Genomic_alterations
  ))



################################################################################
# 3. PROCESS GENERAL MUTATIONS
################################################################################

# Extract mutations (exclude amplifications, deletions, and SVs)
mut_ALL <- mut_CNV_SV2[!mut_CNV_SV2$Genomic_alterations %in% c("Amplification", "Large deletion", "SV"), ]

# Filter out monoallelic mutations in loss-of-function (LoF) driver genes
# Keep only biallelic events for LoF genes
mut_ALL <- mut_ALL[!(mut_ALL$gene %in% cancer_driver_genes[cancer_driver_genes$method_of_action == "LoF", ]$gene & mut_ALL$biallelic_corrected == "No"), ]

# Create binary mutation matrix
mut_ALL2 <- mut_ALL[, c("sample_ID", "gene")]
mut_ALL2 <- unique(mut_ALL2)
mut_ALL2$count <- 1

# Reshape to wide format
mut_ALL3 <- spread(mut_ALL2, gene, count)
mut_ALL3 <- merge(SAMPLES, mut_ALL3, by = "sample_ID", all = TRUE)
mut_ALL3[is.na(mut_ALL3)] <- 0

# Add mutation suffix to column names
colnames(mut_ALL3)[2:ncol(mut_ALL3)] <- paste(
  colnames(mut_ALL3)[2:ncol(mut_ALL3)], 
  "mut", 
  sep = "_"
)

# Create combined BRAF/KRAS mutation indicator
mut_ALL3$BRAF_KRAS_mut <- mut_ALL3$BRAF_mut + mut_ALL3$KRAS_mut

# NOTE: Save mut_ALL3 as data_driver-mutations_ALL-organoids.csv



################################################################################
# 4. PROCESS SPECIFIC DRIVER VARIANTS (KRAS/BRAF)
################################################################################

# Extract specific variant information for KRAS and BRAF
driver_variants_ALL <- mut_ALL[, c("sample_ID", "gene", "variant")]

# Remove variants with unknown information (marked with "?")
driver_variants_ALL <- driver_variants_ALL[-grep("\\?", driver_variants_ALL$variant),]

# Extract amino acid position from variant notation (e.g., G12D -> G12)
driver_variants_ALL$variant2 <- gsub(
  '^([A-Z]{1}[0-9]+)(.+)$', 
  '\\1\\_\\2', 
  driver_variants_ALL$variant
)
driver_variants_ALL$variant2 <- gsub('\\_(.+)$', '', driver_variants_ALL$variant2)

# Focus on KRAS and BRAF variants only
driver_variants_ALL <- driver_variants_ALL[
  driver_variants_ALL$gene %in% c("KRAS", "BRAF"), 
]

# Filter variants that appear in at least 2 samples
driver_variants_ALL_count <- as.data.frame(
  table(driver_variants_ALL$variant, driver_variants_ALL$variant2)
)
driver_variants_ALL_count <- driver_variants_ALL_count[driver_variants_ALL_count$Freq > 0,]
driver_variants_ALL_count <- as.data.frame(table(driver_variants_ALL_count$Var2))

# Reshape data for filtering
driver_variants_ALL <- gather(driver_variants_ALL, "extra", "variant", 3:4)
driver_variants_ALL <- driver_variants_ALL[!driver_variants_ALL$variant %in% driver_variants_ALL_count[driver_variants_ALL_count$Freq < 2, ]$Var1,]

# Create gene_variant identifier
driver_variants_ALL$variant <- paste(
  driver_variants_ALL$gene, 
  driver_variants_ALL$variant, 
  sep = "_"
)
driver_variants_ALL$gene <- NULL
driver_variants_ALL$extra <- NULL
driver_variants_ALL$count <- 1
driver_variants_ALL <- unique(driver_variants_ALL)

# Create binary matrix for specific variants
driver_variants_ALL2 <- spread(driver_variants_ALL, variant, count)
driver_variants_ALL2 <- merge(SAMPLES, driver_variants_ALL2, by = "sample_ID", all = TRUE)
driver_variants_ALL2[is.na(driver_variants_ALL2)] <- 0

# Create "other KRAS" category (non-G12 variants)
driver_variants_ALL2$KRAS_other <- driver_variants_ALL2$KRAS_A146 + 
  driver_variants_ALL2$KRAS_G13 + 
  driver_variants_ALL2$KRAS_Q22 + 
  driver_variants_ALL2$KRAS_Q61

# Cap at maximum of 1 (binary indicator)
driver_variants_ALL2[driver_variants_ALL2$KRAS_other > 1, ]$KRAS_other <- 1

# Remove columns with all zeros
driver_variants_ALL2 <- driver_variants_ALL2[, colSums(driver_variants_ALL2 != 0) > 0]

# NOTE: Save driver_variants_ALL2 as data_driver-variants_ALL-organoids.csv



################################################################################
# 5. PROCESS COPY NUMBER VARIATIONS (CNV)
################################################################################

# Load CNV segment data
CN_seg <- as.data.frame(
  fread('/path_to_folder/CNV-segments_organoids.tsv')
)

# Load gene position information for CNV annotation
position_genes <- as.data.frame(
  fread('/path_to_folder/position_genes.csv')
)

# Filter for organoid models
CN_seg2 <- CN_seg[CN_seg$model == "organoid", ]

# Calculate chromosome arm-level copy numbers
CN_seg3 <- NULL
for (s in unique(CN_seg2$sample_ID)) {
  
  print(s)
  df_temp <- CN_seg2[CN_seg2$sample_ID == s, ]
  s2 <- unique(df_temp$sample)
  
  # Process each chromosome
  for (c in unique(df_temp$chromosome)) {
    
    print(c)
    df_temp2 <- df_temp[df_temp$chromosome == c, ]
    
    # Identify chromosome arm boundaries (p and q arms)
    START_p <- df_temp2[df_temp2$segmentStartSupport == "TELOMERE", ]$start
    END_p <- df_temp2[df_temp2$segmentEndSupport == "CENTROMERE", ]$end
    START_q <- df_temp2[df_temp2$segmentStartSupport == "CENTROMERE", ]$start
    END_q <- df_temp2[df_temp2$segmentEndSupport == "TELOMERE", ]$end
    
    # Assign segments to p or q arm
    df_temp2$arm <- "X"
    df_temp2[df_temp2$start >= START_q & df_temp2$end <= END_q, ]$arm <- "q"
    df_temp2[df_temp2$start >= START_p & df_temp2$end <= END_p, ]$arm <- "p"
    df_temp2$chromosome_arm <- paste(df_temp2$chromosome, df_temp2$arm, sep = "")
    
    # Calculate weighted average copy number for each arm
    for (a in unique(df_temp2$arm)) {
      
      df_temp3 <- df_temp2[df_temp2$arm == a, ]
      MAX <- max(df_temp3$end)
      df_temp3$size_MAX <- df_temp3$size_segment / MAX
      df_temp3$CN_size <- df_temp3$copyNumber * df_temp3$size_MAX
      
      # Calculate copy number adjusted for ploidy
      CN <- sum(df_temp3$CN_size)
      ploidy <- unique(df_temp3$ploidy)
      CNploidy <- round((2^log2(CN / ploidy)) * 2)
      
      CN_seg3 <- rbind(CN_seg3, data.frame(s, c, a, CNploidy))
    }
  }
}

colnames(CN_seg3) <- c("sample_ID", "chromosome", "arm", "copyNumber")

# Categorize copy number states
CN_seg3[CN_seg3$copyNumber >= 4, ]$copyNumber <- "Amp"
CN_seg3[CN_seg3$copyNumber <= 0, ]$copyNumber <- "Del"
CN_seg3[CN_seg3$copyNumber == 1, ]$copyNumber <- "Loss"
CN_seg3[CN_seg3$copyNumber == 2, ]$copyNumber <- "Neutral"
CN_seg3[CN_seg3$copyNumber == 3, ]$copyNumber <- "Soft gain"

CN_seg3 <- unique(CN_seg3)

# Process amplifications
CNV_amp_ALL <- mut_CNV_SV2[mut_CNV_SV2$Genomic_alterations == "Amplification", ]
CNV_amp_ALL <- merge(CNV_amp_ALL, position_genes)

# Create binary matrix for focal amplifications
CNV_amp_ALL2 <- unique(CNV_amp_ALL[, c("sample_ID", "gene")])
CNV_amp_ALL2$count <- 1
CNV_amp_ALL3 <- spread(CNV_amp_ALL2, gene, count)

# Create binary matrix for arm-level amplifications
CNV_amp_ALL_arm <- CN_seg3[CN_seg3$copyNumber == "Amp", ]
CNV_amp_ALL_arm$chromosome_arm <- paste(
  CNV_amp_ALL_arm$chromosome, 
  CNV_amp_ALL_arm$arm, 
  sep = ""
)
CNV_amp_ALL_arm$count <- 1
CNV_amp_ALL_arm$arm <- NULL
CNV_amp_ALL_arm$chromosome <- NULL
CNV_amp_ALL_arm$copyNumber <- NULL
CNV_amp_ALL_arm2 <- spread(CNV_amp_ALL_arm, chromosome_arm, count)

# Merge focal and arm-level amplifications
CNV_amp_ALL3 <- merge(CNV_amp_ALL3, CNV_amp_ALL_arm2, by = "sample_ID", all = TRUE)
CNV_amp_ALL3 <- merge(SAMPLES, CNV_amp_ALL3, by = "sample_ID", all = TRUE)
CNV_amp_ALL3[is.na(CNV_amp_ALL3)] <- 0

# Add amplification suffix to column names
colnames(CNV_amp_ALL3)[2:ncol(CNV_amp_ALL3)] <- paste(
  colnames(CNV_amp_ALL3)[2:ncol(CNV_amp_ALL3)], 
  "amp", 
  sep = "_"
)

# Process deletions (same logic as amplifications)
CNV_del_ALL <- mut_CNV_SV2[mut_CNV_SV2$Genomic_alterations == "Large deletion", ]
CNV_del_ALL <- merge(CNV_del_ALL, position_genes)

CNV_del_ALL2 <- unique(CNV_del_ALL[, c("sample_ID", "gene")])
CNV_del_ALL2$count <- 1
CNV_del_ALL3 <- spread(CNV_del_ALL2, gene, count)

CNV_del_ALL_arm <- CN_seg3[CN_seg3$copyNumber == "Del", ]
CNV_del_ALL_arm$chromosome_arm <- paste(
  CNV_del_ALL_arm$chromosome, 
  CNV_del_ALL_arm$arm, 
  sep = ""
)
CNV_del_ALL_arm$count <- 1
CNV_del_ALL_arm$arm <- NULL
CNV_del_ALL_arm$chromosome <- NULL
CNV_del_ALL_arm$copyNumber <- NULL
CNV_del_ALL_arm2 <- spread(CNV_del_ALL_arm, chromosome_arm, count)

CNV_del_ALL3 <- merge(CNV_del_ALL3, CNV_del_ALL_arm2, by = "sample_ID", all = TRUE)
CNV_del_ALL3 <- merge(SAMPLES, CNV_del_ALL3, by = "sample_ID", all = TRUE)
CNV_del_ALL3[is.na(CNV_del_ALL3)] <- 0

# Add deletion suffix to column names
colnames(CNV_del_ALL3)[2:ncol(CNV_del_ALL3)] <- paste(
  colnames(CNV_del_ALL3)[2:ncol(CNV_del_ALL3)], 
  "del", 
  sep = "_"
)

# Merge amplification and deletion data
CNV_ALL <- merge(CNV_amp_ALL3, CNV_del_ALL3)



################################################################################
# 6. PROCESS RNA-SEQ EXPRESSION DATA
################################################################################

# Load TPM (transcripts per million) expression data
RNASEQ_tpm <- as.data.frame(
  fread('/path_to_folder/supplementary_table_3.6.csv')
)

# Convert to long format
RNASEQ_tpm2 <- gather(RNASEQ_tpm, "sample_ID", "tpm", 3:ncol(RNASEQ_tpm))

# Log2 transform TPM values (adding 1 to avoid log(0))
RNASEQ_tpm2$log2TPM <- log2(RNASEQ_tpm2$tpm + 1)
RNASEQ_tpm2$tpm <- NULL

# Calculate mean expression across all samples
RNASEQ_summary_ALL <- RNASEQ_tpm2 %>% 
  group_by(gene, ENSEMBL_ID) %>% 
  summarise(MEAN_log2TPM = mean(log2TPM))

# Filter out lowly expressed genes (mean log2TPM > 0.1)
RNASEQ_summary_ALL <- RNASEQ_summary_ALL[RNASEQ_summary_ALL$MEAN_log2TPM > 0.1,]

# Keep only expressed genes
RNASEQ_ALL <- RNASEQ_tpm2[RNASEQ_tpm2$ENSEMBL_ID %in% RNASEQ_summary_ALL$ENSEMBL_ID,]
RNASEQ_ALL2 <- RNASEQ_ALL
RNASEQ_ALL2$gene_id <- NULL
RNASEQ_ALL2$cancer_type <- NULL

# Add expression suffix to gene names
RNASEQ_ALL2$gene <- paste(RNASEQ_ALL2$gene, "exp", sep = "_")

# Convert to wide format
RNASEQ_ALL2 <- spread(RNASEQ_ALL2, "gene", "log2TPM")

# NOTE: Save RNASEQ_ALL2 as data_RNAseq_ALL-organoids.csv



################################################################################
# 7. IDENTIFY LOSS-OF-FUNCTION (LOF) EVENTS
################################################################################

# Identify LoF driver genes from reference
LOF_genes <- cancer_driver_genes[cancer_driver_genes$method_of_action == "LoF", ]
ambiguous_genes <- cancer_driver_genes[cancer_driver_genes$method_of_action == "ambiguous", ]

# Extract biallelic LoF mutations from genomic data
drivers_LOF <- mut_CNV_SV2[
  (mut_CNV_SV2$gene %in% LOF_genes$gene & mut_CNV_SV2$biallelic_corrected == "Yes") | 
  (mut_CNV_SV2$gene %in% ambiguous_genes$gene & 
   mut_CNV_SV2$biallelic_corrected == "Yes" & 
   mut_CNV_SV2$Genomic_alterations %in% c("nonsense", "frameshift", "stop_lost", 
                                          "Large deletion", "ess_splice")),]

drivers_LOF <- unique(drivers_LOF[, c("sample_ID", "gene")])
drivers_LOF$count <- 1

# Identify LoF events from low/absent expression
RNASEQ_ALL_LOF <- RNASEQ_ALL[
  (RNASEQ_ALL$gene %in% LOF_genes$gene | 
   RNASEQ_ALL$gene %in% ambiguous_genes$gene) & 
   RNASEQ_ALL$log2TPM < 0.1,]

RNASEQ_ALL_LOF$log2TPM <- NULL
RNASEQ_ALL_LOF$ENSEMBL_ID <- NULL
RNASEQ_ALL_LOF$count <- 1

# Keep only genes present in RNA-seq data
drivers_ALL_LOF <- drivers_LOF[drivers_LOF$gene %in% RNASEQ_ALL$gene,]

# Combine mutation-based and expression-based LoF events
LOF_ALL <- rbind(drivers_ALL_LOF, RNASEQ_ALL_LOF)
LOF_ALL <- LOF_ALL %>% 
  dplyr::group_by(sample_ID, gene) %>% 
  dplyr::summarise(count = sum(count))
LOF_ALL$count <- 1

# Add LoF suffix to gene names
LOF_ALL$gene <- paste(LOF_ALL$gene, "LoF", sep = "_")

# Convert to wide format
LOF_ALL2 <- spread(LOF_ALL, gene, count)
LOF_ALL2 <- merge(SAMPLES, LOF_ALL2, by = "sample_ID", all.x = TRUE)
LOF_ALL2[is.na(LOF_ALL2)] <- 0

# NOTE: Save LOF_ALL2 as data_LOF_ALL-organoids.csv



################################################################################
# 8. IDENTIFY GAIN-OF-FUNCTION (GOF) EVENTS
################################################################################

# Identify activating/oncogenic driver genes
Act_genes <- cancer_driver_genes[cancer_driver_genes$method_of_action == "Act", ]

# Extract GoF mutations (missense, inframe, amplifications)
drivers_GOF <- mut_CNV_SV2[mut_CNV_SV2$gene %in% Act_genes$gene | (mut_CNV_SV2$gene %in% ambiguous_genes$gene & mut_CNV_SV2$Genomic_alterations %in% c("missense", "inframe", "Amplification")),]

# Exclude deletions/SVs for activating genes
drivers_GOF <- drivers_GOF[!(drivers_GOF$gene %in% cancer_driver_genes[cancer_driver_genes$method_of_action == "Act", ]$gene & drivers_GOF$Genomic_alterations %in% c("Large deletion", "SV")),]

drivers_GOF <- unique(drivers_GOF[, c("sample_ID", "gene")])
drivers_GOF$count <- 1

# Identify GoF events from high expression (top 10th percentile per sample)
RNASEQ_GOF <- RNASEQ_ALL[RNASEQ_ALL$gene %in% Act_genes$gene | RNASEQ_ALL$gene %in% ambiguous_genes$gene,]

RNASEQ_GOF2 <- NULL
for (s in unique(RNASEQ_GOF$sample_ID)) {
  
  print(s)
  df_temp <- RNASEQ_GOF[RNASEQ_GOF$sample_ID == s, ]
  
  # Calculate 90th percentile expression threshold per sample
  df_temp$expression_quantile_tpm <- df_temp$log2TPM > quantile(
    df_temp$log2TPM, 
    probs = 0.90, 
    na.rm = TRUE
  )
  
  RNASEQ_GOF2 <- rbind(RNASEQ_GOF2, df_temp)
}

# Filter for overexpressed genes (>90th percentile AND log2TPM > 5)
RNASEQ_GOF2 <- RNASEQ_GOF2[RNASEQ_GOF2$expression_quantile_tpm == TRUE & RNASEQ_GOF2$log2TPM > 5,]
RNASEQ_GOF2$log2TPM <- NULL
RNASEQ_GOF2$ENSEMBL_ID <- NULL
RNASEQ_GOF2$expression_quantile_tpm <- NULL
RNASEQ_GOF2$count <- 1

# Keep only genes present in RNA-seq data
drivers_GOF <- drivers_GOF[drivers_GOF$gene %in% RNASEQ_ALL$gene, ]

# Combine mutation-based and expression-based GoF events
GOF_ALL <- rbind(drivers_GOF, RNASEQ_GOF2)
GOF_ALL <- GOF_ALL %>% 
  dplyr::group_by(sample_ID, gene) %>% 
  dplyr::summarise(count = sum(count))
GOF_ALL$count <- 1

# Add GoF suffix to gene names
GOF_ALL$gene <- paste(GOF_ALL$gene, "GoF", sep = "_")

# Convert to wide format
GOF_ALL2 <- spread(GOF_ALL, gene, count)
GOF_ALL2 <- merge(SAMPLES, GOF_ALL2, all.x = TRUE)
GOF_ALL2[is.na(GOF_ALL2)] <- 0

# NOTE: Save GOF_ALL2 as data_GOF_ALL-organoids.csv



################################################################################
# 9. PROCESS CMS AND CRIS MOLECULAR SUBTYPES
################################################################################

# Load CMS (Consensus Molecular Subtypes) and CRIS classifications
subtypes <- as.data.frame(
  fread('/path_to_folder/supplementary_table_3.7.csv')
)

# Process CMS subtypes
CMS <- subtypes[, c("sample_ID", "CMS_prediction")]
CMS$value <- 1

CMS2 <- spread(CMS, "CMS_prediction", "value")
CMS2[is.na(CMS2)] <- 0
CMS2$`<NA>` <- NULL

# NOTE: Save CMS2 as data_CMS-subtypes_organoids.csv

# Process CRIS subtypes
CRIS <- subtypes[, c("sample_ID", "CRIS_prediction")]
CRIS$value <- 1

CRIS2 <- spread(CRIS, "CRIS_prediction", "value")
CRIS2[is.na(CRIS2)] <- 0
CRIS2$`<NA>` <- NULL

# NOTE: Save CRIS2 as data_CRIS-subtypes_organoids.csv



################################################################################
# 10. PROCESS MUTATIONAL SIGNATURES
################################################################################

# Load mutational signature data
signatures <- as.data.frame(
  fread('/path_to_folder/supplementary_table_3.5.csv')
)

# Filter for organoid samples
signatures2 <- merge(SAMPLES, signatures, by = "sample_ID")
signatures_ORG <- signatures2[signatures2$model == "organoid", ]

# Create binary indicator (signature >5% of total mutations)
signatures_ORG$value <- 0
signatures_ORG[signatures_ORG$percentage_total > 5, ]$value <- 1

# Remove unnecessary columns
signatures_ORG2 <- signatures_ORG
signatures_ORG2$model <- NULL
signatures_ORG2$percentage_total <- NULL
signatures_ORG2$number_mutations <- NULL
signatures_ORG2$annotation <- NULL

# Convert to wide format
signatures_ORG2 <- spread(signatures_ORG2, "signature", "value")
signatures_ORG2[is.na(signatures_ORG2)] <- 0

# NOTE: Save signatures_ORG2 as data_mutational-signatures_ALL-organoids.csv



################################################################################
# 11. PROCESS CLINICAL DATA
################################################################################

# Load clinical metadata
clinical <- as.data.frame(
  fread('/path_to_folder/supplementary_table_2.2.csv')
)

## 11a. Process tumor stages
stages <- clinical[, c("sample_ID", "pTNM_pathological_stage")]
colnames(stages)[2] <- "Stage2"

# Standardize stage nomenclature
stages$Stage <- stages$Stage2
stages$Stage <- gsub("Stage IA", "stage_I", stages$Stage)
stages$Stage <- gsub("Stage IB", "stage_I", stages$Stage)
stages$Stage <- gsub("Stage IIA", "stage_II", stages$Stage)
stages$Stage <- gsub("Stage IIB", "stage_II", stages$Stage)
stages$Stage <- gsub("Stage IIC", "stage_II", stages$Stage)
stages$Stage <- gsub("Stage IIIA2", "stage_III", stages$Stage)
stages$Stage <- gsub("Stage IIIA", "stage_III", stages$Stage)
stages$Stage <- gsub("Stage IIIB", "stage_III", stages$Stage)
stages$Stage <- gsub("Stage IIIC", "stage_III", stages$Stage)
stages$Stage <- gsub("Stage IVA", "stage_IV", stages$Stage)
stages$Stage <- gsub("Stage IVB", "stage_IV", stages$Stage)
stages$Stage <- gsub("Stage $", "NC", stages$Stage)
stages$Stage <- gsub("Stage 0", "stage_I", stages$Stage)
stages$Stage <- gsub("^$", "NC", stages$Stage)
stages$Stage <- gsub("Stage IV", "stage_IV", stages$Stage)
stages$Stage <- gsub("Stage III", "stage_III", stages$Stage)
stages$Stage <- gsub("Stage II", "stage_II", stages$Stage)
stages$Stage <- gsub("Stage I", "stage_I", stages$Stage)

stages$Stage2 <- NULL
stages$value <- 1

# Create binary stage indicators
stages2 <- spread(stages, Stage, value)
stages2[is.na(stages2)] <- 0

# Set unknown stages to NA (not 0)
stages2[stages2$NC == 1, ]$stage_I <- NA
stages2[stages2$NC == 1, ]$stage_II <- NA
stages2[stages2$NC == 1, ]$stage_III <- NA
stages2[stages2$NC == 1, ]$stage_IV <- NA
stages2$NC <- NULL

# NOTE: Save stages2 as data_binary-stages_ALL-organoids.csv

# Create early vs advanced stage indicators
stages3 <- stages2
stages3$advanced_stage <- stages3$stage_III + stages3$stage_IV
stages3$early_stage <- stages3$stage_I + stages3$stage_II
stages3 <- stages3[, c("sample_ID", "advanced_stage", "early_stage")]

# NOTE: Save stages3 as data_binary-stages2_ALL-organoids.csv

## 11b. Process patient age
age <- clinical[, c("sample_ID", "age_at_attendance")]
age <- age[complete.cases(age$age_at_attendance), ]
colnames(age)[2] <- "age"

# Create age categories
age[age$age >= 80, ]$age <- ">= 80"
age[age$age >= 70 & age$age < 80, ]$age <- "70-79"
age[age$age >= 60 & age$age < 70, ]$age <- "60-69"
age[age$age >= 50 & age$age < 60, ]$age <- "50-59"
age[age$age < 50 & age$age > 30, ]$age <- "< 50"
age$value <- 1

# Create binary age matrix
age2 <- spread(age, age, value)
age2[is.na(age2)] <- 0

# Create old vs young categories
age2$old <- age2$`50-59` + age2$`60-69` + age2$`70-79` + age2$`>= 80`
age2$young <- age2$`< 50`
age2 <- age2[, c("sample_ID", "old", "young")]

# NOTE: Save age2 as data_binary-age_ALL-organoids.csv

## 11c. Process patient sex
sex <- clinical[, c("sample_ID", "sex")]
sex$female <- 0
sex[sex$sex == "Female", ]$female <- 1
sex$sex <- NULL

# NOTE: Save sex as data_binary-sex_ALL-organoids.csv

## 11d. Process colorectal tumor anatomical location
side_COLO <- clinical[clinical$primary_tumour_type == "Colorectal", ]
side_COLO <- side_COLO[, c("sample_ID", "anatomic_site_primary")]

# Classify by anatomical location based on ICD-10 codes
side_COLO$anatomical_site <- "Unknown"
side_COLO[grep("C18.[0-4]", side_COLO$anatomic_site_primary), ]$anatomical_site <- "Proximal"
side_COLO[grep("C18.[5-8]", side_COLO$anatomic_site_primary), ]$anatomical_site <- "Distal"
side_COLO[grep("C20", side_COLO$anatomic_site_primary), ]$anatomical_site <- "Rectum"
side_COLO[grep("C19", side_COLO$anatomic_site_primary), ]$anatomical_site <- "Rectum"
side_COLO$anatomic_site_primary <- NULL
side_COLO$value <- 1

# Create binary anatomical location matrix
side_COLO2 <- spread(side_COLO, "anatomical_site", "value")
side_COLO2[is.na(side_COLO2)] <- 0

# Set unknown locations to NA
side_COLO2[side_COLO2$Unknown == 1, ]$Distal <- NA
side_COLO2[side_COLO2$Unknown == 1, ]$Proximal <- NA
side_COLO2[side_COLO2$Unknown == 1, ]$Rectum <- NA
side_COLO2$Unknown <- NULL

# NOTE: Save side_COLO2 as data_binary-sides_COLO-organoids.csv

## 11e. Process metastasis status (colorectal only)
metastasis_COLO <- clinical[clinical$primary_tumour_type == "Colorectal", ]
metastasis_COLO <- metastasis_COLO[, c("sample_ID", "metastasis")]
metastasis_COLO[metastasis_COLO$metastasis == "No", ]$metastasis <- 0
metastasis_COLO[metastasis_COLO$metastasis == "Yes", ]$metastasis <- 1

# NOTE: Save metastasis_COLO as data_binary-metastasis_COLO-organoids.csv

## 11f. Process Barrett's esophagus status (esophageal only)
Barrett_OESO <- clinical[clinical$primary_tumour_type == "Oesophageal", ]
Barrett_OESO <- Barrett_OESO[, c("sample_ID", "diagnosis_barrett's")]
colnames(Barrett_OESO)[2] <- "Barrett"
Barrett_OESO[Barrett_OESO$Barrett == "Yes", ]$Barrett <- 1
Barrett_OESO[Barrett_OESO$Barrett == "No", ]$Barrett <- 0
Barrett_OESO[Barrett_OESO$Barrett %in% c("Unknown", ""), ]$Barrett <- NA

# NOTE: Save Barrett_OESO as data_binary-barrett_OESO-organoids.csv



################################################################################
# 12. PROCESS GENOMIC COVARIATES
################################################################################

# Load additional WGS metrics
other_WGS <- as.data.frame(
  fread('/path_to_folder/supplementary_table_3.1.csv')
)
other_WGS_ORG <- other_WGS[other_WGS$model == "organoid", ]

# Extract ploidy
ploidy <- unique(other_WGS_ORG[, c("sample_ID", "ploidy")])

# Extract microsatellite instability (MSI) status
MSI <- other_WGS_ORG[, c("sample_ID", "msStatus")]

# Combine covariates
covariates_ALL <- merge(ploidy, MSI, by = "sample_ID")

# NOTE: Save covariates_ALL as data_covariates_ALL-organoids.csv

# Add MSI status to CNV data
CNV_ALL <- merge(MSI, CNV_ALL, by = "sample_ID")

# NOTE: Save CNV_ALL as data_CNV-genes_ALL-organoids.csv
