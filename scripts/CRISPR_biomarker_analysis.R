# Load required libraries
library(data.table)
library(tidyverse)
library(digest)
library(scales)
library(dplyr)
library(lmtest)
library(parallel)
library(doParallel)
library(foreach)
library(ggrastr)

################################################################################
# 1. LOAD CLINICAL AND MOLECULAR FEATURES
################################################################################

stages1 <- as.data.frame(fread('/path_to_folder/data_binary-stages_ALL-organoids.csv'))
stages2 <- as.data.frame(fread('/path_to_folder/data_binary-stages2_ALL-organoids.csv'))
age <- as.data.frame(fread('/path_to_folder/data_binary-age_ALL-organoids.csv'))
gender <- as.data.frame(fread('/path_to_folder/data_binary-sex_ALL-organoids.csv'))
localization_COADREAD <- as.data.frame(fread('/path_to_folder/data_binary-sides_COLO-organoids.csv'))
metastasis_COADREAD <- as.data.frame(fread('/path_to_folder/data_binary-metastasis_COLO-organoids.csv'))
barrett_ESCA <- as.data.frame(fread('/path_to_folder/data_binary-barrett_OESO-organoids.csv'))
CMS_COADREAD <- as.data.frame(fread('/path_to_folder/data_CMS-subtypes_organoids.csv'))
CRIS_COADREAD <- as.data.frame(fread('/path_to_folder/data_CRIS-subtypes_organoids.csv'))

clinical <- merge(stages1, stages2, by = "sample_ID")
clinical <- merge(clinical, age, by = "sample_ID", all.x = TRUE)
clinical <- merge(clinical, gender, by = "sample_ID", all.x = TRUE)
clinical <- merge(clinical, localization_COADREAD, by = "sample_ID")
clinical <- merge(clinical, metastasis_COADREAD, by = "sample_ID", all.x = TRUE)
clinical <- merge(clinical, barrett_ESCA, by = "sample_ID", all.x = TRUE)
clinical <- merge(clinical, CMS_COADREAD, by = "sample_ID", all.x = TRUE)
clinical <- merge(clinical, CRIS_COADREAD, by = "sample_ID", all.x = TRUE)

clinical$stage_II <- NULL
clinical$stage_III <- NULL
clinical$old <- NULL
clinical$advanced_stage <- NULL


################################################################################
# 2. LOAD GENOMIC FEATURES
################################################################################

gene_mutations <- as.data.frame(fread('/path_to_folder/data_driver-mutations_ALL-organoids.csv'))
allele_mutations <- as.data.frame(fread('/path_to_folder/data_driver-variants_ALL-organoids.csv'))
mutations <- merge(gene_mutations, allele_mutations, by = "sample_ID", all = TRUE)
SCNA <- as.data.frame(fread('/path_to_folder/data_CNV-genes_ALL-organoids.csv'))
RNASEQ <- as.data.frame(fread('/path_to_folder/data_RNAseq_ALL-organoids.csv'))
signatures <- as.data.frame(fread('/path_to_folder/data_mutational-signatures_ALL-organoids.csv'))
LOF <- as.data.frame(fread('/path_to_folder/data_LOF_ALL-organoids.csv'))
GOF <- as.data.frame(fread('/path_to_folder/data_GOF_ALL-organoids.csv'))
LoF_GoF <- merge(LOF, GOF, by = "sample_ID", all = TRUE)


################################################################################
# 3. LOAD COVARIATES
################################################################################

covariates <- as.data.frame(fread('/path_to_folder/data_covariates_ALL-organoids.csv'))
MSI <- covariates[, c("sample_ID", "msStatus")]


################################################################################
# 4. MERGE ALL BINARY FEATURES
################################################################################

MERGE_binary <- merge(mutations, SCNA, by = "sample_ID", all = TRUE)
MERGE_binary <- merge(MERGE_binary, LoF_GoF, by = "sample_ID", all = TRUE)
MERGE_binary <- merge(MERGE_binary, signatures, by = "sample_ID", all = TRUE)
MERGE_binary <- merge(MERGE_binary, clinical, by = "sample_ID", all = TRUE)
MERGE_binary$msStatus <- NULL


################################################################################
# 5. LOAD CRISPR DEPENDENCY DATA (RESPONSES)
################################################################################

LFC <- as.data.frame(fread('/path_to_folder/supplementary_table_6.csv'))

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

LFC_filtered <- LFC_filtered[, c("sample_ID", "gene", "LFC")]


################################################################################
# 6. LOAD DIFFERENTIAL DEPENDENCY ANALYSIS RESULTS
################################################################################

differential_dep_gastrointestinal <- as.data.frame(fread("/path_to_folder/supplementary_table_7.1.csv"))
differential_dep_gastrointestinal$ct <- "Gastrointestinal"

differential_dep_COADREAD <- as.data.frame(fread("/path_to_folder/supplementary_table_7.2.csv"))
differential_dep_COADREAD$ct <- "Colorectal"

differential_dep_ESCA <- as.data.frame(fread("/path_to_folder/supplementary_table_7.3.csv"))
differential_dep_ESCA$ct <- "Oesophageal"

differential_dep <- rbind(differential_dep_gastrointestinal, differential_dep_COADREAD, differential_dep_ESCA)


################################################################################
# 7. LOAD CANCER TYPE ANNOTATIONS
################################################################################

cancer_type <- as.data.frame(fread('/path_to_folder/supplementary_table_2.2.csv'))
cancer_type <- cancer_type[, c(1, 4)]

cancer_type2 <- cancer_type[cancer_type$primary_tumour_type %in% c("Colorectal", "Oesophageal", "Gastric", "Pancreatic"), ]
cancer_type2$primary_tumour_type <- "Gastrointestinal"
cancer_type <- rbind(cancer_type, cancer_type2)


################################################################################
# 8. LINEAR REGRESSION BIOMARKER ANALYSIS
################################################################################

TYPES <- c("Colorectal", "Oesophageal", "Gastrointestinal")

# Detect available cores and register parallel backend
n_cores <- max(1, detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
message(paste("Running with", n_cores, "parallel cores"))

linear_model_results_RESPONSES <- NULL

for (ct in TYPES) {
  
  message(paste0("[", Sys.time(), "] Processing cancer type: ", ct))
  
  ##########################################################
  ## 8a. Filter data for current cancer type
  ##########################################################
  
  LFC_filtered_ct <- LFC_filtered[
    LFC_filtered$sample_ID %in% cancer_type[cancer_type$primary_tumour_type == ct, ]$sample_ID,
  ]
  LFC_filtered_ct <- spread(LFC_filtered_ct, "gene", "LFC")
  
  MERGE_binary_ct     <- MERGE_binary[MERGE_binary$sample_ID %in% LFC_filtered_ct$sample_ID, ]
  RNASEQ_ct           <- RNASEQ[RNASEQ$sample_ID %in% LFC_filtered_ct$sample_ID, ]
  differential_dep_ct <- differential_dep[differential_dep$ct == ct, ]
  covariates_ct       <- covariates[covariates$sample_ID %in% LFC_filtered_ct$sample_ID, ]
  
  if (ct == "Oesophageal") {
    covariates_ct$msStatus <- NULL
    covariates_ct$tissue   <- NULL
  }
  if (ct == "Colorectal") {
    covariates_ct$tissue <- NULL
  }
  
  ##########################################################
  ## 8b. Filter genes: Keep only differentially dependent genes
  ##########################################################
  
  LFC_filtered_ct <- LFC_filtered_ct[
    , colnames(LFC_filtered_ct) %in% c(
      "sample_ID",
      differential_dep_ct[
        differential_dep_ct$N_not_depleted > 2 & differential_dep_ct$N_depleted > 2,
      ]$g
    )
  ]
  
  ##########################################################
  ## 8c. FILTER BINARY FEATURES
  ##########################################################
  
  message(paste0("[", Sys.time(), "] Filtering binary features..."))
  
  feature_cols <- colnames(MERGE_binary_ct)[2:ncol(MERGE_binary_ct)]
  
  # Filter 1: Remove feature absent or present in less than 3 organoids
  keep_feature <- vapply(feature_cols, function(c) {
    v <- MERGE_binary_ct[[c]]
    v <- v[!is.na(v)]
    s <- sum(v)
    s >= 3 & s <= (length(v) - 3) & s != length(v)
  }, logical(1))
  
  binary_features      <- feature_cols[keep_feature]
  MERGE_binary_ct_temp <- MERGE_binary_ct[, c("sample_ID", binary_features), drop = FALSE]
  
  # Filter 2: Remove features identical to covariates
  remove_columns <- c()
  for (c3 in colnames(MERGE_binary_ct_temp)[2:ncol(MERGE_binary_ct_temp)]) {
    df_temp1 <- MERGE_binary_ct_temp[, c("sample_ID", c3)]
    for (c4 in colnames(covariates_ct)[2:ncol(covariates_ct)]) {
      df_temp2 <- covariates_ct[, c("sample_ID", c4)]
      df_temp3 <- merge(df_temp1, df_temp2, by = "sample_ID")
      df_temp3 <- df_temp3[complete.cases(df_temp3), ]
      if (identical(df_temp3[[2]], df_temp3[[3]])) {
        remove_columns <- c(remove_columns, c3)
      }
    }
  }
  MERGE_binary_ct_temp <- MERGE_binary_ct_temp[
    , !colnames(MERGE_binary_ct_temp) %in% remove_columns, drop = FALSE
  ]
  
  # Filter 3: Remove duplicate features
  nondups <- MERGE_binary_ct_temp[!duplicated(lapply(MERGE_binary_ct_temp, digest))]
  dups    <- MERGE_binary_ct_temp[duplicated(lapply(MERGE_binary_ct_temp, digest))]
  
  # Initialize per cancer type to avoid carry-over across iterations
  duplicated_names <- NULL
  
  if (ncol(dups) > 0) {
    for (i in 1:ncol(nondups)) {
      for (j in 1:ncol(dups)) {
        if (!FALSE %in% (nondups[, i] == dups[, j])) {
          duplicated_names <- rbind(
            duplicated_names,
            data.frame(names(nondups[i]), names(dups[j]))
          )
        }
      }
    }
    colnames(duplicated_names) <- c("Feature", "Feature_duplicated")
    MERGE_binary_ct_temp <- MERGE_binary_ct_temp[
      , !colnames(MERGE_binary_ct_temp) %in% colnames(dups), drop = FALSE
    ]
  }
  
  # Filter 4: Chi-square test for correlated feature pairs
  CHI_SQ  <- NULL
  couples <- c()
  feat_cols_temp <- colnames(MERGE_binary_ct_temp)[2:ncol(MERGE_binary_ct_temp)]
  
  for (c1 in feat_cols_temp) {
    for (c2 in feat_cols_temp) {
      couple_paste  <- ifelse(c1 < c2, paste(c1, c2, sep = "-"), paste(c2, c1, sep = "-"))
      couple_paste2 <- unlist(strsplit(couple_paste, "-"))
      if (c1 != c2 & !(couple_paste %in% couples)) {
        df_temp4 <- MERGE_binary_ct_temp[, c(c1, c2)]
        df_temp4 <- df_temp4[complete.cases(df_temp4), ]
        df_temp4 <- df_temp4 %>% mutate_all(as.numeric)
        if ((sum(df_temp4[, c1]) < nrow(df_temp4) & sum(df_temp4[, c1]) > 0) |
            (sum(df_temp4[, c2]) < nrow(df_temp4) & sum(df_temp4[, c2]) > 0)) {
          chi_sq_table  <- table(df_temp4[, c1], df_temp4[, c2])
          chi_qr_test   <- chisq.test(chi_sq_table)
          chi_qr_result <- data.frame(chi_qr_test$p.value, couple_paste2[1], couple_paste2[2])
          CHI_SQ  <- rbind(CHI_SQ, chi_qr_result)
        } else {
          chi_qr_result_0 <- data.frame(0, couple_paste2[1], couple_paste2[2])
          colnames(chi_qr_result_0)[1] <- "chi_qr_test.p.value"
          CHI_SQ <- rbind(CHI_SQ, chi_qr_result_0)
        }
        couples <- c(couples, paste(c1, c2, sep = "-"), paste(c2, c1, sep = "-"))
      }
    }
  }
  
  CHI_SQ$p.value.adj <- p.adjust(CHI_SQ$chi_qr_test.p.value, method = "BH")
  CHI_SQ$group       <- paste(CHI_SQ$couple_paste2.1., CHI_SQ$couple_paste2.2., sep = "-")
  SIG_group           <- CHI_SQ[CHI_SQ$p.value.adj < 0.5, ]$group
  
  # Compute feature correlation flags
  binary_feature_names <- colnames(MERGE_binary_ct_temp)[2:ncol(MERGE_binary_ct_temp)]
  
  feature_correlation_flags <- data.frame(
    Feature          = binary_feature_names,
    similar_variable = NA_character_,
    stringsAsFactors = FALSE
  )
  
  for (feature1 in sort(binary_feature_names)) {
    other_features <- sort(binary_feature_names[binary_feature_names != feature1])
    for (feature2 in other_features) {
      B <- ifelse(
        feature1 < feature2,
        paste(feature1, feature2, sep = "-"),
        paste(feature2, feature1, sep = "-")
      )
      if (B %in% SIG_group) {
        current <- feature_correlation_flags[feature_correlation_flags$Feature == feature1, "similar_variable"]
        feature_correlation_flags[feature_correlation_flags$Feature == feature1, "similar_variable"] <-
          ifelse(is.na(current), feature2, paste(current, feature2, sep = "/"))
      }
    }
  }
  
  ##########################################################
  ## 8d. FILTER EXPRESSION FEATURES
  ##########################################################
  
  message(paste0("[", Sys.time(), "] Filtering expression features..."))
  
  no_expression <- as.data.frame(colSums(RNASEQ_ct < 1))
  no_expression$gene <- rownames(no_expression)
  no_expression <- no_expression[no_expression[[1]] < nrow(RNASEQ_ct) * 0.95, ]
  RNASEQ_ct_temp <- RNASEQ_ct[, colnames(RNASEQ_ct) %in% no_expression$gene | colnames(RNASEQ_ct) == "sample_ID"]
  
  expr_cols <- colnames(RNASEQ_ct_temp)[2:ncol(RNASEQ_ct_temp)]
  
  keep_expr <- vapply(expr_cols, function(col) {
    g              <- RNASEQ_ct_temp[[col]]
    not_expressed  <- sum(g < 1)
    expressed      <- sum(g >= 1 & g <= 5)
    over_expressed <- sum(g > 5)
    (not_expressed >= 3 & expressed >= 3) |
      (not_expressed >= 3 & over_expressed >= 3) |
      (expressed >= 3 & over_expressed >= 3)
  }, logical(1))
  
  categories     <- expr_cols[keep_expr]
  RNASEQ_ct_temp <- RNASEQ_ct_temp[, c("sample_ID", categories), drop = FALSE]
  
  SD_values_vec <- vapply(
    colnames(RNASEQ_ct_temp)[2:ncol(RNASEQ_ct_temp)],
    function(col) sd(RNASEQ_ct_temp[[col]]),
    numeric(1)
  )
  SD_values <- data.frame(
    col = colnames(RNASEQ_ct_temp)[2:ncol(RNASEQ_ct_temp)],
    SD  = SD_values_vec
  )
  SD_values$Z.SCORE <- (SD_values$SD - mean(SD_values$SD)) / sd(SD_values$SD)
  
  MERGE_numeric_ct_temp <- RNASEQ_ct_temp[
    , colnames(RNASEQ_ct_temp) %in% SD_values[SD_values$Z.SCORE >= 2, ]$col |
      colnames(RNASEQ_ct_temp) %in% paste(colnames(LFC_filtered_ct), "exp", sep = "_") |
      colnames(RNASEQ_ct_temp) == "sample_ID",
    drop = FALSE
  ]
  
  ##########################################################
  ## 8e. SCALE FEATURES TO [0,1]
  ##########################################################
  
  DATA_scale <- merge(MERGE_binary_ct_temp, MERGE_numeric_ct_temp, by = "sample_ID", all = TRUE)
  DATA_scale[2:ncol(DATA_scale)] <- data.frame(
    lapply(DATA_scale[2:ncol(DATA_scale)], function(x) rescale(x, to = c(0, 1)))
  )
  
  ##########################################################
  ## 8f. CLASSIFY FEATURE TYPES
  ##########################################################
  
  message(paste0("[", Sys.time(), "] Classifying feature types..."))
  
  FEATURES_TYPE <- data.frame(
    FEATURES = as.character(colnames(DATA_scale)[2:ncol(DATA_scale)]),
    type     = "Mut",
    stringsAsFactors = FALSE
  )
  
  assign_type <- function(df, source_df, type_label) {
    if (!is.null(source_df)) {
      src_cols <- colnames(source_df)[2:ncol(source_df)]
      df[df$FEATURES %in% src_cols, "type"] <- type_label
    }
    df
  }
  
  FEATURES_TYPE <- assign_type(FEATURES_TYPE, allele_mutations,      "Var")
  FEATURES_TYPE <- assign_type(FEATURES_TYPE, SCNA,                  "CN")
  FEATURES_TYPE <- assign_type(FEATURES_TYPE, LoF_GoF,               "LoF_GoF")
  FEATURES_TYPE <- assign_type(FEATURES_TYPE, signatures,            "Sig")
  FEATURES_TYPE <- assign_type(FEATURES_TYPE, clinical,              "Clin")
  FEATURES_TYPE <- assign_type(FEATURES_TYPE, MERGE_numeric_ct_temp, "Expr")
  
  RESPONSES <- colnames(LFC_filtered_ct)[2:ncol(LFC_filtered_ct)]
  
  ##########################################################
  ## 8g. PRE-BUILD FULL DATA MATRIX 
  ##########################################################
  
  message(paste0("[", Sys.time(), "] Pre-building full data matrix..."))
  
  # Determine whether MSI covariate has variability
  if ("msStatus" %in% colnames(covariates_ct)) {
    msi_sum <- sum(covariates_ct$msStatus, na.rm = TRUE)
    if (msi_sum == 0 | msi_sum == nrow(covariates_ct)) {
      covariates_ct$msStatus <- NULL
    }
  }
  
  # Single join: LFC (wide) + all features + covariates, aligned by sample_ID
  FULL_DATA <- merge(LFC_filtered_ct, DATA_scale,    by = "sample_ID", all = FALSE)
  FULL_DATA <- merge(FULL_DATA,       covariates_ct, by = "sample_ID", all.x = TRUE)
  
  # Record which columns are responses, features, and covariates
  response_cols  <- colnames(LFC_filtered_ct)[2:ncol(LFC_filtered_ct)]
  feature_cols_all <- colnames(DATA_scale)[2:ncol(DATA_scale)]
  covariate_cols <- colnames(covariates_ct)[2:ncol(covariates_ct)]
  
  # Set of binary feature names for effect size branching
  binary_col_names <- colnames(MERGE_binary_ct_temp)[2:ncol(MERGE_binary_ct_temp)]
  
  ##########################################################
  ## 8h. LINEAR REGRESSION: MAIN ANALYSIS LOOP
  ## Parallelized over genes
  ##########################################################
  
  # Progress tracking via log file
  progress_file <- "/nfs/casm/team215mg/users/Carmen/projects/REGRESSIONS/output/script_crispr/o.GI.log"
  message(paste("Progress log:", progress_file))
  n_genes <- length(response_cols)
  # Reset log file for each cancer type
  cat(paste0("=== Cancer type: ", ct, " | Total genes: ", n_genes,
             " | Started: ", Sys.time(), " ===\n"),
      file = progress_file, append = FALSE)
  
  message(paste0("[", Sys.time(), "] Running linear regression for ", ct,
                 " (", n_genes, " genes across ", n_cores, " cores)..."))
  
   # Progress is tracked from the main process via a custom .combine function
  # This avoids nfs concurrent write issues with 63 workers
  n_done_main <- 0L
  combine_with_progress <- function(a, b) {
    n_done_main <<- n_done_main + 1L
    pct <- round(n_done_main / n_genes * 100)
    msg <- paste0("[", Sys.time(), "] Done: ", b$NAME[1],
                  " | Progress: ", n_done_main, "/", n_genes,
                  " (", pct, "%)")
    message(msg)
    cat(msg, "\n", file = progress_file, append = TRUE)
    if (is.null(a)) b else dplyr::bind_rows(a, b)
  }
  
  results_per_gene <- foreach(
    response       = response_cols,
    .combine       = rbind,
    .packages      = c("lmtest"),
    .export        = c(
      "FULL_DATA", "feature_cols_all", "covariate_cols",
      "FEATURES_TYPE", "binary_col_names",
      "feature_correlation_flags", "progress_file", "n_genes"
    )
  ) %dopar% {
    
    # Extract response vector for this gene (drop NAs)
    y_all   <- FULL_DATA[[response]]
    ok_rows <- !is.na(y_all)
    y       <- y_all[ok_rows]
    
    # Covariate matrix
    cov_df  <- FULL_DATA[ok_rows, covariate_cols, drop = FALSE]
    COV_mat <- model.matrix(~ . - 1, data = cov_df)
    
    # Pre-build null model design matrix (intercept + covariates)
    X0 <- cbind(1, COV_mat)
    
    # Store results in a pre-allocated list
    results_list <- vector("list", length(feature_cols_all))
    result_idx   <- 0L
    
    for (feature in feature_cols_all) {
      
      feat_vec <- FULL_DATA[[feature]][ok_rows]
      feat_na  <- is.na(feat_vec)
      
      # Drop rows where feature is NA
      if (any(feat_na)) {
        y_f    <- y[!feat_na]
        X0_f   <- X0[!feat_na, , drop = FALSE]
        feat_f <- feat_vec[!feat_na]
      } else {
        y_f    <- y
        X0_f   <- X0
        feat_f <- feat_vec
      }
      
      # Skip if too few observations
      if (length(y_f) < (ncol(X0_f) + 2)) next
      
      # Full model matrix (intercept + covariates + feature)
      X1_f <- cbind(X0_f, feat_f)
      
      # Linear regressions
      fit0 <- tryCatch(lm.fit(X0_f, y_f), error = function(e) NULL)
      fit1 <- tryCatch(lm.fit(X1_f, y_f), error = function(e) NULL)
      
      if (is.null(fit0) || is.null(fit1)) next
      
      # Extract residuals and degrees of freedom
      rss0 <- sum(fit0$residuals^2)
      rss1 <- sum(fit1$residuals^2)
      n    <- length(y_f)
      p0   <- ncol(X0_f)
      p1   <- ncol(X1_f)
      
      # LRT statistic: 2 * (logLik1 - logLik0) under normality assumption
      logLik0  <- -n / 2 * log(rss0 / n)
      logLik1  <- -n / 2 * log(rss1 / n)
      LRT_stat <- 2 * (logLik1 - logLik0)
      P_val_LRT <- pchisq(LRT_stat, df = p1 - p0, lower.tail = FALSE)
      
      # Feature coefficient and p-value from full model
      # lm.fit() stores coefficients but not std errors -- compute manually
      b_coef <- NA_real_
      P_val_b <- NA_real_
      R_sq    <- NA_real_
      
      tryCatch({
        coef_full <- fit1$coefficients
        b_coef    <- coef_full[length(coef_full)]  # last coef = feature
        
        # Residual std error and coefficient std error
        sigma2   <- rss1 / (n - p1)
        XtX_inv  <- tryCatch(solve(crossprod(X1_f)), error = function(e) NULL)
        if (!is.null(XtX_inv)) {
          se_feat  <- sqrt(sigma2 * XtX_inv[p1, p1])
          t_val    <- b_coef / se_feat
          P_val_b  <- 2 * pt(abs(t_val), df = n - p1, lower.tail = FALSE)
        }
        
        # R-squared of full model
        sst    <- sum((y_f - mean(y_f))^2)
        R_sq   <- 1 - rss1 / sst
      }, error = function(e) NULL)
      
      type <- FEATURES_TYPE[FEATURES_TYPE$FEATURES == feature, "type"]
      
      ###############################################################
      # EFFECT SIZE STATISTICS
      ###############################################################
      
      if (feature %in% binary_col_names) {
        
        g0 <- y_f[feat_f == 0]
        g1 <- y_f[feat_f == 1]
        
        # Safety check: need at least 1 sample in each group
        if (length(g0) < 1 || length(g1) < 1) next
        
        MEAN_g1   <- mean(g0);   MEDIAN_g1 <- median(g0)
        MIN_g1    <- min(g0);    MAX_g1    <- max(g0);    SD_g1 <- sd(g0)
        MEAN_g2   <- mean(g1);   MEDIAN_g2 <- median(g1)
        MIN_g2    <- min(g1);    MAX_g2    <- max(g1);    SD_g2 <- sd(g1)
        
      } else {
        
        q25 <- quantile(feat_f, 0.25)
        q75 <- quantile(feat_f, 0.75)
        lower_idx <- feat_f <= q25
        upper_idx <- feat_f >= q75
        
        # Safety check: need samples in both quartiles
        if (sum(lower_idx) < 1 || sum(upper_idx) < 1) next
        
        g_lower <- y_f[lower_idx]
        g_upper <- y_f[upper_idx]
        
        MEAN_g1   <- mean(g_lower);   MEDIAN_g1 <- median(g_lower)
        MIN_g1    <- min(g_lower);    MAX_g1    <- max(g_lower);    SD_g1 <- sd(g_lower)
        MEAN_g2   <- mean(g_upper);   MEDIAN_g2 <- median(g_upper)
        MIN_g2    <- min(g_upper);    MAX_g2    <- max(g_upper);    SD_g2 <- sd(g_upper)
      }
      
      CORR    <- tryCatch(cor(y_f, feat_f, use = "complete.obs"), error = function(e) NA_real_)
      cohens  <- if (!is.na(CORR) && CORR^2 < 1) (CORR^2) / (1 - CORR^2) else NA_real_
      
      result_idx <- result_idx + 1L
      results_list[[result_idx]] <- data.frame(
        NAME             = response,
        Feature          = feature,
        Feature_type     = type,
        b_coef           = b_coef,
        P_val_b          = P_val_b,
        R_squared_LR     = R_sq,
        P_val_LRT        = P_val_LRT,
        MEAN_group1      = MEAN_g1,
        MEDIAN_group1    = MEDIAN_g1,
        MIN_group1       = MIN_g1,
        MAX_group1       = MAX_g1,
        SD_group1        = SD_g1,
        MEAN_group2      = MEAN_g2,
        MEDIAN_group2    = MEDIAN_g2,
        MIN_group2       = MIN_g2,
        MAX_group2       = MAX_g2,
        SD_group2        = SD_g2,
        MEAN_difference  = MEAN_g2 - MEAN_g1,
        MEDIAN_difference = MEDIAN_g2 - MEDIAN_g1,
        cohens_f2        = cohens,
        stringsAsFactors = FALSE
      )
    }
    
    # Combine results for this gene using do.call
    DF_TYPE <- do.call(rbind, results_list[seq_len(result_idx)])
    
    if (is.null(DF_TYPE) || nrow(DF_TYPE) == 0) return(NULL)
    
    # Apply FDR correction per gene and per feature type
    DF_TYPE$P_val_LRT_BH <- NA_real_
    for (ft in unique(DF_TYPE$Feature_type)) {
      idx <- DF_TYPE$Feature_type == ft
      DF_TYPE$P_val_LRT_BH[idx] <- p.adjust(DF_TYPE$P_val_LRT[idx], method = "BH")
    }
    
    # Join pre-computed correlation flags
    DF_TYPE <- merge(DF_TYPE, feature_correlation_flags, by = "Feature", all.x = TRUE)
    
    # Write progress to log file
    tryCatch({
      n_done <- length(readLines(progress_file)) - 1L
      cat(paste0("[", Sys.time(), "] Done: ", response,
                 " | Progress: ", n_done, "/", n_genes,
                 " (", round(n_done / n_genes * 100), "%)\n"),
          file = progress_file, append = TRUE)
    }, error = function(e) NULL)
    
    DF_TYPE
  }
  
  message(paste0("[", Sys.time(), "] Finished cancer type: ", ct,
                 " -- all ", n_genes, " genes processed."))
  
  results_per_gene$cancer_type <- ct
  
  # Merge duplicated_names per cancer
  if (!is.null(duplicated_names)) {
    results_per_gene <- dplyr::left_join(
      results_per_gene,
      duplicated_names,
      by = "Feature"
    )
  }
  
  # Save results for this cancer type immediately after completion
  out_file <- paste0(
    "/nfs/casm/team215mg/users/Carmen/projects/REGRESSIONS/output/",
    "biomarker_results_", ct, ".csv"
  )
  write.table(
    results_per_gene,
    file      = out_file,
    sep       = ",",
    col.names = TRUE,
    row.names = FALSE,
    quote     = FALSE
  )
  message(paste0("[", Sys.time(), "] Saved results for ", ct, " to: ", out_file))
}

# Stop parallel cluster cleanly
stopCluster(cl)
message(paste0("[", Sys.time(), "] All cancer types complete."))


################################################################################
# 9. MERGE ALL CANCER TYPE RESULTS
################################################################################

message(paste0("[", Sys.time(), "] Merging all cancer type results..."))

output_dir <- "/output_folder/"

linear_model_results_RESPONSES <- dplyr::bind_rows(
  as.data.frame(fread(paste0(output_dir, "biomarker_results_Colorectal.csv"))),
  as.data.frame(fread(paste0(output_dir, "biomarker_results_Oesophageal.csv"))),
  as.data.frame(fread(paste0(output_dir, "biomarker_results_Gastrointestinal.csv")))
)

message(paste0("[", Sys.time(), "] Merged ", nrow(linear_model_results_RESPONSES), " rows across all cancer types."))


################################################################################
# 10. BIOMARKER CLASSIFICATION
################################################################################

message(paste0("[", Sys.time(), "] Classifying biomarkers..."))

linear_model_results_RESPONSES$biomarker_class <- "No"

linear_model_results_RESPONSES[
  linear_model_results_RESPONSES$P_val_LRT_BH < 0.05 &
    linear_model_results_RESPONSES$MEAN_group2 < -0.5 &
    linear_model_results_RESPONSES$MEAN_group1 > linear_model_results_RESPONSES$MEAN_group2,
]$biomarker_class <- "Yes"

linear_model_results_RESPONSES[
  linear_model_results_RESPONSES$P_val_LRT_BH < 0.05 &
    linear_model_results_RESPONSES$MEAN_group1 < -0.5 &
    linear_model_results_RESPONSES$MEAN_group2 > linear_model_results_RESPONSES$MEAN_group1,
]$biomarker_class <- "Yes"

linear_model_results_RESPONSES$association <- paste(
  linear_model_results_RESPONSES$NAME,
  linear_model_results_RESPONSES$Feature,
  sep = "-"
)

write.table(
  linear_model_results_RESPONSES,
  file      = paste0(output_dir, "all_biomarker-analysis_results.csv"),
  sep       = ",",
  col.names = TRUE,
  row.names = FALSE,
  quote     = FALSE
)


################################################################################
# 11. REFINE BIOMARKER CLASSIFICATION BY EFFECT SIZE
################################################################################

message(paste0("[", Sys.time(), "] Refining classification by effect size..."))

significant <- linear_model_results_RESPONSES[
  linear_model_results_RESPONSES$biomarker_class == "Yes",
]

significant <- significant %>%
  mutate(cohens_f2_quantile_90 = abs(cohens_f2) > quantile(abs(cohens_f2), probs = 0.90, na.rm = TRUE))

significant <- significant %>%
  mutate(cohens_f2_quantile_95 = abs(cohens_f2) > quantile(abs(cohens_f2), probs = 0.95, na.rm = TRUE))

significant$ASSOCIATION_EFFECT <- "Decreased Dep."
significant[significant$b_coef < 0, ]$ASSOCIATION_EFFECT <- "Increased Dep."

significant <- significant[significant$cohens_f2_quantile_90 == TRUE, ]


################################################################################
# 12. CLASSIFY BIOMARKERS INTO PRIORITY TIERS
################################################################################

message(paste0("[", Sys.time(), "] Assigning priority tiers (A/B/C)..."))

significant[abs(significant$MEAN_difference) > 0.5, ]$biomarker_class <- "B"

significant[
  significant$cohens_f2_quantile_95 == TRUE &
    abs(significant$MEAN_difference) > 0.75,
]$biomarker_class <- "A"

significant[significant$biomarker_class == "Yes", ]$biomarker_class <- "C"

write.table(
  significant,
  file      = paste0(output_dir, "supplementary_table_8.1.csv"),
  sep       = ",",
  col.names = TRUE,
  row.names = FALSE,
  quote     = FALSE
)


################################################################################
# 13. PRIORITY SCORING SYSTEM
################################################################################

message(paste0("[", Sys.time(), "] Running priority scoring..."))

cancer_driver_genes <- as.data.frame(
  fread('/path_to_folder/cancer_driver_genes.tsv')
)

binary_features_MSI <- MERGE_binary[
  , colnames(MERGE_binary) %in%
    unique(linear_model_results_RESPONSES[
      linear_model_results_RESPONSES$cancer_type == "Gastrointestinal",
    ]$Feature) | colnames(MERGE_binary) == "sample_ID"
]

binary_features_MSI <- binary_features_MSI[
  binary_features_MSI$sample_ID %in% LFC_filtered$sample_ID,
]

binary_features_MSI <- merge(MSI, binary_features_MSI, by = "sample_ID", all.y = TRUE)

MSI_features <- NULL

for (f in colnames(binary_features_MSI)[3:ncol(binary_features_MSI)]) {
  df_temp     <- binary_features_MSI[, c("sample_ID", "msStatus", f)]
  FISHER.TEST <- fisher.test(df_temp[[2]], df_temp[[3]])
  MSI_features <- rbind(MSI_features, data.frame(f, FISHER.TEST$p.value))
}

MSI_features$p.adj <- p.adjust(MSI_features$FISHER.TEST.p.value, method = "BH")


################################################################################
# HELPER: Priority scoring function
################################################################################

run_priority_scoring <- function(significant_ct) {
  
  result <- NULL
  
  for (x in unique(significant_ct$NAME)) {
    
    df_temp <- significant_ct[significant_ct$NAME == x, ]
    df_temp <- df_temp[order(df_temp$biomarker_class), ]
    
    MAX_biomarker <- min(df_temp$biomarker_class)
    df_temp2      <- df_temp[df_temp$biomarker_class == MAX_biomarker, ]
    
    TYPES_count <- length(unique(df_temp2$Feature_type))
    
    if (TYPES_count > 1) {
      TEMP2    <- df_temp2[df_temp2$Feature_type %in% c("Mut", "CN", "Clin", "Var"), ]
      df_temp3 <- if (nrow(TEMP2) > 0) TEMP2 else df_temp2
    } else {
      df_temp3 <- df_temp2[df_temp2$P_val_LRT_BH == min(df_temp2$P_val_LRT_BH), ][1, ]
    }
    
    N_MAX            <- nrow(df_temp3)
    MAX_pvalue       <- min(df_temp3$P_val_LRT_BH)
    MAX_feature      <- paste(df_temp3[df_temp3$P_val_LRT_BH == MAX_pvalue, ]$Feature,      collapse = "-")
    MAX_feature_type <- paste(df_temp3[df_temp3$P_val_LRT_BH == MAX_pvalue, ]$Feature_type, collapse = "-")
    
    result <- rbind(result, data.frame(df_temp3, MAX_feature, MAX_feature_type, N_MAX))
  }
  
  result
}

format_priority_table <- function(priority_df, cols_to_keep) {
  
  priority_df2 <- priority_df[, cols_to_keep]
  
  priority_df3 <- priority_df2 %>%
    group_by(NAME) %>%
    summarise(across(everything(), ~ paste(na.omit(.), collapse = "/")), .groups = "drop")
  
  dedup_cols <- c("biomarker_class", "delta_MEAN", "MEAN_LFC_depleted_group",
                  "MEAN_LFC_not_depleted_group", "N_depleted", "N_not_depleted",
                  "MAX_feature", "MAX_feature_type")
  
  for (col in intersect(dedup_cols, colnames(priority_df3))) {
    priority_df3[[col]] <- sapply(
      priority_df3[[col]],
      function(x) paste(unique(unlist(strsplit(x, "/"))), collapse = "/")
    )
  }
  
  num_cols <- c("delta_MEAN", "MEAN_LFC_depleted_group", "MEAN_LFC_not_depleted_group",
                "N_depleted", "N_not_depleted")
  for (col in intersect(num_cols, colnames(priority_df3))) {
    priority_df3[[col]] <- suppressWarnings(as.numeric(priority_df3[[col]]))
  }
  
  priority_df3
}

shared_cols <- c(
  "NAME", "MAX_feature", "MAX_feature_type", "Feature", "Feature_type",
  "delta_MEAN", "MEAN_LFC_depleted_group", "MEAN_LFC_not_depleted_group",
  "N_depleted", "N_not_depleted", "biomarker_class", "association"
)


################################################################################
# 14. PRIORITY SCORING - GASTROINTESTINAL CANCERS
################################################################################

message(paste0("[", Sys.time(), "] Priority scoring: Gastrointestinal..."))

significant_GAST  <- significant[significant$cancer_type == "Gastrointestinal", ]
significant_GAST2 <- run_priority_scoring(significant_GAST)

significant_GAST2$MSI_feature <- "No"
significant_GAST2[
  significant_GAST2$Feature %in% MSI_features[MSI_features$p.adj < 0.25, ]$f,
]$MSI_feature <- "Yes"

colnames(differential_dep_gastrointestinal)[1] <- "NAME"

priority_score_GAST  <- merge(differential_dep_gastrointestinal, significant_GAST2, by = "NAME", all.x = TRUE)
priority_score_GAST3 <- format_priority_table(priority_score_GAST, c(shared_cols, "MSI_feature"))

priority_score_GAST3$cancer_driver <- ifelse(priority_score_GAST3$NAME %in% cancer_driver_genes$gene, "Yes", "No")
priority_score_GAST3[priority_score_GAST3$biomarker_class == "", ]$biomarker_class <- "No"

write.table(priority_score_GAST3,
            file = paste0(output_dir, "supplementary_table_8.2.csv"),
            sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)


################################################################################
# 15. PRIORITY SCORING - COLORECTAL CANCERS
################################################################################

message(paste0("[", Sys.time(), "] Priority scoring: Colorectal..."))

significant_COLO  <- significant[significant$cancer_type == "Colorectal", ]
significant_COLO2 <- run_priority_scoring(significant_COLO)

significant_COLO2$MSI_feature <- "No"
significant_COLO2[
  significant_COLO2$Feature %in% MSI_features[MSI_features$p.adj < 0.25, ]$f,
]$MSI_feature <- "Yes"

colnames(differential_dep_COADREAD)[1] <- "NAME"

priority_score_COLO  <- merge(differential_dep_COADREAD, significant_COLO2, by = "NAME", all.x = TRUE)
priority_score_COLO3 <- format_priority_table(priority_score_COLO, c(shared_cols, "MSI_feature"))

priority_score_COLO3$cancer_driver <- ifelse(priority_score_COLO3$NAME %in% cancer_driver_genes$gene, "Yes", "No")
priority_score_COLO3[priority_score_COLO3$biomarker_class == "", ]$biomarker_class <- "No"

write.table(priority_score_COLO3,
            file = paste0(output_dir, "supplementary_table_8.3.csv"),
            sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)


################################################################################
# 16. PRIORITY SCORING - OESOPHAGEAL CANCERS
################################################################################

message(paste0("[", Sys.time(), "] Priority scoring: Oesophageal..."))

significant_OESO  <- significant[significant$cancer_type == "Oesophageal", ]
significant_OESO2 <- run_priority_scoring(significant_OESO)

colnames(differential_dep_ESCA)[1] <- "NAME"

priority_score_OESO  <- merge(differential_dep_ESCA, significant_OESO2, by = "NAME", all.x = TRUE)
priority_score_OESO3 <- format_priority_table(priority_score_OESO, shared_cols)

priority_score_OESO3$cancer_driver <- ifelse(priority_score_OESO3$NAME %in% cancer_driver_genes$gene, "Yes", "No")
priority_score_OESO3[priority_score_OESO3$biomarker_class == "", ]$biomarker_class <- "No"

write.table(priority_score_OESO3,
            file = paste0(output_dir, "supplementary_table_8.4.csv"),
            sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

message(paste0("[", Sys.time(), "] All done."))


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
#    - Includes MSI association flags

# 4. priority_score_COLO3 (supplementary_table_8.3)
#    - Priority-scored biomarkers for colorectal cancers
#    - Includes MSI association flags

# 5. priority_score_OESO3 (supplementary_table_8.4)
#    - Priority-scored biomarkers for esophageal cancers