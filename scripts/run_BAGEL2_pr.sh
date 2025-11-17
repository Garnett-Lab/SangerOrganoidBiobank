#!/bin/bash

################################################################################
# Script: run_BAGEL2_pr.sh
# Purpose: Calculate Precision-Recall and FDR values using BAGEL2
# Description: This script runs BAGEL2 'pr' mode (Precision-Recall) to convert
#              Bayes Factors into False Discovery Rates (FDR) and generate
#              precision-recall curves for quality assessment
#
# Prerequisites:
#   - BAGEL2 Bayes Factors must be calculated first (run_BAGEL2_bf.sh)
#   - BAGEL2 installed from: https://github.com/hart-lab/bagel
#
# BAGEL2 Reference:
#   Kim E & Hart T (2021). Improved analysis of CRISPR fitness screens and 
#   reduced off-target effects with the BAGEL2 gene essentiality classifier.
#   Genome Medicine 13, 2. https://doi.org/10.1186/s13073-020-00809-3
#
# Usage:
#   bash run_BAGEL2_pr.sh
#
# Requirements:
#   - Python 3.x
#   - BAGEL2 (BAGEL.py from GitHub)
#   - LSF cluster (bsub)
#   - Input: .bf files from BAGEL2 bf mode
################################################################################

################################################################################
# CONFIGURATION
################################################################################

## Path to BAGEL.py script
# IMPORTANT: Download BAGEL2 from https://github.com/hart-lab/bagel
# Then update this path to point to BAGEL.py in your local installation
MYSCRIPT=/path_to_folder/BAGEL.py

## Input/Output directories
INPUT_FOLDER_BF=/path_to_folder/input_bf
OUT_FOLDER_BF=/path_to_folder/output_bf
logs=/path_to_folder/logs

## Control gene lists
# Essential genes: Should show strong depletion (negative fold changes)
# Used as positive training set for Bayesian classifier
ESSENTIAL=/path_to_folder/curated_BAGEL_essential_genes_organoids_downsampling.txt

# Non-essential genes: Should show no effect (fold changes near zero)
# Used as negative training set for Bayesian classifier
NON_ESSENTIAL=/path_to_folder/curated_BAGEL_non-essential_genes_organoids_downsampling.txt


################################################################################
# PROCESS EACH SAMPLE
################################################################################

## Loop through all .bf files (Bayes Factor files from previous step)
# Each .bf file contains Bayes Factors for all genes in one sample
for PATIENT in $(ls -d $INPUT_FOLDER_PR/*.bf); do
  
  ## Extract sample identifier from filename
  # Example: SAMPLE_001.v1_1.BAGEL2.bf → SAMPLE_001.v1_1
  ID_LIB=$(basename $PATIENT | cut -f1,2 -d'.')
  
  echo "Processing: $ID_LIB"
  
  ## Construct BAGEL2 precision-recall command
  #
  # Parameters:
  #   pr              - Run Precision-Recall mode
  #   -i $PATIENT     - Input file (Bayes Factors from 'bf' mode)
  #   -o ...          - Output file (Precision-Recall with FDR values)
  #   -e $ESSENTIAL   - Essential genes (calculate true positive rate)
  #   -n $NON_ESSENTIAL - Non-essential genes (calculate false positive rate)
  run_PR="python3 ${MYSCRIPT} pr -i $PATIENT -o $OUT_FOLDER_PR/${ID_LIB}.PRECISION-RECALL.PC -e $ESSENTIAL -n $NON_ESSENTIAL"

  echo "  Command: $run_PR"
  
  ## Submit job to LSF cluster
  # Clean up any previous log files for this sample
  rm -f $logs/${ID_LIB}.PR.o.log $logs/${ID_LIB}.PR.e.log
  
  # Submit to cluster with resource requirements
  bsub -q normal \
       -R "select[mem>20000] rusage[mem=20000]" \
       -M 20000 \
       -J ${ID_LIB}.PR \
       -o $logs/${ID_LIB}.PR.o.log \
       -e $logs/${ID_LIB}.PR.e.log \
       "$run_PR"
  
  echo ""
done


################################################################################
# NEXT STEPS
################################################################################

# After BAGEL2 analysis completes, import results into R (CRISPR_pipeline.R).

