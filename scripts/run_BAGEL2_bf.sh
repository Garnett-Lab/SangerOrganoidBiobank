#!/bin/bash

################################################################################
# Script: run_BAGEL2_bf.sh
# Purpose: Calculate Bayes Factors for CRISPR screen data using BAGEL2
# Description: This script runs BAGEL2 (Bayesian Analysis of Gene EssentiaLity)
#              to calculate Bayes Factors for each gene, which quantify the
#              likelihood that a gene is essential vs non-essential
#
# BAGEL2 Reference:
#   Kim E & Hart T (2021). Improved analysis of CRISPR fitness screens and 
#   reduced off-target effects with the BAGEL2 gene essentiality classifier.
#   Genome Medicine 13, 2. https://doi.org/10.1186/s13073-020-00809-3
#
# BAGEL2 Repository:
#   https://github.com/hart-lab/bagel
#
# Installation:
#   git clone https://github.com/hart-lab/bagel.git
#   # Then update MYSCRIPT path below to point to BAGEL.py
#
# Usage:
#   bash run_BAGEL2_bf.sh
#
# Requirements:
#   - Python 3.x
#   - BAGEL2 (BAGEL.py from GitHub repository)
#   - LSF cluster (bsub commands)
#   - Input files: Gene-level fold changes in TSV format
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

## Loop through all TSV files in input folder
# Each file contains gene-level fold changes for one sample
for PATIENT in $(ls -d $INPUT_FOLDER_BF/*tsv); do
  
  ## Extract sample identifier from filename
  # Example: SAMPLE_001.v1_1.gene_LFC.tsv → SAMPLE_001.v1_1
  ID_LIB=$(basename $PATIENT | cut -f1,2 -d'.')
  
  echo "Processing: $ID_LIB"
  
  ## Construct BAGEL2 command
  #
  # Parameters:
  #   bf              - Run Bayes Factor calculation mode
  #   -i $PATIENT     - Input file (gene-level fold changes)
  #   -o ...          - Output file (Bayes Factors for each gene)
  #   -e $ESSENTIAL   - Essential genes (positive training set)
  #   -n $NON_ESSENTIAL - Non-essential genes (negative training set)
  #   -c 1            - Column containing fold changes (0-indexed, so 1 = 2nd column)
  run_BF="python3 ${MYSCRIPT} bf -i $PATIENT -o $OUT_FOLDER_BF/${ID_LIB}.BAGEL2.bf -e $ESSENTIAL -n $NON_ESSENTIAL -c 1"

  echo "  Command: $run_BF"
  
  ## Submit job to LSF cluster
  # Clean up any previous log files for this sample
  rm -f $logs/${ID_LIB}.BF.o.log $logs/${ID_LIB}.BF.e.log
  
  # Submit to cluster with resource requirements
  bsub -q normal \
       -R "select[mem>20000] rusage[mem=20000]" \
       -M 20000 \
       -J ${ID_LIB}.BF \
       -o $logs/${ID_LIB}.BF.o.log \
       -e $logs/${ID_LIB}.BF.e.log \
       "$run_BF"
  
  echo ""
done


################################################################################
# NEXT STEPS
################################################################################

# After all jobs complete, run BAGEL2 'pr' mode to calculate precision-recall.

