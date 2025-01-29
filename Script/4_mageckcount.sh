#!/bin/bash

# First, install GNU parallel if you haven't:
# sudo apt-get install parallel

# Create a file list
echo "T0_rep1.fastq.gz
T0_rep2.fastq.gz
T12_rep1_DMSO.fastq.gz
T12_rep2_DMSO.fastq.gz
T12_rep1_MYCi.fastq.gz
T12_rep2_MYCi.fastq.gz" > fastq_files.txt

# Set number of cores to use (adjust based on your system)
N_CORES=6

# Run MAGeCK count with parallel processing
mageck count \
-l lib.csv \
-n MYCi_screen \
--sample-label "T0_rep1,T0_rep2,T12_DMSO_rep1,T12_DMSO_rep2,T12_MYCi_rep1,T12_MYCi_rep2" \
--fastq $(cat fastq_files.txt | parallel -j $N_CORES echo {}) \
--norm-method median \
--pdf-report \
--count-n \
--unmapped-to-file \

