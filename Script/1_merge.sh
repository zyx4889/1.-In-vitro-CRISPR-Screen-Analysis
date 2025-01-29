#!/bin/bash

# Create analysis directory for merged files
mkdir -p ./analysis

# Loop through all samples - note the full experiment design:
# T0_POP_S1: Initial population sample
# T12_DMSO_POP_S2: 12-day DMSO control population
# T12_MYCi_POP_S3: 12-day MYCi975 treated population
# T0_CLONE_S4: Initial clonal sample
# T12_CLONE_S5: 12-day clonal control
# T12_CLONE_S6: 12-day clonal treatment
for sample in T0_POP_S1 T12_DMSO_POP_S2 T12_MYCi_POP_S3 T0_CLONE_S4 T12_CLONE_S5 T12_CLONE_S6; do
    # Merge all lanes (L001-L004) for each sample using zcat for gzipped files
    # The * wildcard captures all L00X lane files for each sample
    # Output is re-compressed with gzip (-c flag writes to stdout)
    zcat ./BaseCalls/${sample}*.fastq.gz | gzip -c > ./analysis/${sample}.fastq.gz && \
    echo "✓ Merged ${sample}"
done