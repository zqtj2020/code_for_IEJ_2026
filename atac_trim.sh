#!/usr/bin/env bash

ADAPTERS=$(find $CONDA_PREFIX -name "NexteraPE-PE.fa" | head -n 1)

for r1 in *_1.fq.gz; do
    r2="${r1%_1.fq.gz}_2.fq.gz"
    sample="${r1%_1.fq.gz}"
    
    trimmomatic PE -Xmx16G -threads 12 -phred33 \
        "$r1" "$r2" \
        "${sample}_1_trimmed_paired.fq.gz" "${sample}_1_trimmed_unpaired.fq.gz" \
        "${sample}_2_trimmed_paired.fq.gz" "${sample}_2_trimmed_unpaired.fq.gz" \
        ILLUMINACLIP:"${ADAPTERS}":2:30:10:8:TRUE \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:7
done
