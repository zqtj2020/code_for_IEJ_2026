#!/usr/bin/env bash

for r1 in *_1_trimmed_paired.fq.gz; do
    r2="${r1%_1_trimmed_paired.fq.gz}_2_trimmed_paired.fq.gz"
    sample="${r1%_1_trimmed_paired.fq.gz}"

    bowtie2 -x /home/graceatac/genome/mm10_index/mm10 \
        -1 "$r1" \
        -2 "$r2" \
        -p 12 -X 1500 \
        --rg-id "$sample" --rg "SM:$sample" --rg "LB:ATAC-seq" --rg "PL:ILLUMINA" \
        2> "${sample}_bowtie2.log" | \
        samtools sort -@ 6 -m 2G -o "${sample}_sorted.bam" -

    samtools index "${sample}_sorted.bam"
done
