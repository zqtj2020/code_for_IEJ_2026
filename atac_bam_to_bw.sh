#!/usr/bin/env bash
mkdir -p ../bw
mkdir -p ../bed

for sample in *_m_d.bam; do
    describer=${sample%_m_d.bam}
    
    bamCoverage --bam "${sample}" \
                -o "${describer}.bw" \
                --binSize 10 \
                --normalizeUsing RPGC \
                --effectiveGenomeSize 2308125349 \
                --ignoreForNormalization chrX chrM \
                --extendReads \
                -p 16

    bedtools bamtobed -i "$sample" > "../bed/${describer}.bed"
done

mv *.bw ../bw