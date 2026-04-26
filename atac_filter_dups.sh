#!/usr/bin/env bash

PICARD_JAR=$(find $CONDA_PREFIX -name "picard.jar" | head -n 1)
MEMORY="-Xmx20g"

mkdir -p ./mapped

for file in *_sorted.bam; do
    i=${file%_sorted.bam}
    
    samtools view -h "${file}" | grep -E '^@|chr[0-9XY]+' | grep -v -E 'chrM|random|chrUn' | samtools view -b - > "${i}_m.bam"

    java $MEMORY -jar $PICARD_JAR MarkDuplicates \
        I="${i}_m.bam" \
        O="${i}_m_d.bam" \
        M="${i}_dups_metrics.txt" \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=SILENT

    samtools index "${i}_m_d.bam"
    
    mv "${i}_m_d.bam" ./mapped/
    mv "${i}_m_d.bam.bai" ./mapped/
    rm -f "${i}_m.bam"
done