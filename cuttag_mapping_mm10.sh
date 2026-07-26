#!/usr/bin/env bash

set -e
set -o pipefail

PROJ_PATH="$(cd "$(dirname "${BASH_SOURCE}")" && pwd)"
REF_INDEX="${PROJ_PATH}/genome/mm10_index/mm10"
CORES=16

CONDA_PATH=$(conda info --base)
source "$CONDA_PATH/etc/profile.d/conda.sh"

mkdir -p "${PROJ_PATH}/alignment/bam" "${PROJ_PATH}/alignment/log"

conda activate map_env

cd "${PROJ_PATH}/clean_data"

for r1 in *_clean_R1.fq.gz; do
    [ -e "$r1" ] || continue
    
    r2="${r1/_clean_R1.fq.gz/_clean_R2.fq.gz}"
    sampleName=$(basename "${r1}" _clean_R1.fq.gz)
    
    output_bam="${PROJ_PATH}/alignment/bam/${sampleName}_sorted.bam"
    log_file="${PROJ_PATH}/alignment/log/${sampleName}_bowtie2.txt"
    
    if [ ! -f "$output_bam" ]; then
        bowtie2 --local --very-sensitive --no-mixed --no-discordant \
            --phred33 -I 10 -X 700 -p "$CORES" \
            -x "$REF_INDEX" -1 "$r1" -2 "$r2" 2> "$log_file" | \
            samtools sort -@ 8 -o "$output_bam" -

        samtools index "$output_bam"
    fi
done
