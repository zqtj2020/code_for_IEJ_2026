#!/usr/bin/env bash

set -e
set -o pipefail

PROJ_PATH="$(cd "$(dirname "${BASH_SOURCE}")" && pwd)"
SPIKE_REF="${PROJ_PATH}/genome/dm6_index/dm6"
CORES=16

CONDA_PATH=$(conda info --base)
source "$CONDA_PATH/etc/profile.d/conda.sh"

mkdir -p "${PROJ_PATH}/alignment/spikein_bam" "${PROJ_PATH}/alignment/spikein_log"

conda activate map_env

cd "${PROJ_PATH}/clean_data"

for r1 in *_clean_R1.fq.gz; do
    [ -e "$r1" ] || continue
    
    r2="${r1/_clean_R1.fq.gz/_clean_R2.fq.gz}"
    sampleName=$(basename "${r1}" _clean_R1.fq.gz)
    
    output_bam="${PROJ_PATH}/alignment/spikein_bam/${sampleName}_spikein_sorted.bam"
    log_file="${PROJ_PATH}/alignment/spikein_log/${sampleName}_spikein_bowtie2.txt"
    
    if [ ! -f "$output_bam" ]; then
        bowtie2 --local --very-sensitive --no-mixed --no-discordant \
            --no-overlap --no-dovetail \
            --phred33 -I 10 -X 700 -p "$CORES" \
            -x "$SPIKE_REF" -1 "$r1" -2 "$r2" 2> "$log_file" | \
            samtools sort -@ 8 -o "$output_bam" -

        samtools index "$output_bam"
    fi
done
