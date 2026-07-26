#!/usr/bin/env bash

set -e
set -o pipefail

PROJ_PATH="$(cd "$(dirname "${BASH_SOURCE}")" && pwd)"
CORES=16

CONDA_PATH=$(conda info --base)
source "$CONDA_PATH/etc/profile.d/conda.sh"

conda activate map_env

mkdir -p "${PROJ_PATH}/alignment/bigwig"

cd "${PROJ_PATH}/alignment/bam"

for bam in *_sorted.bam; do
    [ -e "$bam" ] || continue
    
    sample=$(basename "$bam" _sorted.bam)
    log_file="${PROJ_PATH}/alignment/spikein_log/${sample}_spikein_bowtie2.txt"
    output_bw="${PROJ_PATH}/alignment/bigwig/${sample}_spikein_norm.bw"
    
    if [ ! -f "$log_file" ]; then
        echo "Warning: Spike-in log not found for ${sample}, skipping..." >&2
        continue
    fi
    
    if [ ! -f "$output_bw" ]; then
        align_1=$(grep "aligned concordantly exactly 1 time" "$log_file" | awk '{print $1}') || true
        align_gt1=$(grep "aligned concordantly >1 times" "$log_file" | awk '{print $1}') || true
        
        align_1=${align_1:-0}
        align_gt1=${align_gt1:-0}
        
        spike_reads=$((align_1 + align_gt1))
        
        if [ "$spike_reads" -eq 0 ]; then
            echo "Warning: Spike-in reads for ${sample} is 0, skipping..." >&2
            continue
        fi
        
        scale_factor=$(echo "scale=6; 100000 / $spike_reads" | bc)
        
        bamCoverage -b "$bam" \
                    -o "$output_bw" \
                    --scaleFactor "${scale_factor}" \
                    --binSize 1 \
                    --extendReads \
                    -p "$CORES"
    fi
done
