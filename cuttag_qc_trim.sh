#!/usr/bin/env bash

set -e
set -o pipefail

PROJ_PATH="$(cd "$(dirname "${BASH_SOURCE}")" && pwd)"
DATA_DIR="${PROJ_PATH}/raw_data"
CORES=16

CONDA_PATH=$(conda info --base)
source "$CONDA_PATH/etc/profile.d/conda.sh"

mkdir -p "${PROJ_PATH}/clean_data" "${PROJ_PATH}/reports"

conda activate qc_env

if [ -d "$DATA_DIR" ]; then
    cd "$DATA_DIR"
fi

for r1 in *.R1.fastq.gz; do
    [ -e "$r1" ] || continue
    
    r2="${r1/.R1./.R2.}"
    sample=$(basename "$r1" .R1.fastq.gz)
    
    clean_r1="${PROJ_PATH}/clean_data/${sample}_clean_R1.fq.gz"
    clean_r2="${PROJ_PATH}/clean_data/${sample}_clean_R2.fq.gz"
    
    if [ ! -f "$clean_r1" ]; then
        fastp -i "$r1" -I "$r2" \
              -o "$clean_r1" \
              -O "$clean_r2" \
              -h "${PROJ_PATH}/reports/${sample}.html" \
              -j "${PROJ_PATH}/reports/${sample}.json" \
              --thread "$CORES" --detect_adapter_for_pe --trim_poly_g
    fi
done
