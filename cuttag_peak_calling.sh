#!/usr/bin/env bash

set -e
set -o pipefail

PROJ_PATH="$(cd "$(dirname "${BASH_SOURCE}")" && pwd)"
CORES=16

CONDA_PATH=$(conda info --base)
source "$CONDA_PATH/etc/profile.d/conda.sh"

conda activate seacr_env

mkdir -p "${PROJ_PATH}/seacr/bedgraph" "${PROJ_PATH}/seacr/peaks"

# ==========================================
# 1. Convert BAM to BedGraph (bedtools)
# ==========================================
cd "${PROJ_PATH}/alignment/bam"

for bam in *_sorted.bam; do
    [ -e "$bam" ] || continue
    sample=$(basename "$bam" _sorted.bam)
    output_bg="${PROJ_PATH}/seacr/bedgraph/${sample}.bedgraph"
    
    if [ ! -f "$output_bg" ]; then
        bedtools genomecov -ibam "$bam" -bg -pc > "$output_bg"
    fi
done

# ==========================================
# 2. Peak Calling (SEACR) - H2A.Z & IgG Only
# ==========================================
declare -A pairs
pairs=(
    ["NC_H2A.Z_1"]="NC_IgG" 
    ["NC_H2A.Z_2"]="NC_IgG" 
    ["NC_H2A.Z_3"]="NC_IgG" 
    ["NC_H2A.Z_4"]="NC_IgG"
    
    ["si_H2A.Z_1"]="si_IgG" 
    ["si_H2A.Z_2"]="si_IgG" 
    ["si_H2A.Z_3"]="si_IgG" 
    ["si_H2A.Z_4"]="si_IgG"
)

bg_dir="${PROJ_PATH}/seacr/bedgraph"
out_dir="${PROJ_PATH}/seacr/peaks"

if command -v SEACR_1.3.sh &> /dev/null; then
    SEACR_CMD="SEACR_1.3.sh"
elif command -v seacr.sh &> /dev/null; then
    SEACR_CMD="seacr.sh"
else
    echo "Error: SEACR command not found in the current environment." >&2
    exit 1
fi

for treat in "${!pairs[@]}"; do
    control="${pairs[$treat]}"
    
    if [ ! -f "${bg_dir}/${treat}.bedgraph" ] || [ ! -f "${bg_dir}/${control}.bedgraph" ]; then
        echo "Warning: Missing bedgraph for ${treat} or ${control}, skipping..." >&2
        continue
    fi
    
    $SEACR_CMD "${bg_dir}/${treat}.bedgraph" "${bg_dir}/${control}.bedgraph" non relaxed "${out_dir}/${treat}_vs_${control}_seacr"
done
