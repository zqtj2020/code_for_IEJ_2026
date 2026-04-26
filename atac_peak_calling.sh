#!/usr/bin/env bash

cd ../bed
mkdir -p ../shift
mkdir -p ../peak

for sample in *.bed; do
    describer=${sample%.bed}
    
    awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 4, $3 + 4, $4, $5, $6; else print $1, $2 - 5, $3 - 5, $4, $5, $6}' "$sample" > "../shift/${describer}_shifted.bed"
done
mv *shifted.bed ./shift

cd ../shift

for foo in *_shifted.bed; do
    bar=${foo%_shifted.bed}
    
    macs2 callpeak --nomodel \
        -t "$foo" \
        -n "$bar" \
        --nolambda \
        --gsize 1.87e9 \
        --keep-dup all \
        --slocal 10000
done

mv *_summits.bed ./peak
mv *.narrowPeak ./peak
mv *.xls ./peak