#!/usr/bin/env bash
mkdir -p ../clean_qc
fastqc *trimmed_paired.fq.gz -o ../clean_qc
multiqc ../clean_qc -n clean_multiqc_report -o ../clean_qc
