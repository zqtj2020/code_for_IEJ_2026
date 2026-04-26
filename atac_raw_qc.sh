#!/usr/bin/env bash
mkdir -p ../raw_qc
fastqc *.fq.gz -o ../raw_qc
multiqc ../raw_qc -n raw_multiqc_report -o ../raw_qc