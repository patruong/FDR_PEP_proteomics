#!/bin/bash

nextflow run nf-core/quantms --input sdrf_cleaned_local.csv --database 2023-11-01-decoys-reviewed-contam-UP000005640-UP000464024.fas --FDR_level psm-level-fdrs --outdir './results' -work-dir /hd2/tmp -profile docker --max_memory 30GB --max_cpus 8 -resume 
