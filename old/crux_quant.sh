#!/bin/bash

FASTA=2023-11-01-decoys-reviewed-contam-UP000005640-UP000464024.fas
WORKDIR=reanalysis

mkdir -p $WORKDIR/mzML/bullseye

for file in $WORKDIR/mzML/*; do
  echo "$file"
done

# Crux bullseye
crux bullseye $WORKDIR/mzML/NEG1.mzML $WORKDIR/mzML/NEG1.mzML --output-dir $WORKDIR