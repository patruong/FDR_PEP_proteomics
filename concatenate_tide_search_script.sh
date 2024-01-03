#!/bin/bash

echo "Concatenating tide-search results"

FILEDIR=reanalysis/tide-search
FOLDERS=($FILEDIR/*)

echo "Concatenate tide-search results: ${FOLDERS[0]} into combined.tsv"
cat "${FOLDERS[0]}/tide-search.txt" > combined.tsv

for folder in "${FOLDERS[@]:1}"; do
    echo "Concatenate tide-search results: ${folder} into combined.tsv"
    tail -n +2 "$folder/tide-search.txt" >> combined.tsv
done
