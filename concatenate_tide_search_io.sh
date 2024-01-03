#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_directory output_file"
    exit 1
fi

# Get the input directory and output file from the command line arguments
FILEDIR=$1
OUTFILE=$2

echo "Concatenating tide-search results"

# Get the list of folders in the input directory
FOLDERS=($FILEDIR/*)

# Use the first file to create or overwrite the output file
echo "Concatenate tide-search results: ${FOLDERS[0]} into $OUTFILE"
cat "${FOLDERS[0]}/tide-search.txt" > $OUTFILE

# Loop over the rest of the folders
for folder in "${FOLDERS[@]:1}"; do
    echo "Concatenate tide-search results: ${folder} into $OUTFILE"
    tail -n +2 "$folder/tide-search.txt" >> $OUTFILE
done