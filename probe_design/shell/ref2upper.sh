#!/bin/bash

# Convert contents of all files in the data/ref/ directory to uppercase
# CAUTION: fasta files in the data/ref/ must be uppercase 
# otherwise the pipeline will not work properly

# Directory path    
dir="data/ref"

# Loop over all files in the directory
for file in "$dir"/*; do
    # Check if it is a file (not a directory)
    if [ -f "$file" ]; then
        # Convert the file content to uppercase and overwrite the file
        tr '[:lower:]' '[:upper:]' < "$file" > temp && mv temp "$file"
    fi
done