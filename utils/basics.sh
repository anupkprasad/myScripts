#!/bin/bash

for file in *.pdb; do
    # Extract the filename without extension
    base="${file%.pdb}"

    echo "Processing $file"
    # Run the reduce command
    reduce -NOFLIP "${file}" > "${base}_H.pdb"
done