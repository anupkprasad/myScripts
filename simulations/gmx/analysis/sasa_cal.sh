#!/bin/bash

# Stop on any error, but don't exit completely
set -u
set -o pipefail

# Input files
file_suffix="r_395_r_479"
GRO_FILE="../system.pdb"  # Adjust if needed
XTC_FILE="./md.xtc"
NDX_FILE="index_$file_suffix.ndx"
OUTPUT_FILE="sasa_$file_suffix.xvg"

# Prompt user for residue selection
echo "Enter the residue selection for SASA (e.g., resid 395 to 479):"
read RESID_SELECTION

# Generate index file for selected residues
echo -e "$RESID_SELECTION\nq" | gmx make_ndx -f "$GRO_FILE" -o "$NDX_FILE"

# Ensure the index file was created
if [ ! -f "$NDX_FILE" ]; then
    echo "Error: Index file $NDX_FILE was not created. Check your residue selection."
    exit 1
fi

# Extract the group number assigned by GROMACS
GROUP_NUMBER=$(grep -B1 "$RESID_SELECTION" "$NDX_FILE" | head -1 | awk '{print $1}')

if [ -z "$GROUP_NUMBER" ]; then
    echo "Error: Could not determine the correct index group for SASA calculation."
    exit 1
fi

# Run SASA calculation with the correct group number
gmx sasa -s "$GRO_FILE" -f "$XTC_FILE" -n "$NDX_FILE" -o "$OUTPUT_FILE" -surface 1 -output "$GROUP_NUMBER"

echo "SASA calculation complete. Output saved in $OUTPUT_FILE."
