#!/bin/bash
module load GROMACS/2023.3-foss-2023a-CUDA-12.1.1-PLUMED-2.9.0  # Load software module
### xtc conversion: every frame is 1ns, here dt 1000ps == 1ns
### In gmx trjconv,  input 1 select protein for centering and 0  select system to be output



# List of main directories
main_dirs=( "dclk3_separatedLobes")

# Iterate over each main directory
for main_dir in "${main_dirs[@]}"; do
    # Iterate over each subdirectory (rep1, rep2, rep3)
    for rep in rep1 rep2 rep3; do
        # Define the full path
        sub_dir="${main_dir}/${rep}"

        # Check if the directory exists before running commands
        if [ -d "$sub_dir" ]; then
            echo "Entering directory: $sub_dir"
            cd "$sub_dir" || { echo "Failed to enter $sub_dir"; exit 1; }

            # Run gmx trjconv command and wait for it to complete
            echo "Running gmx trjconv in $sub_dir..."
	    source /home/akp66103/gmx_script/conv_sim.sh  <<EOF
1
EOF
            wait  # Ensures gmx trjconv finishes before proceeding

            echo "Completed gmx trjconv in $sub_dir"
            cd - > /dev/null  # Go back to the original directory
        else
            echo "Directory $sub_dir does not exist. Skipping..."
        fi
    done
done

echo "All tasks completed!"

