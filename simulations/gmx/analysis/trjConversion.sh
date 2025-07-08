#!/bin/bash
module load GROMACS/2023.3-foss-2023a-CUDA-12.1.1-PLUMED-2.9.0  # Load software module
echo "0 0" | gmx trjconv -s md.tpr -f md.xtc -o mdWrap1ns.xtc -dt 1000 -pbc mol -center
### here input 0 (first) selects protein for centering, 0 (second)  select systems for output
echo "#####################done########################"
