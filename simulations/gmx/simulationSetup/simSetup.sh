#!/bin/bash

# Usage check
if [ -z "$1" ]; then
    echo "Usage: $0 input.pdb"
    exit 1
fi

# Variables
pdb_name="$1"
pad=1.5
ion_conc=0.15
pos_ion="NA"
neg_ion="CL"

echo "=============================================="
echo "System setup with PDB: $pdb_name"
echo "Padding: $pad nm"
echo "Ion concentration: $ion_conc M with ions: $pos_ion/$neg_ion"
echo "=============================================="

# Prepare working directory
mkdir -p equilibration
cp "$pdb_name" ./equilibration/
cd ./equilibration || exit 1

# Copy force field and MDP files
cp -r ~/myScripts/simulations/gmx/simulationSetup/charmm36-jul2022.ff ./
cp -r ~/myScripts/simulations/gmx/simulationSetup/mdp ./

# Stepwise GROMACS setup
gmx pdb2gmx -f "$pdb_name" -o processed.gro -ff charmm36-jul2022 -water tip3p
gmx editconf -f processed.gro -o boxed.gro -c -d "$pad" -bt cubic
gmx solvate -cp boxed.gro -p topol.top -o pdb_solv.gro
gmx grompp -f ./mdp/min.mdp -c pdb_solv.gro -p topol.top -o ion.tpr -maxwarn 1

# Ion addition with custom ion names
echo SOL | gmx genion -s ion.tpr -o pdb_ion.gro -p topol.top -conc "$ion_conc" -neutral -pname "$pos_ion" -nname "$neg_ion"

# Preprocessing for energy minimization
gmx grompp -f ./mdp/min.mdp -c pdb_ion.gro -p topol.top -o min.tpr -maxwarn 1


echo "=============================================="
echo "System setup with PDB: $pdb_name"
echo "Padding: $pad nm"
echo "Ion concentration: $ion_conc M with ions: $pos_ion/$neg_ion"
echo "Done: Check the file pdb_ion.gro in VMD or GROMACS"
echo "=============================================="

