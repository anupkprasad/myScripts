#!/bin/bash

# Variables
gro_name="prot_lig_ion.gro"
pad=1.2
ion_conc=0.15
pos_ion="NA"
neg_ion="CL"


## convert pdbtogro and merge gro files all components
#gmx insert-molecules -f protein_processed.gro -ci ligand_h_GMX.gro -nmol 1 -o prot_lig.gro
gmx editconf -f ion.pdb -o ion.gro
#gmx insert-molecules -f prot_lig.gro -ci mn.gro -nmol 1 -o prot_lig_mn.gro

source /home/anup/workspace/project1/env-md/bin/activate
python /home/anup/myScripts/simulations/mergeGroFiles.py



# Prepare working directory
mkdir -p equilibration
cp "$gro_name" ./equilibration/
cp "topol.top" ./equilibration/
cp "ligand_h_GMX.itp" ./equilibration/
cp -r ./amber99sb-ildn.ff ./equilibration/
cp "posre.itp" ./equilibration/
cd ./equilibration || exit 1
cp -r ~/myScripts/simulations/gmx/simulationSetup/mdp ./


### solvation and system setup
gmx editconf -f "$gro_name" -o boxed.gro -c -d "$pad" -bt cubic
gmx solvate -cp boxed.gro -p topol.top -o pdb_solv.gro

gmx grompp -f ./mdp/min.mdp -c pdb_solv.gro -p topol.top -o ion.tpr -maxwarn 1

# Ion addition with custom ion names
echo SOL | gmx genion -s ion.tpr -o pdb_ion.gro -p topol.top -conc "$ion_conc" -neutral -pname "$pos_ion" -nname "$neg_ion"

# Preprocessing for energy minimization
gmx grompp -f ./mdp/min.mdp -c pdb_ion.gro -p topol.top -o min.tpr -maxwarn 1


python /home/anup/myScripts/simulations/createTopology.py topol.top pdb_ion.gro
deactivate