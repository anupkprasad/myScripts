#!/bin/bash

# Variables
gro_name="ligand_h_GMX.gro"
top_name="ligand_h_GMX.top"
itp_name="ligand_h_GMX.itp"
pad=1.2
ion_conc=0.15
pos_ion="NA"
neg_ion="CL"


# Prepare working directory
mkdir -p equilibration
cp "$gro_name" ./equilibration/
cp "$top_name" ./equilibration/
cp "$itp_name" ./equilibration/
cd ./equilibration || exit 1
cp -r ~/myScripts/simulations/gmx/simulationSetup/mdp_lig/mdp ./
cp ~/myScripts/simulations/gmx/simulationSetup/equil_spl2.sh ./

### rename topology file to topol.top for pipeline compatibility
mv "$top_name" topol.top
top_name="topol.top"


### solvation and system setup
gmx editconf -f "$gro_name" -o boxed.gro -c -d "$pad" -bt cubic
gmx solvate -cp boxed.gro -p "$top_name" -o pdb_solv.gro
gmx grompp -f ./mdp/min.mdp -c pdb_solv.gro -p "$top_name" -o ion.tpr -maxwarn 1

# Ion addition with custom ion names
echo SOL | gmx genion -s ion.tpr -o pdb_ion.gro -p "$top_name" -conc "$ion_conc" -neutral -pname "$pos_ion" -nname "$neg_ion"

# Preprocessing for energy minimization
gmx grompp -f ./mdp/min.mdp -c pdb_ion.gro -p "$top_name" -o min.tpr -maxwarn 1

### generate topology file for AMBER
source /home/anup/workspace/project1/env-md/bin/activate
python /home/anup/myScripts/simulations/createTopology.py "$top_name" "pdb_ion.gro"
deactivate