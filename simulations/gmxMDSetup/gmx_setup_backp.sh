#!/bin/bash
echo "Enter pdb file name"
read pdb_name
gmx pdb2gmx -f $pdb_name -o processed.gro -ff charmm36-jul2022 -water tip3p
gmx editconf -f processed.gro -o boxed.gro -c -d 1.5 -bt cubic
gmx solvate -p topol.top -cp boxed.gro -o pdb_solv.gro
gmx grompp -p topol.top -c pdb_solv.gro -o ion.tpr -f ~/scripts/md_setup/mdp/min.mdp -maxwarn 1
echo SOL | gmx genion -s ion.tpr -o pdb_ion.gro -p topol.top -conc 0.15 -neutral
gmx grompp -p topol.top -c pdb_ion.gro -o min.tpr -f ~/scripts/md_setup/mdp/min.mdp -maxwarn 1

echo "Done: check the pdb_ion file in VMD"
