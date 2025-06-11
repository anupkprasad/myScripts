#!/usr/bin/tclsh
# @2022, VMD script for psf genration, solvation and ionisation of system from a pdb
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy

package require autopsf
mol new 381.pdb; ## input pdb file
autopsf -mol top
solvate 381_autopsf.psf 381_autopsf.pdb -t 5 -o water
autoionize -psf water.psf -pdb water.pdb -sc 0.15 -cation SOD -anion CLA -o ionized