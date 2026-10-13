#!/usr/bin/tclsh
# @2020, VMD script to apply constraints
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy

set all [atomselect top all] 
set p [atomselect top "protein and backbone"] 
$all set beta 0 
$p set beta 1
$all writepdb fixed.pdb 

set all [atomselect top all] 
set q [atomselect top "protein and name CA"] 
$all set beta 0 
$q set beta 0.5
$all writepdb constrained.pdb
