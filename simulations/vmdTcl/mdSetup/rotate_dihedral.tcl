#!/usr/bin/tclsh
# @2024, VMD script to apply constraints
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy

## want to rotate about lat say "a-b" bond, then need to select all atoms of one side let say $moveselection
# to rotate along mentioned bond

set a [atomselect top "index index_of_a"]
set b [atomselect top "index index_of_b"]
$moveselection move [trans bond [lindex [$a1 get {x y z}] 0] [lindex [$a2 get {x y z}] 0] 10 deg]

#### take images of molecule rotating about a bond
set movesel [atomselect top "index 1 2 5 6 7 8 9"]
for {set i 0} {$i <360 } {incr i 10}
{
$movesel move [trans bond [lindex [$a1 get {x y z}] 0] [lindex [$a2 get {x y z}] 0] 10 deg]
render TachyonInternal [format "rotate%03d.tga" $i]
}

