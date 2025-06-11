#!/usr/bin/tclsh
#### Functions to check clashes within moluclular atoms  VMD @2022
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy

proc get_clash { {molid top} {cutoff 1.5} {rep "No"} }  {
    set molecule_number [molinfo $molid get id]
    set clash_atoms [measure contacts $cutoff [atomselect $molid noh]]
    set atom_num [llength [lindex $clash_atoms 0]]
    if {$atom_num != 0} {
        puts "Total clash_atoms =  $atom_num"
        puts "following atom indices: $clash_atoms"
    } else {
        puts "Total clash_atoms =  No Clash in Atoms"
    }

    if {($rep == "Yes") && ($atom_num > 0)} {
        puts $atom_num
        mol addrep $molecule_number
        set repid [expr [molinfo $molid get numreps] - 1]
        mol modselect $repid $molecule_number "index [lsort -unique [lindex $clash_atoms 0]]"
        mol modstyle $repid $molecule_number vdw
        mol modcolor $repid $molecule_number Chain

        mol addrep $molecule_number
        set repid [expr [molinfo $molid get numreps] - 1]
        mol modselect $repid $molecule_number "index [lsort -unique [lindex $clash_atoms 1]]"
        mol modstyle $repid $molecule_number vdw
        mol modcolor $repid $molecule_number Name

    }
    return $atom_num
}
