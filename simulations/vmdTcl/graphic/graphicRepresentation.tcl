# Hide all representations of the given molecule(s)
# Usage: hideAllRep [molid|"all"|"top"] - defaults to "all"
proc hideAllRep {{mol all}} {
    if {$mol eq "all"} {
        # Get all molecule IDs
        set mollist [molinfo list]
        if {[llength $mollist] == 0} {
            puts "No molecules loaded"
            return
        }
        foreach molid $mollist {
            set numrep [molinfo $molid get numreps]
            for {set i 0} {$i < $numrep} {incr i} {
                mol showrep $molid $i off
            }
            puts "All reps hidden for molecule $molid"
        }
        puts "Completed hiding all representations for all [llength $mollist] molecules"
    } elseif {$mol eq "top"} {
        set molid [molinfo top]
        if {$molid == -1} {
            puts "No top molecule available"
            return
        }
        set numrep [molinfo $molid get numreps]
        for {set i 0} {$i < $numrep} {incr i} {
            mol showrep $molid $i off
        }
        puts "All reps hidden for top molecule $molid"
    } else {
        # Specific molecule ID provided
        if {[catch {molinfo $mol get numreps} numrep]} {
            puts "Error: Invalid molecule ID $mol"
            return
        }
        for {set i 0} {$i < $numrep} {incr i} {
            mol showrep $mol $i off
        }
        puts "All reps hidden for molecule $mol"
    }
}

# Show all representations of the given molecule(s)
# Usage: showAllRep [molid|"all"|"top"] - defaults to "all"
proc showAllRep {{mol all}} {
    if {$mol eq "all"} {
        # Get all molecule IDs
        set mollist [molinfo list]
        if {[llength $mollist] == 0} {
            puts "No molecules loaded"
            return
        }
        foreach molid $mollist {
            set numrep [molinfo $molid get numreps]
            for {set i 0} {$i < $numrep} {incr i} {
                mol showrep $molid $i on
            }
            puts "All reps shown for molecule $molid"
        }
        puts "Completed showing all representations for all [llength $mollist] molecules"
    } elseif {$mol eq "top"} {
        set molid [molinfo top]
        if {$molid == -1} {
            puts "No top molecule available"
            return
        }
        set numrep [molinfo $molid get numreps]
        for {set i 0} {$i < $numrep} {incr i} {
            mol showrep $molid $i on
        }
        puts "All reps shown for top molecule $molid"
    } else {
        # Specific molecule ID provided
        if {[catch {molinfo $mol get numreps} numrep]} {
            puts "Error: Invalid molecule ID $mol"
            return
        }
        for {set i 0} {$i < $numrep} {incr i} {
            mol showrep $mol $i on
        }
        puts "All reps shown for molecule $mol"
    }
}


proc hideAllMol {} {
    set mollist [molinfo list]
    if {[llength $mollist] == 0} {
        puts "No molecules loaded"
        return
    }
    foreach molid $mollist {
        mol off $molid
    }
}

proc showAllMol {} {
    set mollist [molinfo list]
    if {[llength $mollist] == 0} {
        puts "No molecules loaded"
        return
    }
    foreach molid $mollist {
        mol on $molid
    }
}


proc showAllComponents {{molid top}} {
    # set color_index 1
    # # Protein representation
    # mol modselect 0 $molid "protein"
    # mol modcolor 0 $molid ColorID $color_index        
    # Ligand representation
    set sel [atomselect $molid "resname AMP or resname ADP or resname ATP or resname GTP"]
    if {[$sel num] > 0} {
        mol selection "resname AMP or resname ADP or resname ATP or resname GTP"
        mol representation Licorice
        mol color ColorID 4
        mol addrep $molid
    }
    
    # MG ions
    set sel [atomselect $molid "resname MG"]
    if {[$sel num] > 0} {
        mol selection "resname MG"
        mol representation VDW
        mol color ColorID 7
        mol addrep $molid
    }

    # ions representation
    set sel [atomselect $molid "ions"]
    if {[$sel num] > 0} {
        mol selection "ions"
        mol representation VDW
        mol color name
        mol addrep $molid
    }

    # Water representation
    set sel [atomselect $molid "water"]
    if {[$sel num] > 0} {
        mol selection "water"
        mol representation lines
        mol color name
        mol addrep $molid
    }

}


proc checkSimulationSetup { {molid top} } {
    # show all components
    showAllComponents $molid
    pbc box
    getSysInfo $molid
}


proc drawCOMofSelection { {molid top} selection } {
    set sel [atomselect $molid $selection]
    set com [measure center $sel weight mass]
    graphics top color red
    graphics top sphere $com radius 1.0
    puts "Center of mass for selection \"$selection\" in molecule $molid: $com"
}