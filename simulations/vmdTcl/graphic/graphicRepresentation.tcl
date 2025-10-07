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