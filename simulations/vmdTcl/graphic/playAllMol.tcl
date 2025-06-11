#!/usr/bin/tclsh
#### Basic functions to apply on trajectory/frames in VMD @2022
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy

proc playFrame {} {
    if {$::frame < $::nframes} {
        animate goto $::frame
        incr ::frame
        after 100 update
        after 500 playFrame
    } else {
        set molid [lindex $::mol_list $::mol_index]
        puts "Finished molecule $molid"
        mol off $molid
        incr ::mol_index
        after 100 playEachMolecule
    }
}

proc playEachMolecule {} {
    if {$::mol_index >= [llength $::mol_list]} {
        puts "Done playing all molecules."
        return
    }

    set molid [lindex $::mol_list $::mol_index]
    mol top $molid
    mol on $molid

    set nreps [molinfo $molid get numreps]
    for {set i 0} {$i < $nreps} {incr i} {
        mol showrep $molid $i on
    }

    set ::frame 0
    set ::nframes [molinfo $molid get numframes]

    if {$::nframes < 2} {
        puts "Skipping molecule $molid (no trajectory)"
        incr ::mol_index
        after 100 playEachMolecule
        return
    }

    puts "Playing trajectory for molecule $molid..."

    playFrame
}


proc playAllMol {} {
	# Global control variables
	set ::mol_index 0
	set ::mol_list [molinfo list]
	# Start the playback
	playEachMolecule
}