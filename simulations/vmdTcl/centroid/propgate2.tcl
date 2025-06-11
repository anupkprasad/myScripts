#!/usr/bin/tclsh
#### calculates propgating structures from list of pdbs @2023
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy
proc delete_dockedfiles { } {
    if {[catch {glob -directory ./ *docked.pdb} err opts] == 0} {
    	puts "test"
    	set docked_pdbs [glob -directory ./ *docked.pdb]
    	foreach pdb $docked_pdbs {
    	exec rm $pdb	
    	}
    }
    puts "dele done"
}


proc mkdirifnotexist { dirname } {
	if {[file exist $dirname]} {
    # check that it's a directory
    if {! [file isdirectory $dirname]} {
        puts "$dirname exists, but it's a file"
    }
	} else {
    file mkdir $dirname
	}
}


source "$pathVmdScript/centroid/propgate1.tcl"
proc getPropgatingPdbs {} {
	delete_dockedfiles; # delete unwanted pdb files like docked.pdb
	mkdirifnotexist prop_str
	set pdb_list [glob *pdb]
	foreach pdb $pdb_list {
		set pdbprefix [lindex [split $pdb "."] 0]
		puts $pdb
		getPropPdb $pdbprefix 8 0
	}

	cd ./prop_str/
	set pdb_list [glob *pdb]
	foreach pdb $pdb_list {
		mol new $pdb
	}
}


