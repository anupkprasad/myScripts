#!/usr/bin/tclsh
#### Basic functions to apply on trajectory/frames in VMD @2022
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy

#### Load simulation data  #########

proc trajLoadDCD {{datapath} {pattern "None"} {step 10} {wrap "No"} {sscache "NO"}} {
	set tclscriptpath /media/anupkumar/BackupPlus/scripts/tcl/VMD_scripts
	set psf_name [lindex [glob -directory $datapath ionized.psf] {0}]
	puts $psf_name
	mol new $psf_name
	set prod_list [glob -directory $datapath -type d *{prod,Prod}*]
	if {$pattern != "None"} {
		set prod_list [glob -directory $datapath -type d *{$pattern}*]
	}
	## sorting the prod_list according ASCII
	set prod_list [lsort -dictionary $prod_list]
	foreach prod $prod_list {
		puts $prod
		set filename [lindex [glob -directory $prod/ *.dcd] {0}]
		puts $filename
		mol addfile $filename first 0 last -1 step $step waitfor all
	}
	if {$wrap != "No"} {
		pbc wrap -center com -centersel "protein" -compound fragment -all
	}
	if {$sscache != "No"} {
		source $tclscriptpath/sscache.tcl
	}
	puts "Trajectory is loaded"
}


proc trajLoadXtc {{filename "md_frm1ns.xtc" } {step 1}} {
	mol new ../equilibration/min.gro
	mol addfile ./$filename first 0 last -1 step $step waitfor all
	puts "Trajectory is loaded"
}



###### Load AlphaFold structure ###########

proc loadPdbAll {} {
	set pdb_list [glob *pdb]
	set i 0
	foreach pdb $pdb_list {
		if {$i == 0} {
			mol new $pdb
		} else {
			mol addfile $pdb
		}
		incr i
	}
}

proc loadCifAll {} {
	set pdb_list [glob *cif]
	set i 0
	foreach pdb $pdb_list {
		if {$i == 0} {
			mol new $pdb
		} else {
			mol addfile $pdb
		}
		incr i
	}
}

proc loadCifAF3 {} {
	set dir [lsort -increasing [glob seed*]]
	set i 1
	foreach d $dir {
		set filecif [glob ./$d/*.cif]
		puts $d
		if {$i ==1} {
			mol new $filecif
		} else {
				mol addfile $filecif
		}
		incr i
		}
}


proc loadCifAF3Multijobs {{prefix ""} {refStr "${refStructures}/DCLK1_KD_tail.pdb"} {KDResids "0 to 258"}} {
	### load and represent refernce structure 
	mol new $refStr
	mol modselect 0 top "protein and resid ${KDResids}"
	mol modcolor 0 top ColorID 8
    # Create a new representation for non-KD
	mol selection "protein and not resid ${KDResids}"
	mol representation NewCartoon
	mol color ColorID 2
	mol addrep top

    # Set the base path
    set base_path [pwd]
    set base_path [file normalize "."]
    puts "Base path is: $base_path"

    # Get the list of all entries in the directory
    set entries [glob -nocomplain -directory $base_path *]
    # Filter only directories with the given prefix
    set directories {}
    foreach entry $entries {
        if {[file isdirectory $entry]} {
            set dir_name [file tail $entry]
            if {[string match ${prefix}* $dir_name]} {
                lappend directories $entry
            }
        }
    }

    # Sort for deterministic order (optional)
    set directories [lsort -dictionary $directories]

    # Loop over all matched directories
    set i 1
    foreach dir $directories {
        puts "Loading directory: $dir"
        cd $dir
        if {$i == 1} {
            loadCifAF3
        } else {
            loadCifAF3
        }
        trajAlignTwoMol 0 top "backbone and resid ${KDResids}"
        #mol showrep top 0 off
        ### change default represention to KD selection only
        mol modselect 0 top "protein and resid ${KDResids}"
        # Create a new representation for non-KD
		mol selection "protein and not resid ${KDResids}"
		mol representation NewCartoon
		mol color ColorID 4
		mol addrep top
		### diplay of
        mol off top
        cd $base_path
        incr i
    }
}

proc addCifAF3 {{rename ""}} {
    set dir [lsort -increasing [glob seed*]]
    foreach d $dir {
        set filecif [glob ./$d/*.cif]
        puts $d
        mol addfile $filecif
    }

    set molid [molinfo top]

    ;# Only rename if a non-empty string was passed
    if {$rename ne ""} {
        mol rename $molid $rename
    }
}
