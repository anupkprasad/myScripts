#!/usr/bin/tclsh
#### Basic functions to apply on trajectory/frames in VMD @2022
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy

proc trajLoad {{datapath} {pattern "None"} {step 10} {wrap "No"} {sscache "NO"}} {
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


proc trajWritepdb {{molid top} {sel "protein"} {frames "None"}} {
	## write pdb for trajectory
	set sel [atomselect $molid $sel]
	set n [molinfo $molid get numframes]
	for { set i 0 } { $i < $n } { incr i } {

		$sel frame $i
		$sel writepdb $i.pdb
	}
}



proc trajAlign {{molid top} {sel "protein"}} {
	## sel: segment (Kinase domain) which are required to align, allsel: transformation matrix to entire molecules (Kinase domaain +tail)
	set sel0 [atomselect $molid $sel]
	set n [molinfo $molid get numframes]
	for { set i 0 } { $i < $n } { incr i } {
		set sel1 [atomselect $molid $sel]
		set allsel [atomselect $molid "protein" frame $i]
		$sel1 frame $i
		set M [measure fit $sel1 $sel0]
		$allsel move $M
	}
}


proc TrajAlignTwoMol {{molidRef top} {molid top} {sel "protein and backbone"}} {
	## sel: segment (Kinase domain) which are required to align, allsel: transformation matrix to entire molecules (Kinase domaain +tail)
	set sel0 [atomselect $molidRef $sel]
	set n [molinfo $molid get numframes]
	for { set i 0 } { $i < $n } { incr i } {
		set sel1 [atomselect $molid $sel]
		set allsel [atomselect $molid all frame $i]
		$sel1 frame $i
		set M [measure fit $sel1 $sel0]
		$allsel move $M
	}
}



proc alignMultimolid {{molid_refid 0} {molid_align top} {sel "protein and backbone"}} {
	### mol_align : the molids which are to be aligned
	set sel0 [atomselect $molid_refid $sel]
	foreach al $molid_align {
		set sel1 [atomselect $al $sel]
		set M [measure fit $sel1 $sel0]
		$sel1 move $M
	}
}


proc getAminoseq {{molid top} {sel "protein"}} {
	set aa_triple { ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR}
	set aa_single { A C D E F G H I K L M N P Q R S T V W Y }
	set sel [atomselect $molid $sel]
	set uresid [lsort -unique [$sel get resid]]
	## sorting the prod_list according ASCII
	set uresid [lsort -dictionary $uresid]
	set amino_acids ""
	foreach rid $uresid {
		set sel_rid [atomselect $molid "resid $rid"]
	    set rsnm [lsort -unique [$sel_rid get resname]]
	    #puts $rsnm
	    set aa_in [lsearch -exact $aa_triple $rsnm]
	    append amino_acids [lindex $aa_single $aa_in]
	    $sel_rid delete
	}
	set cter [atomselect $molid "resid [lindex $uresid end]"]
	set cter_atms [$cter get name]
	if {[lsearch -exact $cter_atms "NT"] >= 0} {
    append amino_acids "-NH2"
	}
	puts "Amino Acid Sequence: $amino_acids"
	$sel delete
	return $amino_acids
}



proc getResnameResid { {molid top} {selection protein}} {
	set sel [atomselect $molid $selection]
	set uresid [lsort -unique [$sel get {resid resname}]]
	## sorting the prod_list according ASCII
	set uresid [lsort -dictionary $uresid]
	foreach rid $uresid {
	    puts "$rid"
	}
	$sel delete
}




proc assign_uniqueid {pdbfile} {
	mol new $pdbfile.pdb
	pbc set {0 0 0} 
	set proteins [atomselect top "protein"]
	set outfile [open ${pdbfile}_u.pdb w+]
	## pfrag does not work if protein fragments clash and mistacally gives -ve value
	set chains [lsort -unique [$proteins get pfrag]]
	foreach chain $chains {
		set sel [atomselect top "pfrag $chain"];
		$sel set segname P$chain
		$sel writepdb frag_temp.pdb
		set temp_data [split [read [open frag_temp.pdb]] "\n"]
		if {[lindex $temp_data 0 0] == "CRYST1"} {
			puts "Box info removed from pdb of pfrag $chain"
			puts $outfile [join [lrange $temp_data 1 end-1] "\n"]
		} else {
			puts "No box info in given pdb"
			puts $outfile [join [lrange $temp_data 0 end-1] "\n"]
		}
	}
	close $outfile
	mol delete top
	assign_ENDtoTER_pdb ${pdbfile}_u
}


proc assign_ENDtoTER_pdb {filename} {
    mol new $filename.pdb
    set myfile [open ${filename}.pdb]
    set lines [split [read $myfile] "\n"]
    close $myfile;
    set outfile [open ${filename}.pdb w+]
    foreach line [lrange $lines 0 end-3] { ;# loop upto end-3 to skip last "END"
        # do something with each line...
        if {$line == "END"} {
            puts $outfile "TER"
        } elseif {[string range $line 0 3] eq "ATOM"} {
            puts $outfile $line
        }
	}
    puts $outfile END
    close $outfile
    mol delete top
}


proc trajWrapSscacheCterrep {{syswrap No} {sysstr Yes}} {
	set sel [atomselect top "protein"]
	mol representation CPK
	mol color ColorID 4
	set cter_index [lindex [$sel get resid] end]; ### get resid of cter residue
	mol selection "protein and resid ${cter_index} and name CA"
	[atomselect top "protein and resid ${cter_index} and name CA"] set radius 7
	mol addrep top
	if {$syswrap !=  No} {
		pbc wrap -center com -centersel "protein" -compound fragment -all
	}

	if {$sysstr != No} {
		global pathVmdScript
		source $pathVmdScript/sscache.tcl
		start_sscache top
		animate forward
	}
	puts "done"
}




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


proc loadCifAF3Multijobs {{prefix ""} {refStr ${refStructures}/DCLK1_KD_tail.pdb} {KDResids "0 to 258"}} {
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
        TrajAlignTwoMol 0 top "backbone and resid ${KDResids}"
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


proc trajSave {{outfile traj.trr}} {
	set sel [atomselect top all]
	$sel writepdb frm0.pdb
	# Write trajectory
	animate write trr $outfile beg 0 end -1 skip 1 waitfor all top

}

proc transparencySet {new_alpha} {
        # This will always get the correct colors even if VMD
        # gains new definitions in the future
        set color_start [colorinfo num]
        set color_end [expr $color_start * 2]
        # Go through the list of colors (by index) and 
        # change their transp. value
        for {set color $color_start} {$color < $color_end} {incr color} {
                color change alpha $color $new_alpha
        }
}


# get the frames in the movie
proc snapshot {{molid top} {dim 2000} {filename ""}} {
	display depthcue off
	if {$filename eq ""} {
		set frm [molinfo $molid get frame]
		set f_name [lindex [split [molinfo $molid get name] "."] 0]
		set filename ${f_name}_${frm}.pov
	}
	puts $filename
	#render POV3 $filename  povray +W$dim +H$dim -I%s -O%s.tga +X +A +FT +UA +WT8
	render POV3 $filename \
    povray +W$dim +H$dim -I%s -O%s.tga \
    +A0.3 +R2 +UA +FT +WT12
	exec convert $filename.tga $filename.png
	exec rm -f $filename $filename.tga
	puts "redering finished"
}



proc snapframes {start step stop} {
    for {set i $start} {$i <= $stop} {incr i $step} {
        puts "Processing frame: $i"
        animate goto $i
        display update ui
        snapshot top 1000
    }
}



# change resid for each chain and then set sel [atomselect top all]; $sel writepdb name.pdb
proc changeResidOfChain {{chainId} {start_resid}} {
	set sel [atomselect top "chain $chainId"]
	set resid_list [lsort -integer [lsort -unique [$sel get resid]]]
	set new_resid $start_resid
	foreach old_resid $resid_list {
	    [atomselect top "chain $chainId and resid $old_resid"] set resid $new_resid
	    puts "$old_resid to $new_resid"
	    incr new_resid
	}
}


