#!/usr/bin/tclsh
#### Basic functions to apply on trajectory/frames in VMD @2022
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy

proc alignMultiMolid {{molid_refid 0} {molid_align top} {sel "protein and backbone"}} {
	### mol_align : the molids which are to be aligned
	set sel0 [atomselect $molid_refid $sel]
	foreach al $molid_align {
		set sel1 [atomselect $al $sel]
		set M [measure fit $sel1 $sel0]
		$sel1 move $M
	}
}

proc getSequence {{molid top} {selection "protein"} {NH2_check "No"}} {
	set aa_triple { ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR HSE} ; ### HSE alias for HIS
	set aa_single { A C D E F G H I K L M N P Q R S T V W Y H}; # last H for HSE
	set sel [atomselect $molid $selection]
	set uresid [lsort -unique [$sel get resid]]
	## sorting the prod_list according ASCII
	set uresid [lsort -dictionary $uresid]
	set amino_acids ""
	foreach rid $uresid {
		set sel_rid [atomselect $molid "resid $rid and protein"]
	    set rsnm [lsort -unique [$sel_rid get resname]]
	    #puts $rsnm
	    set aa_in [lsearch -exact $aa_triple $rsnm]
	    append amino_acids [lindex $aa_single $aa_in]
	    $sel_rid delete
	}
	$sel delete
	set output $amino_acids
	## check for NH2 at C-terminal
	if {$NH2_check != "No"} {
		set cter [atomselect $molid "resid [lindex $uresid end]"]
		set cter_atms [$cter get name]
		set is_NH2 "No"
		if {[lsearch -exact $cter_atms "NT"] >= 0} {
			append amino_acids "-NH2"
			set is_NH2 "Yes"
		}
		$cter delete
		set output [list $amino_acids $is_NH2]
	}
	return $output
}



proc getResnameResid { {molid top} {selection protein}} {
	set sel [atomselect $molid $selection]
	set uresid [lsort -unique [$sel get {resid resname}]]
	## sorting the prod_list according ASCII
	set uresid [lsort -dictionary $uresid]
	foreach rid $uresid {
	    puts "$rid"
	}
	# Print first and last elements
    set a [lindex $uresid 0]
    set b [lindex $uresid end]
    puts "First: $a"
    puts "Last:  $b"
	$sel delete
}



proc getResid { {molid top} {selection protein}} {
	set sel [atomselect $molid $selection]
	set uresid [lsort -unique [$sel get resid]]
	## sorting the prod_list according ASCII
	set uresid [lsort -dictionary $uresid]
	# Print first and last elements
    set a [lindex $uresid 0]
    set b [lindex $uresid end]
    puts "First: $a"
    puts "Last:  $b"
	$sel delete
	return "$a to $b"
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




proc setUniqueID {pdbfile} {
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

proc setENDtoTERinPdb {filename} {
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



proc changeTransparency {new_alpha} {
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


proc findClash { {molid top} {cutoff 1.5} {rep "No"} }  {
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


proc getCellDimention {{molid top}} {
	 set all [atomselect $molid all]
	 set minmax [measure minmax $all]
	 set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]]
	 set a "cellBasisVector1 [lindex $vec 0] 0 0"
	 set b "cellBasisVector2 0 [lindex $vec 1] 0"
	 set c "cellBasisVector3 0 0 [lindex $vec 2]"
	 set center "cellOrigin [measure center $all]"
	 $all delete
	 return [list $a\n$b\n$c\n$center]
}


proc getCharge {{sel "all"}} {
    set charge [eval "vecadd [[atomselect top $sel] get charge]"]
    return $charge
}



proc getResname { {molid top} {selection all} } {
    set sel [atomselect $molid $selection]
    set resnames [lsort -unique [$sel get resname]]
    $sel delete
    return $resnames
}


proc getSysInfo {{molid top} {write_file "No"}} {
	set aa_seqs "Amino acid Sequence: [lindex [getSequence $molid] 0]"
	if {[lindex [getSequence $molid] 1] == "Yes"} {
		set aa_num "Peptide length: [expr [llength [split [lindex [getSequence $molid] 0] ""]] - 4]"
	} else {
		set aa_num "Peptide length: [llength [split [lindex [getSequence $molid] 0] ""]]"
	}
	set cell_dimention "system dimention: [getCellDimention $molid]"
	set frame "Frame_num: [molinfo top get numframes]"
	set Schar "System charge: [eval "vecadd [[atomselect $molid all] get charge]"]"
	set Pchar "Protein charge: [eval "vecadd [[atomselect $molid "protein"] get charge]"]"
	set all "all_atoms_num: [[atomselect top all] num]"
	set wat "Water_num: [expr double ([[atomselect top "water"] num])/3]"
	set pocket_charge "Pocket_charge: [getCharge "protein and (within 5.0 of resname AMP ADP ATP GTP MG)"]"
	#set TFE "TFE: [expr double ([[atomselect top "resname TFE"] num])/9]"
	#set Wat_TFE "Water:TFE= [expr double ([expr [[atomselect top "water"] num]/3])/[expr [[atomselect top "resname TFE"] num]/9]]"
	set ion "Ions: [[atomselect top "ions"] num]"
 	set all_info [list $aa_seqs\n$aa_num\n$Pchar\n$cell_dimention\n$Schar\n$frame\n$all\n$wat\n$ion\n$pocket_charge]
 	if {$write_file != "No"} {
		set outfile [open system_info.dat w]
		puts $outfile $all_info
		puts $outfile "Total clash_atoms =  [findClash]"
		close $outfile
 	}
 	
 	return $all_info
}


proc generatePdbsIdentifiers {n} {
    set identifiers {}
    set count 0
    for {set first 0} {$first <= 9} {incr first} {
        for {set second 65} {$second <= 90} {incr second} {
            for {set third 65} {$third <= 90} {incr third} {
                for {set fourth 0} {$fourth <= 9} {incr fourth} {
                    set id "${first}[format %c $second][format %c $third]$fourth"
                    lappend identifiers $id
                    incr count
                    if {$count >= $n} {
                        return $identifiers
                    }
                }
            }
        }
    }
    return $identifiers
}

proc writePdbSel {{sel all} {outfile ""}} {
	set selection [atomselect top $sel]
	if {$outfile eq ""} {
		set frame [molinfo top get frame]
		set outfile "${frame}.pdb"
	}
	$selection writepdb $outfile
	$selection delete
	puts "PDB file written to $outfile"
}




proc writePdbWithActualIndex {{sel "protein"} {outfile} {molID top}} {
	set selection $sel
	set sel [atomselect $molID $sel]
	set keep_serials [$sel list]
	# Convert the list into a set for fast lookup
	array set keep_atoms {}
	foreach s $keep_serials {
		set s [expr {$s + 1}]
		set keep_atoms($s) 1
	}

	### write temp pdb for all atoms
	set selall [atomselect $molID "all"]
	$selall writepdb temp.pdb

	# Open files
	set fin [open temp.pdb r]
	set fout [open $outfile w]

	# Write PDB header
	#puts $fout "HEADER    Selected Atoms from Molecule ID $molID"
	puts $fout "TITLE     Selection: $selection"
	puts $fout "REMARK    Original indices preserved in this PDB file"

	# Loop through each line of the PDB
	while {[gets $fin line] >= 0} {
		if {[string match "ATOM*" $line] || [string match "HETATM*" $line]} {
			# Extract the serial number (columns 7-11, 1-based in PDB)
			set serial [string trim [string range $line 6 10]]
			if {[info exists keep_atoms($serial)]} {
				puts $fout $line
			}
		} elseif {[string match "TER*" $line] || [string match "END*" $line]} {
			puts $fout $line
		}
	}

	# Close files
	close $fin
	close $fout
	$sel delete
	$selall delete
	exec rm -f temp.pdb
	puts "PDB file with actual indices written to $outfile"
}


proc writeFasta {{molIDs top} {sel "protein"} {outfile ""}} {
	if {$outfile eq ""} {
		set outfile "sequence.fasta"
	}
	set fout [open $outfile w]
	foreach molID $molIDs {
		set name [molinfo $molID get name]
		puts $fout ">Molecule_$molID : $name"
		set seqs [getSequence $molID $sel]
		puts $fout $seqs
	}
	close $fout
	puts "FASTA file written to $outfile"
}


proc getDistBetTwoSelections {{molid top} {sel1 "name CA"} {sel2 "name CA"}} {
    # Create atom selections
    set selection1 [atomselect $molid $sel1]
    set selection2 [atomselect $molid $sel2]

    # Get positions
    set pos1 [$selection1 get {x y z}]
    set pos2 [$selection2 get {x y z}]

    # Compute geometric center for sel1
    set com1 {0.0 0.0 0.0}
    set n1 [llength $pos1]
    foreach p $pos1 {
        set com1 [vecadd $com1 $p]
    }
    set com1 [vecscale [expr {1.0 / $n1}] $com1]

    # Compute geometric center for sel2
    set com2 {0.0 0.0 0.0}
    set n2 [llength $pos2]
    foreach p $pos2 {
        set com2 [vecadd $com2 $p]
    }
    set com2 [vecscale [expr {1.0 / $n2}] $com2]

    # Compute distance between geometric centers
    set dist [veclength [vecsub $com1 $com2]]

    # Clean up
    $selection1 delete
    $selection2 delete

    return $dist
}


proc writeTwoSelectionsToOnePDB {molid1 seltext1 molid2 seltext2 outfile} {
    # Create selections
    set sel1 [atomselect $molid1 $seltext1]
    set sel2 [atomselect $molid2 $seltext2]

    # Temporary files to hold each selection
    set tmp1 "tmp_sel1.pdb"
    set tmp2 "tmp_sel2.pdb"

    # Write selections to temporary PDBs
    $sel1 writepdb $tmp1
    $sel2 writepdb $tmp2

    # Open output file
    set fp [open $outfile w]

    # Write selection 1
    set f1 [open $tmp1 r]
    while {[gets $f1 line] >= 0} {
        if {[string match "ATOM*" $line] || [string match "HETATM*" $line]} {
            puts $fp $line
        }
    }
    close $f1

    # Write selection 2
    set f2 [open $tmp2 r]
    while {[gets $f2 line] >= 0} {
        if {[string match "ATOM*" $line] || [string match "HETATM*" $line]} {
            puts $fp $line
        }
    }
    close $f2
    puts $fp "END"

    # Close merged file
    close $fp

    # Clean up temporary files
    file delete $tmp1 $tmp2
    $sel1 delete
    $sel2 delete

    puts "Merged PDB written to $outfile (MODEL 1: molid $molid1, MODEL 2: molid $molid2)"
}


