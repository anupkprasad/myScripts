#!/usr/bin/tclsh
#### Basic functions to apply on trajectory/frames in VMD @2022
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy



######## Functions for Alignment  ############


proc trajAlign {{sel "protein"} {Refframe 0}} {
	set molid top
	## sel: segment (Kinase domain) which are required to align, allsel: transformation matrix to entire molecules (Kinase domaain +tail)
	set sel0 [atomselect $molid $sel]
	$sel0 frame $Refframe
	set n [molinfo $molid get numframes]
	for { set i 0 } { $i < $n } { incr i } {
		set sel1 [atomselect $molid $sel]
		#set allsel [atomselect $molid "protein" frame $i]
		set allsel [atomselect $molid all frame $i]
		$sel1 frame $i
		set M [measure fit $sel1 $sel0]
		$allsel move $M
	}
}


# Helper function to adjust residue ranges to have equal number of atoms
proc adjustRangeForAlignment {mol1 range1 mol2 range2} {
    # Parse ranges - handle both "X-Y" and "X to Y" formats
    if {[string match "*to*" $range1]} {
        set r1 [split $range1]
        set start1 [lindex $r1 0]
        set end1 [lindex $r1 2]
    } else {
        set r1 [split $range1 "-"]
        set start1 [lindex $r1 0]
        set end1 [lindex $r1 1]
    }
    
    if {[string match "*to*" $range2]} {
        set r2 [split $range2]
        set start2 [lindex $r2 0]
        set end2 [lindex $r2 2]
    } else {
        set r2 [split $range2 "-"]
        set start2 [lindex $r2 0]
        set end2 [lindex $r2 1]
    }

    ### Check number of CA atoms in both ranges
    set sel1 [atomselect $mol1 "resid $start1 to $end1 and name CA"]
    set sel2 [atomselect $mol2 "resid $start2 to $end2 and name CA"]
    
    set num1 [$sel1 num]
    set num2 [$sel2 num]
    
    puts "Range comparison: $range1 ($num1 atoms) vs $range2 ($num2 atoms)"
    
    $sel1 delete
    $sel2 delete
    
    if { $num1 == $num2 } {
        return [list "$start1 to $end1" "$start2 to $end2"]
    }
    
    # Use minimum length
    set len1 [expr $end1 - $start1 + 1]
    set len2 [expr $end2 - $start2 + 1]
    set min_len [expr min($len1, $len2)]
    
    set adj_range1 "$start1 to [expr $start1 + $min_len - 1]"
    set adj_range2 "$start2 to [expr $start2 + $min_len - 1]"

    puts "Adjusted ranges: $adj_range1, $adj_range2"
    return [list $adj_range1 $adj_range2]
}



proc trajAlignTwoMol {{molidRef top} {molid top} {selRef "protein and backbone"} {sel ""}} {
    ## it aligns each frame of molid to first frame of molidRef based on selRef atoms
    
    # selRef: reference selection (on molidRef)
    # sel: mobile selection to be aligned (default = selRef)
    # molidRef: reference molid (usually static)
    # molid: trajectory molid (to be aligned)

    # Reference selection fixed at frame 0
    set sel0 [atomselect $molidRef $selRef frame 0]
    puts "Reference atoms: [$sel0 num]"

    if {$sel eq ""} {
        set sel $selRef
    }

    set n [molinfo $molid get numframes]

    for {set i 0} {$i < $n} {incr i} {
        set sel1 [atomselect $molid $sel frame $i]
        if {$i == 0} {
            puts "Selected atoms: [$sel1 num]"
        }

        set allsel [atomselect $molid all frame $i]
        set M [measure fit $sel1 $sel0]
        $allsel move $M

        # cleanup per frame
        $sel1 delete
        $allsel delete
    }
    $sel0 delete
}



proc alignTwoMol {{molidRef top} {molid top} {selRef "protein and backbone"} {sel ""}} {
    ## it aligns first frame of molid to first frame of molidRef based on selRef atoms, then aligns all other frames of molid to its first frame based on sel atoms.
    # selRef: reference selection (on molidRef)
    # sel: mobile selection to be aligned (default = selRef)
    # molidRef: reference molid (usually static)
    # molid: trajectory molid (to be aligned)

    # Reference selection fixed at frame 0
    set sel0 [atomselect $molidRef $selRef frame 0]
    puts "Reference atoms: [$sel0 num]"

    if {$sel eq ""} {
        set sel $selRef
    }

    set n [molinfo $molid get numframes]

    for {set i 0} {$i < $n} {incr i} {
        set allsel [atomselect $molid all frame $i]
        set sel1 [atomselect $molid $sel frame $i]
        if {$i == 0} {
            set sel1_f0 [atomselect $molid $sel frame $i]
            puts "Selected atoms: [$sel1 num]"
            set M [measure fit $sel1_f0 $sel0]
            $allsel move $M
            continue
        }
        
        set M [measure fit $sel1 $sel1_f0]
        $allsel move $M

        # cleanup per frame
        $sel1 delete
        $allsel delete
    }
    $sel0 delete
    $sel1_f0 delete
}


proc trajAlignTwoMolByAlignment {{molidRef top} {seqNameRef} {molid top} {seqName} {alignmentFile}} {
    # Read alignment file and align structures based on sequence alignment
    # molidRef: reference molecule ID
    # seqNameRef: sequence name in alignment file for reference structure
    # molid: mobile molecule ID to be aligned
    # seqName: sequence name in alignment file for mobile structure
    # alignmentFile: path to alignment file (.xlsx exported as tab-delimited or CSV)
    
    puts "Reading alignment file: $alignmentFile"
    
    # Read the alignment file
    set fp [open $alignmentFile r]
    set content [read $fp]
    close $fp
    
    # Handle different file formats (tab-delimited, CSV, or space-separated)
    set lines [split $content "\n"]
    
    # Parse alignment data - find lines starting with sequence names
    set refLine ""
    set mobLine ""
    set refResids {}
    set mobResids {}
    
    foreach line $lines {
        set line [string trim $line]
        if {$line eq ""} continue
        
        # Try different delimiters (tab, comma, or multiple spaces)
        set lineData {}
        if {[string match "*\t*" $line]} {
            set lineData [split $line "\t"]
        } elseif {[string match "*,*" $line]} {
            set lineData [split $line ","]
        } else {
            # Handle space-separated values
            set lineData [regexp -all -inline {\S+} $line]
        }
        
        if {[llength $lineData] < 2} continue
        
        set seqName [lindex $lineData 0]
        set residues [lrange $lineData 1 end]
        
        # Check if this line contains the reference sequence
        if {[string match "*$seqNameRef*" $seqName] || $seqName eq $seqNameRef} {
            set refLine $line
            set refResids $residues
            puts "Found reference sequence: $seqNameRef ([llength $residues] residues)"
        }
        
        # Check if this line contains the mobile sequence
        if {[string match "*$seqName*" $seqName] || $seqName eq $seqName} {
            set mobLine $line
            set mobResids $residues
            puts "Found mobile sequence: $seqName ([llength $residues] residues)"
        }
    }
    
    if {$refLine eq "" || $mobLine eq ""} {
        puts "Available sequences in file:"
        foreach line $lines {
            set line [string trim $line]
            if {$line eq ""} continue
            
            # Extract first column (sequence name)
            if {[string match "*\t*" $line]} {
                set seqName [lindex [split $line "\t"] 0]
            } elseif {[string match "*,*" $line]} {
                set seqName [lindex [split $line ","] 0]
            } else {
                set seqName [lindex [regexp -all -inline {\S+} $line] 0]
            }
            puts "  $seqName"
        }
        error "Could not find both sequences ($seqNameRef and $seqName) in alignment file"
    }
    
    puts "Reference residues: [llength $refResids] entries"
    puts "Mobile residues: [llength $mobResids] entries"
    
    # Ensure both lists have the same length
    set minLength [expr min([llength $refResids], [llength $mobResids])]
    if {[llength $refResids] != [llength $mobResids]} {
        puts "Warning: Residue lists have different lengths. Using minimum: $minLength"
        set refResids [lrange $refResids 0 [expr $minLength - 1]]
        set mobResids [lrange $mobResids 0 [expr $minLength - 1]]
    }
    
    # Find positions where both sequences have valid (non-NaN) residues
    set refAlignedResids {}
    set mobAlignedResids {}
    
    for {set i 0} {$i < $minLength} {incr i} {
        set refRes [lindex $refResids $i]
        set mobRes [lindex $mobResids $i]
        
        # Clean up residue values (remove any quotes or extra whitespace)
        set refRes [string trim $refRes "\"' "]
        set mobRes [string trim $mobRes "\"' "]
        
        # Only include positions where both sequences have valid (non-NaN) residues
        if {$refRes ne "NaN" && $mobRes ne "NaN" && $refRes ne "" && $mobRes ne "" && [string is integer $refRes] && [string is integer $mobRes]} {
            lappend refAlignedResids $refRes
            lappend mobAlignedResids $mobRes
        }
    }
    
    puts "Found [llength $refAlignedResids] aligned residue pairs"
    
    if {[llength $refAlignedResids] == 0} {
        error "No valid aligned residues found"
    }
    
    # Create selection strings
    set refSelStr "resid [join $refAlignedResids " "] and name CA"
    set mobSelStr "resid [join $mobAlignedResids " "] and name CA"
    
    puts "Reference selection: $refSelStr"
    puts "Mobile selection: $mobSelStr"
    
    # Create selections
    set refSel [atomselect $molidRef $refSelStr frame 0]
    set refAtoms [$refSel num]
    puts "Reference selection atoms: $refAtoms"
    
    if {$refAtoms == 0} {
        $refSel delete
        error "No atoms found for reference selection. Check residue numbering."
    }
    
    # Get number of frames in mobile molecule
    set numFrames [molinfo $molid get numframes]
    puts "Aligning $numFrames frames"
    
    # Align each frame
    for {set frame 0} {$frame < $numFrames} {incr frame} {
        set mobSel [atomselect $molid $mobSelStr frame $frame]
        set mobAtoms [$mobSel num]
        
        if {$frame == 0} {
            puts "Mobile selection atoms: $mobAtoms"
            if {$mobAtoms == 0} {
                $mobSel delete
                $refSel delete
                error "No atoms found for mobile selection. Check residue numbering."
            }
            if {$refAtoms != $mobAtoms} {
                puts "Warning: Atom count mismatch: ref=$refAtoms, mobile=$mobAtoms"
                puts "This may indicate missing residues or different numbering schemes."
            }
        }
        
        # Get all atoms for transformation
        set allSel [atomselect $molid "all" frame $frame]
        
        # Calculate transformation matrix
        if {$mobAtoms > 0 && $refAtoms > 0} {
            set transMatrix [measure fit $mobSel $refSel]
            # Apply transformation to all atoms
            $allSel move $transMatrix
        }
        
        # Clean up selections for this frame
        $mobSel delete
        $allSel delete
        
        if {$frame % 10 == 0 || $frame == [expr $numFrames - 1]} {
            puts "Aligned frame $frame/[expr $numFrames - 1]"
        }
    }
    
    # Clean up reference selection
    $refSel delete
    
    puts "Trajectory alignment completed based on sequence alignment"
    puts "Used [llength $refAlignedResids] CA atoms for alignment"
}





##################### Functions for visualization ############

proc trajMoveMol {{molid top } {sel all} {disXYZ {0 0 0}}} {
	## sel: segment (Kinase domain) which are required to align, allsel: transformation matrix to entire molecules (Kinase domaain +tail)
	set sel0 [atomselect $molid $sel]
	set n [molinfo $molid get numframes]
	for { set i 0 } { $i < $n } { incr i } {
		puts "Moving frame $i by $disXYZ"
		$sel0 frame $i
		$sel0 moveby $disXYZ
	}
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
		global pathAnalysis
		source [file join $pathAnalysis sscache.tcl]
		start_sscache top
		animate forward
	}
	puts "done"
}


####################### writing trajectory #######################
proc trajWritepdb {{sel "protein"} {frames ""} {outfiles ""}} {
	## frames: list of frames to write (ni: first frame, nf: last frame)
	if {$frames eq ""} {
		set ni 0
		set nf [molinfo top get numframes]
	} else {
		set ni [lindex $frames 0]
		set nf [lindex $frames end]
	}
	set framesW {}
	for {set i $ni} {$i < $nf} {incr i} {
		lappend framesW $i
	}


	if {$outfiles eq ""} {
		set outfiles $framesW
	}
	set sel [atomselect top $sel]
	for { set i 0 } { $i < [llength $framesW] } { incr i } {
		set frameIndex [lindex $framesW $i]
		$sel frame $frameIndex
		puts "Writing frame $frameIndex to file: [lindex $outfiles $i].pdb"
		set outfile [lindex $outfiles $i]
		$sel writepdb $outfile.pdb
	}
}	




proc trajSave {{outfile traj.trr} {sel "protein"}} {
	set sel [atomselect top $sel]
	$sel frame 0
	$sel writepdb frm0.pdb
	# Write trajectory
	animate write trr $outfile beg 0 end -1 skip 1 waitfor all top
	puts "Trajectory saved to $outfile and first frame saved to frm0.pdb"
}
