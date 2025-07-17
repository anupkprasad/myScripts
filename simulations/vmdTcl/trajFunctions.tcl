#!/usr/bin/tclsh
#### Basic functions to apply on trajectory/frames in VMD @2022
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy

######## Functions for Trajectories  ############

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


proc trajAlignTwoMol {{molidRef top} {molid top} {selRef "protein and backbone"} {sel ""}} {
	## sel: segment (Kinase domain) which are required to align, allsel: transformation matrix to entire molecules (Kinase domaain +tail)
	set sel0 [atomselect $molidRef $selRef]
	puts "Reference atoms: [$sel0 num]"
	set n [molinfo $molid get numframes]
	for { set i 0 } { $i < $n } { incr i } {
		if {$sel eq ""} {
			set sel $selRef
		}

		set sel1 [atomselect $molid $sel]
		if {$i == 0} {
			puts "Selected atoms: [$sel1 num]"
		}
		set allsel [atomselect $molid all frame $i]
		$sel1 frame $i
		set M [measure fit $sel1 $sel0]
		$allsel move $M
	}
}



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
		global pathVmdScript
		source $pathVmdScript/simulations/vmdTcl/sscache.tcl
		start_sscache top
		animate forward
	}
	puts "done"
}

proc trajSave {{outfile traj.trr} {sel "protein"}} {
	set sel [atomselect top $sel]
	$sel frame 0
	$sel writepdb frm0.pdb
	# Write trajectory
	animate write trr $outfile beg 0 end -1 skip 1 waitfor all top
	puts "Trajectory saved to $outfile and first frame saved to frm0.pdb"
}
