#!/usr/bin/tclsh
#### Basic functions to apply on trajectory/frames in VMD @2022
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy

######## Functions for Trajectories  ############

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


proc trajAlignTwoMol {{molidRef top} {molid top} {selRef "protein and backbone"} {sel "protein and backbone"}} {
	## sel: segment (Kinase domain) which are required to align, allsel: transformation matrix to entire molecules (Kinase domaain +tail)
	set sel0 [atomselect $molidRef $selRef]
	puts "Reference atoms: [$sel0 num]"
	set n [molinfo $molid get numframes]
	for { set i 0 } { $i < $n } { incr i } {
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


proc trajSave {{outfile traj.trr}} {
	set sel [atomselect top all]
	$sel writepdb frm0.pdb
	# Write trajectory
	animate write trr $outfile beg 0 end -1 skip 1 waitfor all top

}


proc trajSnapshot {start step stop} {
    for {set i $start} {$i <= $stop} {incr i $step} {
        puts "Processing frame: $i"
        animate goto $i
        display update ui
        snapshot top 1000
    }
}


