proc mergePdbs {pdbList} {
	set outfile [open merge.pdb w+]
	foreach pdb $pdbfiles {
	  set myfile [open $pdb]
	  set data [read $myfile]
	  puts $outfile $data
	  close $myfile;
	}
	close $outfile
	assign_ENDtoTER_pdb prot_lip
}


proc protMoveOnLipSurface {{pdbfile prot_lip.pdb}} {
	mol new $pdbfile
	set prot [atomselect top "protein"]
	moveSysOrigin
	$prot moveby {0 0 33.34}  ; # Z of lip surface is 23.34 and +10 pad
	set sel [atomselect top all]
	$sel writepdb prot_lip.pdb
}

