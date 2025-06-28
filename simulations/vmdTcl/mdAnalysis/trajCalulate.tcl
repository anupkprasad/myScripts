# get the secondary structure of protein CA atoms, each line represent frame
proc calSecStr {{molid top} {CA_sel "protein and name CA"}} {
	set sel [atomselect $molid $CA_sel]
	set n [molinfo $molid get numframes]
	set outfile [open ss_traj.dat w]

	for { set i 0 } { $i < $n } { incr i } {
	 animate goto $i
	 display update off
	 $sel frame $i
	 mol ssrecalc $molid
	 set secstruct [$sel get structure]
	 puts $outfile "$secstruct"
	}
	display update on
	close $outfile
}