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


proc calRmsd {{molid top} {sel "protein and backbone"} {atomselectRef ""}} {
	if {$atomselectRef eq ""} {
                set reference [atomselect $molid $sel frame 0]
        } else {
                set reference $atomselectRef
        }
	set compare [atomselect $molid $sel]
        puts "RMSD of Reference [$reference num] and selection: [$compare num]"
	set outfile [open rmsd.dat w]

                set num_steps [molinfo $molid get numframes]
                for {set frame 0} {$frame < $num_steps} {incr frame} {
                        # get the correct frame
                        $compare frame $frame
                        # compute the transformation
                        set trans_mat [measure fit $compare $reference]
                        # do the alignment
                        $compare move $trans_mat
                        # compute the RMSD
                        set rmsd [measure rmsd $compare $reference]
                        # print the RMSD
                        #puts "RMSD of $frame is $rmsd"
                        puts $outfile "$frame $rmsd"
                }
                close $outfile
}


