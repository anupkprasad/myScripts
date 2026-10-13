# and the first timestep for the given molecule id (default: top)
proc trajRmsd {{molid top} {selection "protein"}} {
        #frame 0 for the reference
        set reference [atomselect $molid $selection frame 0]
        # the frame being seld
        set sel [atomselect $molid $selection]
        set outfile [open trajRmsd.dat w]
        set frames [molinfo $molid get numframes]
        for {set frame 0} {$frame < $frames} {incr frame} {
                # get the correct frame
                $sel frame $frame
                # compute the transformation
                set trans_mat [measure fit $sel $reference]
                # do the alignment
                $sel move $trans_mat
                # compute the RMSD
                set rmsd [measure rmsd $sel $reference]
                # write in outfile
                puts $outfile "$frame $rmsd"
        }
        close $outfile
}



proc trajRmsdAllRes {{molid top} {selection "protein"}} {
        #frame 0 for the reference
        set reference [atomselect $molid $selection frame 0]
        # the frame being sel
        set sel [atomselect $molid $selection]
        set outfile [open trajRmsd.dat w]
        set frames [molinfo $molid get numframes]
        for {set frame 0} {$frame < $frames} {incr frame} {
                # get the correct frame
                $sel frame $frame
                # compute the transformation
                set trans_mat [measure fit $sel $reference]
                # do the alignment
                $sel move $trans_mat
                # compute the RMSD
                set rmsd [measure rmsd $sel $reference]
                # write in outfile
                puts $outfile "$frame $rmsd"
        }
        close $outfile
}

# RMSD between two molecules (molid_ref fixed at frame 0 vs molid_traj over all frames)
# Example: molRmsd 0 "protein and name CA" 1 "protein and name CA"
proc molRmsd {{molid_ref 0} {sel1 "protein and name CA"} {molid_traj 1} {sel2 "protein and name CA"} {outfile "molRmsd.dat"}} {
        set reference [atomselect $molid_ref $sel1 frame 0]
        set compare [atomselect $molid_traj $sel2]

        if {[$reference num] != [$compare num]} {
                puts "Error: atom count mismatch (ref [$reference num] vs traj [$compare num])"
                $reference delete
                $compare delete
                return
        }

        set out [open $outfile w]
        set frames [molinfo $molid_traj get numframes]
        for {set frame 0} {$frame < $frames} {incr frame} {
                $compare frame $frame
                set trans_mat [measure fit $compare $reference]
                $compare move $trans_mat
                set rmsd [measure rmsd $compare $reference]
                puts $out "$frame $rmsd"
        }
        close $out
        $reference delete
        $compare delete
        puts "Wrote $outfile ($frames frames)"
}
