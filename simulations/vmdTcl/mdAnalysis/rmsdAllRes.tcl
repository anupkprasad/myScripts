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
