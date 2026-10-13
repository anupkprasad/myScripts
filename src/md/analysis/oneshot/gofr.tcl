#loads a trajectory, then calcualtes and writes out a RDF
#does not write the integrated values correctly


#set up the atom selections
set sel1 [atomselect top "$sel_own"]
set sel2 [atomselect top "$sel_own2"]

#calculate g(r)
set gr [measure gofr $sel1 $sel2 delta .1 rmax 20 usepbc 1 selupdate 1 first $first_frame last $last_frame step 1]
# delta = bin value (default = .1)

#set up the outfile and write out the data

## r = radial distance, gr2 = g(r), igr = number integral for all pairs of atoms in the two selections.
set outfile [open gofr_${sel_own}${sel_own2}${first_frame}to${last_frame}.dat w]

set r [lindex $gr 0]
set gr2 [lindex $gr 1]
set igr [lindex $gr 2]

set i 0
foreach j $r k $gr2 l $igr {
   puts $outfile "$j $k $l"
}

close $outfile
