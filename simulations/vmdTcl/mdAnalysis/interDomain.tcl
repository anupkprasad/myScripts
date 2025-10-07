# Procedure to compute interdomain distance
# Arguments:
#   molid     → molecule id
#   sel1  → atom selection text for domain 1
#   sel2  → atom selection text for domain 2
#   outfile   → output filename

proc interDomainDistance {molid sel1 sel2 outfile} {
    # Open output file
    set out [open $outfile w]
    puts $out "# Frame   Distance(Å)"
    # Define selections
    set sel1 [atomselect $molid $sel1]
    set sel2 [atomselect $molid $sel2]

    # Loop over trajectory frames
    set nframes [molinfo $molid get numframes]
    for {set i 0} {$i < $nframes} {incr i} {
        animate goto $i
        $sel1 update
        $sel2 update

        # Compute COMs
        set com1 [measure center $sel1 weight mass]
        set com2 [measure center $sel2 weight mass]

        # Distance between COMs
        set dist [veclength [vecsub $com1 $com2]]
        puts $out "$i   $dist"
    }

    # Clean up
    close $out
    $sel1 delete
    $sel2 delete
    puts "Interdomain distance written to $outfile"
}



proc getInterface {molid selA selB cutoff outfile} {
    set a [atomselect $molid $selA]
    set b [atomselect $molid $selB]
    set coordsA [$a get index]
    set coordsB [$b get index]
    set out [open $outfile w]
    puts $out "# resA resB atomA atomB dist"
    set atomA [$a get index]
    set atomB [$b get index]
    foreach i $atomA {
        foreach j $atomB {
            set sel_pair [atomselect $molid "index $i or index $j"]
            set d [veclength [vecsub [lindex [$sel_pair get {x y z}] 0] [lindex [$sel_pair get {x y z}] 1]]]
            $sel_pair delete
            if {$d <= $cutoff} {
                set resA [lindex [$a get resid] [lsearch -exact $atomA $i]]
                set resB [lindex [$b get resid] [lsearch -exact $atomB $j]]
                puts $out "$resA $resB $i $j $d"
            }
        }
    }
    close $out
    $a delete; $b delete
}
# Example usage:
# get_interface 0 "segid A and resid 100 to 220 and name CA" "segid B and resid 400 to 520 and name CA" 4.5 "interface_pairs.txt"



##Draw principal axes on the structure

proc drawPrincipalAxes { {molid top} {sel "protein"}} {
    # Select your domain (example: kinase domain)
    set sel1 [atomselect top $sel]
    # Compute principal axes from inertia tensor
    set pa [lindex [measure inertia $sel1] 1]
    puts "Principal axes: $pa"

    # Get center of mass of selection
    set com [measure center $sel1 weight mass]

    # Draw 3 principal axes; red=1, green=2, blue=3
    foreach axis $pa color {red green blue} {
        set end [vecadd $com [vecscale 20 $axis]] ;# length=20 Å
        graphics top color $color
        graphics top line $com $end width 3
    }
}



proc interDomainAnalysis {molid selA selB {outfile ""}} {
    # molid   : molecule id
    # selA    : atom selection for domain A
    # selB    : atom selection for domain B
    # outfile : output filename (optional)

    set nframes [molinfo $molid get numframes]
    puts "Analyzing interdomain geometry for $nframes frames"

    # Open output file if specified
    set output_file ""
    if {$outfile ne ""} {
        set output_file [open $outfile w]
        puts $output_file "# Frame, COM_distance(Å), Angle(deg), Dihedral(deg)"
    }

    for {set f 0} {$f < $nframes} {incr f} {
        # Select domains for this frame
        set selAframe [atomselect $molid $selA frame $f]
        set selBframe [atomselect $molid $selB frame $f]

        # Centers of mass
        set comA [measure center $selAframe weight mass]
        set comB [measure center $selBframe weight mass]
        set distCOM [veclength [vecsub $comA $comB]]

        # Principal axes (take eigenvectors)
        set paA [lindex [lindex [measure inertia $selAframe] 1] 0]
        set paB [lindex [lindex [measure inertia $selBframe] 1] 0]

        # --- Angle between domain axes ---
        set dot [vecdot $paA $paB]
        set norm [expr {[veclength $paA] * [veclength $paB]}]
        set cosang [expr {$dot / $norm}]
        if {$cosang > 1} {set cosang 1}
        if {$cosang < -1} {set cosang -1}
        set angle [expr {acos($cosang) * 180 / acos(-1)}]

        # --- Dihedral (torsional) angle ---
        set hinge [vecsub $comB $comA]
        set projA [vecsub $paA [vecscale [expr {[vecdot $paA $hinge]/[veclength2 $hinge]}] $hinge]]
        set projB [vecsub $paB [vecscale [expr {[vecdot $paB $hinge]/[veclength2 $hinge]}] $hinge]]
        set normA [vecnorm $projA]
        set normB [vecnorm $projB]
        set cross [veccross $normA $normB]
        set sign [expr {[vecdot $cross $hinge] >= 0 ? 1 : -1}]
        set cosphi [vecdot $normA $normB]
        if {$cosphi > 1} {set cosphi 1}
        if {$cosphi < -1} {set cosphi -1}
        set dihedral [expr {$sign * acos($cosphi) * 180 / acos(-1)}]

        # Write to output file if specified
        if {$output_file ne ""} {
            puts $output_file [format "%d, %.3f, %.3f, %.3f" $f $distCOM $angle $dihedral]
        }

        $selAframe delete
        $selBframe delete
    }

    # Close output file
    if {$output_file ne ""} {
        close $output_file
        puts "Interdomain analysis results written to: $outfile"
    }
}



# Example usage for interdomain geometric analysis:
# interDomainAnalysis top "resid 54 to 311" "resid 437 to 623" "interdomain_geometry.dat"

