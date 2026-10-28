# interDomain.tcl
# Two-domain analysis on one VMD molecule, over a trajectory.
#
# Load in the VMD Tk console:
#   source ~/myScripts/src/md/analysis/interDomain.tcl
#   InterDomain help
#
# Output files: lines starting with # are metadata. The next line is the
# column header. Remaining lines are whitespace-separated values.
#   pandas.read_csv(path, sep=r"\s+", comment="#")
#
# Clustering concept:
#   geometcluster — k-means on geometry numbers (COM distance, axis angle,
#     twist/dihedral). Groups frames with similar domain–domain geometry,
#     not by atomic RMSD. Use when you care about domain separation,
#     orientation, or twist states.
#   rmsdcluster — average-linkage hierarchical clustering on RMSD of sel2
#     after fitting sel1. Groups frames by how similar domain 2 looks once
#     overall motion of domain 1 is removed. Use when you care about the
#     structural pose of one domain relative to the other.

proc InterDomain {args} {
    if {[llength $args] == 0} {
        _id_help
        return
    }
    set cmd [lindex $args 0]
    set rest [lrange $args 1 end]
    if {$cmd in {help --help -help -h}} {
        _id_help
        return
    }
    switch -exact -- $cmd {
        distance  { _id_cmd_distance $rest }
        com       { _id_cmd_com $rest }
        angle     { _id_cmd_angle $rest }
        dihedral  { _id_cmd_dihedral $rest }
        geometry  { _id_cmd_geometry $rest }
        interface { _id_cmd_interface $rest }
        geometcluster { _id_cmd_geometcluster $rest }
        rmsdcluster { _id_cmd_rmsdcluster $rest }
        draw      { _id_cmd_draw $rest }
        default {
            if {$cmd eq "cluster"} {
                puts "InterDomain: 'cluster' was renamed to 'geometcluster'"
            } else {
                puts "InterDomain: unknown command '$cmd'"
            }
            _id_help
        }
    }
}

proc _id_help {} {
    puts {
InterDomain — distance, COM, principal-axis angle, and interface of two selections

Usage:
  InterDomain <command> [options]
  InterDomain help

Commands:
  distance     COM–COM distance for each frame
  com          Center of mass of sel1 and sel2 for each frame
  angle        Angle between the first principal axes (0–90 deg)
  dihedral     Signed twist of those axes about the COM–COM vector
  geometry     distance, both COMs, angle, and dihedral in one pass
  interface    Residue pairs within a cutoff (default 5 A)
  geometcluster Group frames by selected geometry features, write PDBs
  rmsdcluster  Fit sel1, cluster frames by RMSD of sel2, write PDBs
  draw         Draw principal axes of the current frame in the VMD window
  help         Show this message

Common options (trajectory commands):
  -molid <id|top>   molecule id (default: top)
  -sel1 <text>      VMD selection for domain 1 (required)
  -sel2 <text>      VMD selection for domain 2 (required)
  -out <file>       output file (each command has its own default name)
  -first <frame>    first frame, inclusive (default: 0)
  -last <frame|end> last frame, inclusive (default: end)
  -stride <n>       use every n-th frame (default: 1)

interface also accepts:
  -cutoff <A>       contact distance in Angstrom (default: 5.0)

draw also accepts:
  -length <A>       axis length drawn from the COM (default: 20)
  -sel2 <text>      optional second domain

geometcluster also accepts:
  -features <list>  geometry features to cluster on (default: com,dihedral)
                    Names: com|distance, angle, dihedral
                    Examples: "com,angle"  "angle,dihedral"  "com,angle,dihedral"
  -ncluster <k>     number of groups (default: 3)
  -outdir <dir>     output directory (default: geomet_clusters)
  -data <file>      existing geometry file; if omitted, InterDomain_geometry.dat
                    is calculated automatically and then used
  -maxiter <n>      k-means iterations (default: 50)

rmsdcluster also accepts one of:
  -ncluster <k>     cut the tree into k groups
  -cutoff <A>       cut the tree where the average RMSD exceeds this value
  -outdir <dir>     output directory (default: rmsd_clusters)
  sel1 is fitted. sel2 is compared. Writes an N x N RMSD matrix, then
  average-linkage hierarchical clustering.

Selections must be quoted if they contain spaces.
Both selections are in the same molecule. Add "and noh" to ignore hydrogens.
Angle uses the first principal axis from "measure inertia".
The reported angle is the smaller one (0–90), so an axis flip does not jump the curve.
Dihedral sign is fixed to the first frame, then axes are kept continuous.

Python:
  df = pandas.read_csv("InterDomain_distance.dat", sep=r"\s+", comment="#")

Examples:
  InterDomain distance -molid 0 -sel1 "resid 54 to 311" -sel2 "resid 437 to 623" -out dist.dat
  InterDomain com -sel1 "resid 54 to 311 and name CA" -sel2 "resid 437 to 623 and name CA"
  InterDomain angle -sel1 "resid 54 to 311" -sel2 "resid 437 to 623" -first 0 -last 500 -stride 10
  InterDomain geometry -sel1 "resid 54 to 311" -sel2 "resid 437 to 623" -out geometry.dat
  InterDomain interface -sel1 "resid 54 to 311" -sel2 "resid 437 to 623" -cutoff 5 -out interface.dat
  InterDomain draw -sel1 "resid 54 to 311" -sel2 "resid 437 to 623" -length 25
  InterDomain geometcluster -sel1 "resid 54 to 311" -sel2 "resid 437 to 623" -features "com,angle" -ncluster 4
  InterDomain geometcluster -data InterDomain_geometry.dat -features "angle,dihedral" -ncluster 3
  InterDomain rmsdcluster -sel1 "resid 54 to 311" -sel2 "resid 437 to 623" -ncluster 5 -outdir rmsd_clusters

interface writes two files:
  <out>              residue pairs, occupancy, and minimum distance
  <out>_frames.dat   contact counts per frame (if <out> ends in .dat, the
                     stem is used: interface.dat -> interface_frames.dat)
}
}

proc _id_args {argv allowed defaults} {
    array set opt $defaults
    set n [llength $argv]
    for {set i 0} {$i < $n} {incr i} {
        set tok [lindex $argv $i]
        if {![string match -* $tok]} {
            error "InterDomain: unexpected '$tok'. Use -flag value. Type: InterDomain help"
        }
        set name [string range $tok 1 end]
        if {$name in {h help -help}} {
            _id_help
            set opt(_help) 1
            return [array get opt]
        }
        if {[lsearch -exact $allowed $name] < 0} {
            error "InterDomain: unknown option -$name. Allowed: -[join $allowed " -"]"
        }
        incr i
        if {$i >= $n} {
            error "InterDomain: -$name needs a value"
        }
        set opt($name) [lindex $argv $i]
    }
    return [array get opt]
}

proc _id_need {optName args} {
    upvar $optName opt
    foreach key $args {
        if {![info exists opt($key)] || $opt($key) eq ""} {
            error "InterDomain: missing -$key. Type: InterDomain help"
        }
    }
}

proc _id_frames {molid first last stride} {
    set nframes [molinfo $molid get numframes]
    if {$nframes < 1} {
        error "InterDomain: molecule $molid has no frames"
    }
    if {$last eq "end" || $last eq ""} {
        set last [expr {$nframes - 1}]
    }
    if {![string is integer -strict $first] || ![string is integer -strict $last] || ![string is integer -strict $stride]} {
        error "InterDomain: -first, -last, and -stride must be integers (-last may be end)"
    }
    if {$first < 0 || $last >= $nframes || $first > $last || $stride < 1} {
        error "InterDomain: bad range first=$first last=$last stride=$stride (nframes=$nframes)"
    }
    set nanalyzed [expr {($last - $first) / $stride + 1}]
    return [list $first $last $stride $nframes $nanalyzed]
}

proc _id_open_pair {molid sel1text sel2text minatoms} {
    # VMD deletes an atomselect when the procedure that created it returns.
    # Create sel1/sel2 in the caller so they stay valid for the frame loop.
    upvar 1 sel1 sel1
    upvar 1 sel2 sel2
    set sel1 [uplevel 1 [list atomselect $molid $sel1text]]
    set sel2 [uplevel 1 [list atomselect $molid $sel2text]]
    set n1 [$sel1 num]
    set n2 [$sel2 num]
    if {$n1 < $minatoms || $n2 < $minatoms} {
        error "InterDomain: selection too small (sel1=$n1 atoms, sel2=$n2 atoms, need >= $minatoms)"
    }
    set both [atomselect $molid "($sel1text) and ($sel2text)"]
    set noverlap [$both num]
    $both delete
    if {$noverlap > 0} {
        puts "InterDomain: warning: sel1 and sel2 share $noverlap atoms"
    }
}

proc _id_header {fh title meta} {
    puts $fh "# $title"
    foreach {key value} $meta {
        puts $fh "# $key: $value"
    }
    puts $fh "# pandas: read_csv(path, sep=r'\\s+', comment='#')"
}

proc _id_fmt {x} {
    if {$x eq "nan"} {
        return "nan"
    }
    return [format "%.4f" $x]
}

proc _id_field {x} {
    if {$x eq ""} {
        return "-"
    }
    return [string map {" " _ "\t" _} $x]
}

proc _id_progress {step total frame} {
    if {$step == 1 || $step == $total || $step % 100 == 0} {
        puts "InterDomain: frame $frame ($step/$total)"
    }
}

proc _id_first_axis {sel} {
    set axes [lindex [measure inertia $sel] 1]
    set axis [lindex $axes 0]
    set norm [veclength $axis]
    if {$norm < 1e-8} {
        error "InterDomain: principal axis is undefined; the selection may be too small or linear"
    }
    return [vecscale [expr {1.0 / $norm}] $axis]
}

proc _id_continuous {axis prevName} {
    upvar $prevName prev
    if {$prev ne "" && [vecdot $axis $prev] < 0.0} {
        set axis [vecscale -1.0 $axis]
    }
    set prev $axis
    return $axis
}

proc _id_axis_angle {axis1 axis2} {
    set cosine [expr {abs([vecdot $axis1 $axis2])}]
    if {$cosine > 1.0} { set cosine 1.0 }
    return [expr {acos($cosine) * 180.0 / acos(-1.0)}]
}

proc _id_dihedral {axis1 axis2 com1 com2} {
    set hinge [vecsub $com2 $com1]
    set h2 [veclength2 $hinge]
    if {$h2 < 1e-8} {
        return "nan"
    }
    set proj1 [vecsub $axis1 [vecscale [expr {[vecdot $axis1 $hinge] / $h2}] $hinge]]
    set proj2 [vecsub $axis2 [vecscale [expr {[vecdot $axis2 $hinge] / $h2}] $hinge]]
    if {[veclength2 $proj1] < 1e-10 || [veclength2 $proj2] < 1e-10} {
        return "nan"
    }
    set n1 [vecnorm $proj1]
    set n2 [vecnorm $proj2]
    set sign [expr {[vecdot [veccross $n1 $n2] $hinge] >= 0.0 ? 1.0 : -1.0}]
    set cosine [vecdot $n1 $n2]
    if {$cosine > 1.0} { set cosine 1.0 }
    if {$cosine < -1.0} { set cosine -1.0 }
    return [expr {$sign * acos($cosine) * 180.0 / acos(-1.0)}]
}

proc _id_frame_step {sel1 sel2 frame} {
    $sel1 frame $frame
    $sel2 frame $frame
    $sel1 update
    $sel2 update
}

proc _id_com_pair {sel1 sel2} {
    set com1 [measure center $sel1 weight mass]
    set com2 [measure center $sel2 weight mass]
    set dist [veclength [vecsub $com1 $com2]]
    return [list $com1 $com2 $dist]
}

proc _id_axes_pair {sel1 sel2 prev1Name prev2Name} {
    upvar $prev1Name prev1 $prev2Name prev2
    set axis1 [_id_continuous [_id_first_axis $sel1] prev1]
    set axis2 [_id_continuous [_id_first_axis $sel2] prev2]
    return [list $axis1 $axis2]
}

proc _id_common {argv extraDefaults extraAllowed} {
    set allowed [concat {molid sel1 sel2 out first last stride} $extraAllowed]
    set defaults [concat {
        molid top
        sel1 {}
        sel2 {}
        first 0
        last end
        stride 1
    } $extraDefaults]
    array set opt [_id_args $argv $allowed $defaults]
    if {[info exists opt(_help)]} {
        return -code return
    }
    _id_need opt sel1 sel2 out
    return [array get opt]
}

proc _id_cmd_distance {argv} {
    array set opt [_id_common $argv {out InterDomain_distance.dat} {}]
    lassign [_id_frames $opt(molid) $opt(first) $opt(last) $opt(stride)] first last stride nframes nanalyzed
    _id_open_pair $opt(molid) $opt(sel1) $opt(sel2) 1
    set fh [open $opt(out) w]
    _id_header $fh "InterDomain distance" [list \
        molid $opt(molid) sel1 $opt(sel1) sel2 $opt(sel2) \
        frames $nframes frames_analyzed $nanalyzed \
        first $first last $last stride $stride units Angstrom]
    puts $fh "frame distance_A"
    set step 0
    for {set f $first} {$f <= $last} {incr f $stride} {
        incr step
        _id_progress $step $nanalyzed $f
        _id_frame_step $sel1 $sel2 $f
        lassign [_id_com_pair $sel1 $sel2] com1 com2 dist
        puts $fh "$f [_id_fmt $dist]"
    }
    close $fh
    $sel1 delete
    $sel2 delete
    puts "InterDomain: wrote $opt(out)"
    puts "InterDomain: columns: frame distance_A"
}

proc _id_cmd_com {argv} {
    array set opt [_id_common $argv {out InterDomain_com.dat} {}]
    lassign [_id_frames $opt(molid) $opt(first) $opt(last) $opt(stride)] first last stride nframes nanalyzed
    _id_open_pair $opt(molid) $opt(sel1) $opt(sel2) 1
    set fh [open $opt(out) w]
    _id_header $fh "InterDomain COM" [list \
        molid $opt(molid) sel1 $opt(sel1) sel2 $opt(sel2) \
        frames $nframes frames_analyzed $nanalyzed \
        first $first last $last stride $stride units Angstrom]
    puts $fh "frame com1_x com1_y com1_z com2_x com2_y com2_z"
    set step 0
    for {set f $first} {$f <= $last} {incr f $stride} {
        incr step
        _id_progress $step $nanalyzed $f
        _id_frame_step $sel1 $sel2 $f
        lassign [_id_com_pair $sel1 $sel2] com1 com2 dist
        puts $fh "$f [_id_fmt [lindex $com1 0]] [_id_fmt [lindex $com1 1]] [_id_fmt [lindex $com1 2]] [_id_fmt [lindex $com2 0]] [_id_fmt [lindex $com2 1]] [_id_fmt [lindex $com2 2]]"
    }
    close $fh
    $sel1 delete
    $sel2 delete
    puts "InterDomain: wrote $opt(out)"
    puts "InterDomain: columns: frame com1_x com1_y com1_z com2_x com2_y com2_z"
}

proc _id_cmd_angle {argv} {
    array set opt [_id_common $argv {out InterDomain_angle.dat} {}]
    lassign [_id_frames $opt(molid) $opt(first) $opt(last) $opt(stride)] first last stride nframes nanalyzed
    _id_open_pair $opt(molid) $opt(sel1) $opt(sel2) 3
    set fh [open $opt(out) w]
    _id_header $fh "InterDomain principal-axis angle" [list \
        molid $opt(molid) sel1 $opt(sel1) sel2 $opt(sel2) \
        frames $nframes frames_analyzed $nanalyzed \
        first $first last $last stride $stride \
        axis "first principal axis of measure inertia" \
        angle "smaller angle between axes, degrees, range 0-90"]
    puts $fh "frame axis_angle_deg"
    set prev1 ""
    set prev2 ""
    set step 0
    for {set f $first} {$f <= $last} {incr f $stride} {
        incr step
        _id_progress $step $nanalyzed $f
        _id_frame_step $sel1 $sel2 $f
        lassign [_id_axes_pair $sel1 $sel2 prev1 prev2] axis1 axis2
        puts $fh "$f [_id_fmt [_id_axis_angle $axis1 $axis2]]"
    }
    close $fh
    $sel1 delete
    $sel2 delete
    puts "InterDomain: wrote $opt(out)"
    puts "InterDomain: columns: frame axis_angle_deg"
}

proc _id_cmd_dihedral {argv} {
    array set opt [_id_common $argv {out InterDomain_dihedral.dat} {}]
    lassign [_id_frames $opt(molid) $opt(first) $opt(last) $opt(stride)] first last stride nframes nanalyzed
    _id_open_pair $opt(molid) $opt(sel1) $opt(sel2) 3
    set fh [open $opt(out) w]
    _id_header $fh "InterDomain principal-axis dihedral" [list \
        molid $opt(molid) sel1 $opt(sel1) sel2 $opt(sel2) \
        frames $nframes frames_analyzed $nanalyzed \
        first $first last $last stride $stride \
        dihedral "signed twist about the COM-COM vector, degrees, range -180 to 180" \
        note "nan if an axis is parallel to the COM-COM vector"]
    puts $fh "frame dihedral_deg"
    set prev1 ""
    set prev2 ""
    set step 0
    for {set f $first} {$f <= $last} {incr f $stride} {
        incr step
        _id_progress $step $nanalyzed $f
        _id_frame_step $sel1 $sel2 $f
        lassign [_id_com_pair $sel1 $sel2] com1 com2 dist
        lassign [_id_axes_pair $sel1 $sel2 prev1 prev2] axis1 axis2
        puts $fh "$f [_id_fmt [_id_dihedral $axis1 $axis2 $com1 $com2]]"
    }
    close $fh
    $sel1 delete
    $sel2 delete
    puts "InterDomain: wrote $opt(out)"
    puts "InterDomain: columns: frame dihedral_deg"
}

proc _id_write_geometry_file {path molid sel1text sel2text first last stride} {
    # Write the same InterDomain_geometry.dat format used by "geometry",
    # and return compact rows: frame distance_A axis_angle_deg dihedral_deg.
    lassign [_id_frames $molid $first $last $stride] first last stride nframes nanalyzed
    _id_open_pair $molid $sel1text $sel2text 3
    set fh [open $path w]
    _id_header $fh "InterDomain geometry" [list \
        molid $molid sel1 $sel1text sel2 $sel2text \
        frames $nframes frames_analyzed $nanalyzed \
        first $first last $last stride $stride units Angstrom \
        axis_angle_deg "smaller angle between first principal axes, 0-90" \
        dihedral_deg "signed twist about COM-COM, -180 to 180"]
    puts $fh "frame com1_x com1_y com1_z com2_x com2_y com2_z distance_A axis_angle_deg dihedral_deg"
    set rows {}
    set prev1 ""
    set prev2 ""
    set step 0
    for {set f $first} {$f <= $last} {incr f $stride} {
        incr step
        _id_progress $step $nanalyzed $f
        _id_frame_step $sel1 $sel2 $f
        lassign [_id_com_pair $sel1 $sel2] com1 com2 dist
        lassign [_id_axes_pair $sel1 $sel2 prev1 prev2] axis1 axis2
        set angle [_id_axis_angle $axis1 $axis2]
        set dihedral [_id_dihedral $axis1 $axis2 $com1 $com2]
        puts $fh [format "%d %s %s %s %s %s %s %s %s %s" \
            $f \
            [_id_fmt [lindex $com1 0]] [_id_fmt [lindex $com1 1]] [_id_fmt [lindex $com1 2]] \
            [_id_fmt [lindex $com2 0]] [_id_fmt [lindex $com2 1]] [_id_fmt [lindex $com2 2]] \
            [_id_fmt $dist] [_id_fmt $angle] [_id_fmt $dihedral]]
        if {$dihedral eq "nan"} {
            puts "InterDomain: frame $f has no dihedral; skipped for clustering"
            continue
        }
        lappend rows [list $f $dist $angle $dihedral]
    }
    close $fh
    $sel1 delete
    $sel2 delete
    puts "InterDomain: wrote $path"
    if {[llength $rows] < 1} {
        error "InterDomain: no usable geometry rows written to $path"
    }
    return $rows
}

proc _id_cmd_geometry {argv} {
    array set opt [_id_common $argv {out InterDomain_geometry.dat} {}]
    _id_write_geometry_file $opt(out) $opt(molid) $opt(sel1) $opt(sel2) $opt(first) $opt(last) $opt(stride)
    puts "InterDomain: columns: frame com1_x com1_y com1_z com2_x com2_y com2_z distance_A axis_angle_deg dihedral_deg"
}

proc _id_frames_file {out} {
    if {[string match *.dat $out]} {
        return "[string range $out 0 end-4]_frames.dat"
    }
    return "${out}_frames.dat"
}

proc _id_resmap {sel} {
    set residue [dict create]
    set coord [dict create]
    foreach idx [$sel get index] xyz [$sel get {x y z}] \
            ri [$sel get resid] rn [$sel get resname] \
            ch [$sel get chain] sg [$sel get segname] {
        dict set residue $idx [list $ri $rn $ch $sg]
        dict set coord $idx $xyz
    }
    return [list $residue $coord]
}

proc _id_cmd_interface {argv} {
    array set opt [_id_common $argv {out InterDomain_interface.dat cutoff 5.0} {cutoff}]
    if {![string is double -strict $opt(cutoff)] || $opt(cutoff) <= 0} {
        error "InterDomain: -cutoff must be a positive number (Angstrom)"
    }
    lassign [_id_frames $opt(molid) $opt(first) $opt(last) $opt(stride)] first last stride nframes nanalyzed
    _id_open_pair $opt(molid) $opt(sel1) $opt(sel2) 1
    set pairFile $opt(out)
    set frameFile [_id_frames_file $opt(out)]
    set pf [open $pairFile w]
    set ff [open $frameFile w]
    set meta [list \
        molid $opt(molid) sel1 $opt(sel1) sel2 $opt(sel2) \
        cutoff_A $opt(cutoff) frames $nframes frames_analyzed $nanalyzed \
        first $first last $last stride $stride \
        occupancy "n_frames / frames_analyzed" \
        contact "at least one atom pair within the cutoff"]
    _id_header $pf "InterDomain interface residue pairs" $meta
    _id_header $ff "InterDomain interface per frame" $meta
    puts $pf "resid1 resname1 chain1 segname1 resid2 resname2 chain2 segname2 n_frames occupancy min_distance_A"
    puts $ff "frame n_atom_contacts n_residue_pairs n_res_sel1 n_res_sel2"

    set tally [dict create]
    set step 0
    for {set f $first} {$f <= $last} {incr f $stride} {
        incr step
        _id_progress $step $nanalyzed $f
        _id_frame_step $sel1 $sel2 $f
        lassign [_id_resmap $sel1] map1 coord1
        lassign [_id_resmap $sel2] map2 coord2
        set hit [measure contacts $opt(cutoff) $sel1 $sel2]
        set atoms1 [lindex $hit 0]
        set atoms2 [lindex $hit 1]
        set nAtom [llength $atoms1]
        set seen [dict create]
        set res1on [dict create]
        set res2on [dict create]
        foreach i $atoms1 j $atoms2 {
            if {![dict exists $map1 $i] || ![dict exists $map2 $j]} {
                continue
            }
            set r1 [dict get $map1 $i]
            set r2 [dict get $map2 $j]
            set key [join [concat $r1 $r2] "|"]
            set dist [veclength [vecsub [dict get $coord1 $i] [dict get $coord2 $j]]]
            if {![dict exists $seen $key]} {
                dict set seen $key $dist
            } elseif {$dist < [dict get $seen $key]} {
                dict set seen $key $dist
            }
            dict set res1on [join $r1 "|"] 1
            dict set res2on [join $r2 "|"] 1
        }
        dict for {key dist} $seen {
            if {![dict exists $tally $key]} {
                dict set tally $key [list 1 $dist]
            } else {
                lassign [dict get $tally $key] count mind
                incr count
                if {$dist < $mind} { set mind $dist }
                dict set tally $key [list $count $mind]
            }
        }
        puts $ff "$f $nAtom [dict size $seen] [dict size $res1on] [dict size $res2on]"
    }

    set rows {}
    dict for {key stats} $tally {
        lassign $stats count mind
        lassign [split $key "|"] ri1 rn1 ch1 sg1 ri2 rn2 ch2 sg2
        set occ [expr {double($count) / double($nanalyzed)}]
        set sortkey [format "%09d %s" [expr {100000000 - $count}] $key]
        lappend rows [list $sortkey $ri1 $rn1 $ch1 $sg1 $ri2 $rn2 $ch2 $sg2 $count $occ $mind]
    }
    set rows [lsort -index 0 $rows]
    foreach row $rows {
        puts $pf [format "%s %s %s %s %s %s %s %s %d %s %s" \
            [_id_field [lindex $row 1]] [_id_field [lindex $row 2]] \
            [_id_field [lindex $row 3]] [_id_field [lindex $row 4]] \
            [_id_field [lindex $row 5]] [_id_field [lindex $row 6]] \
            [_id_field [lindex $row 7]] [_id_field [lindex $row 8]] \
            [lindex $row 9] [_id_fmt [lindex $row 10]] [_id_fmt [lindex $row 11]]]
    }
    close $pf
    close $ff
    $sel1 delete
    $sel2 delete
    puts "InterDomain: wrote $pairFile"
    puts "InterDomain: wrote $frameFile"
    puts "InterDomain: pair columns: resid1 resname1 chain1 segname1 resid2 resname2 chain2 segname2 n_frames occupancy min_distance_A"
    puts "InterDomain: frame columns: frame n_atom_contacts n_residue_pairs n_res_sel1 n_res_sel2"
}

proc _id_draw_axes {molid sel length width colors} {
    set com [measure center $sel weight mass]
    set axes [lindex [measure inertia $sel] 1]
    foreach axis $axes color $colors {
        set norm [veclength $axis]
        if {$norm < 1e-8} { continue }
        set end [vecadd $com [vecscale [expr {$length / $norm}] $axis]]
        graphics $molid color $color
        graphics $molid line $com $end width $width
    }
}

proc _id_cmd_draw {argv} {
    array set opt [_id_args $argv {molid sel1 sel2 length} {
        molid top
        sel1 {}
        sel2 {}
        length 20
    }]
    if {[info exists opt(_help)]} {
        return
    }
    _id_need opt sel1
    if {![string is double -strict $opt(length)] || $opt(length) <= 0} {
        error "InterDomain: -length must be a positive number"
    }
    set sel1 [atomselect $opt(molid) $opt(sel1)]
    if {[$sel1 num] < 3} {
        set n [$sel1 num]
        $sel1 delete
        error "InterDomain: sel1 has $n atoms; drawing an axis needs at least 3"
    }
    _id_draw_axes $opt(molid) $sel1 $opt(length) 3 {red green blue}
    $sel1 delete
    puts "InterDomain: drew sel1 axes (red green blue) on molid $opt(molid)"
    if {$opt(sel2) ne ""} {
        set sel2 [atomselect $opt(molid) $opt(sel2)]
        if {[$sel2 num] < 3} {
            set n [$sel2 num]
            $sel2 delete
            error "InterDomain: sel2 has $n atoms; drawing an axis needs at least 3"
        }
        _id_draw_axes $opt(molid) $sel2 $opt(length) 2 {orange cyan magenta}
        $sel2 delete
        puts "InterDomain: drew sel2 axes (orange cyan magenta) on molid $opt(molid)"
    }
}

proc _id_scan_geometry {molid sel1 sel2 first last stride nanalyzed} {
    # Each row: frame distance_A axis_angle_deg dihedral_deg
    set rows {}
    set prev1 ""
    set prev2 ""
    set step 0
    for {set f $first} {$f <= $last} {incr f $stride} {
        incr step
        _id_progress $step $nanalyzed $f
        _id_frame_step $sel1 $sel2 $f
        lassign [_id_com_pair $sel1 $sel2] com1 com2 dist
        lassign [_id_axes_pair $sel1 $sel2 prev1 prev2] axis1 axis2
        set angle [_id_axis_angle $axis1 $axis2]
        set dih [_id_dihedral $axis1 $axis2 $com1 $com2]
        if {$dih eq "nan"} {
            puts "InterDomain: frame $f has no dihedral; skipped"
            continue
        }
        lappend rows [list $f $dist $angle $dih]
    }
    return $rows
}

proc _id_read_geometry {path} {
    # Each row: frame distance_A axis_angle_deg dihedral_deg
    set fh [open $path r]
    set header {}
    set i_frame -1
    set i_dist -1
    set i_ang -1
    set i_dih -1
    set rows {}
    while {[gets $fh line] >= 0} {
        set line [string trim $line]
        if {$line eq "" || [string match "#*" $line]} {
            continue
        }
        if {$header eq ""} {
            set header $line
            set i_frame [lsearch -exact $header frame]
            set i_dist [lsearch -exact $header distance_A]
            set i_ang [lsearch -exact $header axis_angle_deg]
            set i_dih [lsearch -exact $header dihedral_deg]
            if {$i_frame < 0 || $i_dist < 0 || $i_ang < 0 || $i_dih < 0} {
                close $fh
                error "InterDomain: $path needs columns frame, distance_A, axis_angle_deg, dihedral_deg"
            }
            continue
        }
        set dih [lindex $line $i_dih]
        set ang [lindex $line $i_ang]
        if {$dih eq "nan" || $ang eq "nan"} {
            continue
        }
        lappend rows [list [lindex $line $i_frame] [lindex $line $i_dist] $ang $dih]
    }
    close $fh
    if {[llength $rows] < 1} {
        error "InterDomain: no usable rows in $path"
    }
    return $rows
}

proc _id_parse_features {text} {
    # Accept "com,angle", "angle dihedral", or "com, angle, dihedral"
    set cleaned [string map {, " "} $text]
    set names {}
    foreach tok $cleaned {
        set tok [string tolower [string trim $tok]]
        if {$tok eq ""} { continue }
        switch -exact -- $tok {
            com - distance - distance_a { set name com }
            angle - axis_angle - axis_angle_deg { set name angle }
            dihedral - dihedral_deg { set name dihedral }
            default {
                error "InterDomain: unknown feature '$tok'. Use com, angle, dihedral"
            }
        }
        if {[lsearch -exact $names $name] < 0} {
            lappend names $name
        }
    }
    if {[llength $names] < 1} {
        error "InterDomain: -features needs at least one of: com, angle, dihedral"
    }
    return $names
}

proc _id_sqdist {a b} {
    set s 0.0
    foreach x $a y $b {
        set d [expr {$x - $y}]
        set s [expr {$s + $d * $d}]
    }
    return $s
}

proc _id_vec_zero {ndim} {
    set v {}
    for {set i 0} {$i < $ndim} {incr i} {
        lappend v 0.0
    }
    return $v
}

proc _id_vec_add {a b} {
    set out {}
    foreach x $a y $b {
        lappend out [expr {$x + $y}]
    }
    return $out
}

proc _id_vec_scale {v s} {
    set out {}
    foreach x $v {
        lappend out [expr {$x * $s}]
    }
    return $out
}

proc _id_cluster_features {rows feature_names} {
    # Build a z-scored feature vector from selected geometry columns.
    # com -> distance_A
    # angle -> axis_angle_deg (0-90)
    # dihedral -> sin/cos of dihedral so -180 and 180 stay neighbors
    set n [llength $rows]
    set pi [expr {acos(-1.0)}]
    set packed {}
    foreach row $rows {
        set dist [lindex $row 1]
        set ang [lindex $row 2]
        set dih [lindex $row 3]
        set vec {}
        foreach name $feature_names {
            switch -exact -- $name {
                com { lappend vec $dist }
                angle { lappend vec $ang }
                dihedral {
                    set th [expr {$dih * $pi / 180.0}]
                    lappend vec [expr {sin($th)}] [expr {cos($th)}]
                }
            }
        }
        lappend packed $vec
    }
    set ndim [llength [lindex $packed 0]]
    set mean [_id_vec_zero $ndim]
    foreach item $packed {
        set mean [_id_vec_add $mean $item]
    }
    set mean [_id_vec_scale $mean [expr {1.0 / double($n)}]]
    set var [_id_vec_zero $ndim]
    foreach item $packed {
        for {set d 0} {$d < $ndim} {incr d} {
            set diff [expr {[lindex $item $d] - [lindex $mean $d]}]
            lset var $d [expr {[lindex $var $d] + $diff * $diff}]
        }
    }
    set scale {}
    for {set d 0} {$d < $ndim} {incr d} {
        set s [expr {sqrt([lindex $var $d] / double($n))}]
        if {$s < 1e-8} { set s 1.0 }
        lappend scale $s
    }
    set features {}
    foreach item $packed {
        set vec {}
        for {set d 0} {$d < $ndim} {incr d} {
            lappend vec [expr {([lindex $item $d] - [lindex $mean $d]) / [lindex $scale $d]}]
        }
        lappend features $vec
    }
    return $features
}

proc _id_kmeans {features k maxiter} {
    set n [llength $features]
    set ndim [llength [lindex $features 0]]
    set centers [list [lindex $features 0]]
    for {set c 1} {$c < $k} {incr c} {
        set best_i 0
        set best_d -1.0
        for {set i 0} {$i < $n} {incr i} {
            set mind 1e300
            foreach center $centers {
                set d [_id_sqdist [lindex $features $i] $center]
                if {$d < $mind} { set mind $d }
            }
            if {$mind > $best_d} {
                set best_d $mind
                set best_i $i
            }
        }
        lappend centers [lindex $features $best_i]
    }
    set labels {}
    for {set i 0} {$i < $n} {incr i} { lappend labels 0 }
    for {set iter 0} {$iter < $maxiter} {incr iter} {
        set newlabels {}
        set sums {}
        set counts {}
        for {set c 0} {$c < $k} {incr c} {
            lappend sums [_id_vec_zero $ndim]
            lappend counts 0
        }
        set changed 0
        for {set i 0} {$i < $n} {incr i} {
            set point [lindex $features $i]
            set best 0
            set bestd 1e300
            for {set c 0} {$c < $k} {incr c} {
                set d [_id_sqdist $point [lindex $centers $c]]
                if {$d < $bestd} {
                    set bestd $d
                    set best $c
                }
            }
            if {[lindex $labels $i] != $best} { set changed 1 }
            lappend newlabels $best
            lset counts $best [expr {[lindex $counts $best] + 1}]
            lset sums $best [_id_vec_add [lindex $sums $best] $point]
        }
        for {set c 0} {$c < $k} {incr c} {
            set cnt [lindex $counts $c]
            if {$cnt == 0} {
                set far 0
                set fard -1.0
                for {set i 0} {$i < $n} {incr i} {
                    set d [_id_sqdist [lindex $features $i] [lindex $centers 0]]
                    if {$d > $fard} { set fard $d; set far $i }
                }
                lset centers $c [lindex $features $far]
                set changed 1
                continue
            }
            lset centers $c [_id_vec_scale [lindex $sums $c] [expr {1.0 / double($cnt)}]]
        }
        set labels $newlabels
        if {$iter > 0 && !$changed} {
            break
        }
    }
    return $labels
}

proc _id_circular_mean_deg {values} {
    set pi [expr {acos(-1.0)}]
    set sum_s 0.0
    set sum_c 0.0
    foreach v $values {
        set th [expr {$v * $pi / 180.0}]
        set sum_s [expr {$sum_s + sin($th)}]
        set sum_c [expr {$sum_c + cos($th)}]
    }
    return [expr {atan2($sum_s, $sum_c) * 180.0 / $pi}]
}

# geometcluster: k-means on z-scored geometry features (com / angle / dihedral).
# Groups frames by domain–domain geometry states, not by atomic RMSD.
proc _id_cmd_geometcluster {argv} {
    array set opt [_id_args $argv {molid sel1 sel2 first last stride ncluster outdir data maxiter features} {
        molid top
        sel1 {}
        sel2 {}
        first 0
        last end
        stride 1
        ncluster 3
        outdir geomet_clusters
        data {}
        maxiter 50
        features "com,dihedral"
    }]
    if {[info exists opt(_help)]} {
        return
    }
    if {![string is integer -strict $opt(ncluster)] || $opt(ncluster) < 2} {
        error "InterDomain: -ncluster must be an integer >= 2"
    }
    if {![string is integer -strict $opt(maxiter)] || $opt(maxiter) < 1} {
        error "InterDomain: -maxiter must be a positive integer"
    }
    set feature_names [_id_parse_features $opt(features)]
    set sel1 ""
    set sel2 ""
    if {$opt(data) eq ""} {
        _id_need opt sel1 sel2
        set opt(data) "InterDomain_geometry.dat"
        puts "InterDomain: no -data given; calculating $opt(data)"
        set rows [_id_write_geometry_file $opt(data) $opt(molid) $opt(sel1) $opt(sel2) $opt(first) $opt(last) $opt(stride)]
    } else {
        set rows [_id_read_geometry $opt(data)]
    }
    puts "InterDomain: geometcluster using data file $opt(data)"
    set n [llength $rows]
    if {$n < $opt(ncluster)} {
        error "InterDomain: $n frames are fewer than -ncluster $opt(ncluster)"
    }
    puts "InterDomain: geometcluster features = [join $feature_names {, }]"
    set features [_id_cluster_features $rows $feature_names]
    set raw [_id_kmeans $features $opt(ncluster) $opt(maxiter)]

    array set count {}
    foreach lab $raw { incr count($lab) }
    set ranked {}
    foreach lab [array names count] {
        lappend ranked [list $count($lab) $lab]
    }
    set ranked [lsort -integer -decreasing -index 0 $ranked]
    array set rank {}
    set cid 1
    foreach item $ranked {
        set rank([lindex $item 1]) $cid
        incr cid
    }

    array set members {}
    for {set i 0} {$i < $n} {incr i} {
        set id $rank([lindex $raw $i])
        lappend members($id) $i
    }

    file mkdir $opt(outdir)
    foreach old [glob -nocomplain -directory $opt(outdir) -types d cluster*] {
        if {[regexp {cluster[0-9]+$} $old]} {
            file delete -force $old
        }
    }

    set summary [open [file join $opt(outdir) summary.dat] w]
    _id_header $summary "InterDomain geometry clusters" [list \
        method "k-means on z-scored selected geometry features" \
        features [join $feature_names {,}] \
        ncluster $opt(ncluster) nframes $n outdir $opt(outdir)]
    puts $summary "cluster n_frames mean_distance_A mean_angle_deg mean_dihedral_deg representative_frame"
    set membership [open [file join $opt(outdir) frames.dat] w]
    _id_header $membership "InterDomain geometry cluster membership" [list \
        ncluster $opt(ncluster) features [join $feature_names {,}]]
    puts $membership "frame cluster distance_A axis_angle_deg dihedral_deg"

    set writer [atomselect $opt(molid) all]
    set ndim [llength [lindex $features 0]]
    foreach id [lsort -integer [array names members]] {
        set idxs $members($id)
        set folder [file join $opt(outdir) cluster$id]
        file mkdir $folder
        set dists {}
        set angs {}
        set dihs {}
        set centroid [_id_vec_zero $ndim]
        foreach i $idxs {
            set centroid [_id_vec_add $centroid [lindex $features $i]]
            set row [lindex $rows $i]
            lappend dists [lindex $row 1]
            lappend angs [lindex $row 2]
            lappend dihs [lindex $row 3]
        }
        set nmem [llength $idxs]
        set centroid [_id_vec_scale $centroid [expr {1.0 / double($nmem)}]]
        set best_i [lindex $idxs 0]
        set best_d 1e300
        foreach i $idxs {
            set d [_id_sqdist [lindex $features $i] $centroid]
            if {$d < $best_d} {
                set best_d $d
                set best_i $i
            }
        }
        set rep [lindex [lindex $rows $best_i] 0]
        set mean_d 0.0
        set mean_a 0.0
        foreach d $dists { set mean_d [expr {$mean_d + $d}] }
        foreach a $angs { set mean_a [expr {$mean_a + $a}] }
        set mean_d [expr {$mean_d / double($nmem)}]
        set mean_a [expr {$mean_a / double($nmem)}]
        set mean_th [_id_circular_mean_deg $dihs]
        puts $summary "$id $nmem [_id_fmt $mean_d] [_id_fmt $mean_a] [_id_fmt $mean_th] $rep"

        foreach i [lsort -integer $idxs] {
            set row [lindex $rows $i]
            set frame [lindex $row 0]
            puts $membership "$frame $id [_id_fmt [lindex $row 1]] [_id_fmt [lindex $row 2]] [_id_fmt [lindex $row 3]]"
            $writer frame $frame
            $writer writepdb [file join $folder frame_${frame}.pdb]
        }
        puts "InterDomain: cluster$id has $nmem frames, representative frame $rep"
    }
    $writer delete
    close $summary
    close $membership
    puts "InterDomain: wrote PDBs under $opt(outdir)/cluster1 .. cluster$opt(ncluster)"
    puts "InterDomain: wrote [file join $opt(outdir) summary.dat] and [file join $opt(outdir) frames.dat]"
}

proc _id_transform_point {matrix xyz} {
    set x [lindex $xyz 0]
    set y [lindex $xyz 1]
    set z [lindex $xyz 2]
    set r0 [lindex $matrix 0]
    set r1 [lindex $matrix 1]
    set r2 [lindex $matrix 2]
    return [list \
        [expr {[lindex $r0 0] * $x + [lindex $r0 1] * $y + [lindex $r0 2] * $z + [lindex $r0 3]}] \
        [expr {[lindex $r1 0] * $x + [lindex $r1 1] * $y + [lindex $r1 2] * $z + [lindex $r1 3]}] \
        [expr {[lindex $r2 0] * $x + [lindex $r2 1] * $y + [lindex $r2 2] * $z + [lindex $r2 3]}]]
}

proc _id_rmsd_coords {a b} {
    set n [llength $a]
    if {$n < 1 || $n != [llength $b]} {
        error "InterDomain: RMSD selections have different atom counts"
    }
    set s 0.0
    foreach p $a q $b {
        set dx [expr {[lindex $p 0] - [lindex $q 0]}]
        set dy [expr {[lindex $p 1] - [lindex $q 1]}]
        set dz [expr {[lindex $p 2] - [lindex $q 2]}]
        set s [expr {$s + $dx * $dx + $dy * $dy + $dz * $dz}]
    }
    return [expr {sqrt($s / double($n))}]
}

proc _id_pair_rmsd {matrix i j} {
    if {$i == $j} { return 0.0 }
    if {$i > $j} {
        set tmp $i
        set i $j
        set j $tmp
    }
    return [lindex [lindex $matrix $i] [expr {$j - $i - 1}]]
}

proc _id_avg_link {matrix group_a group_b} {
    set sum 0.0
    set n 0
    foreach i $group_a {
        foreach j $group_b {
            set sum [expr {$sum + [_id_pair_rmsd $matrix $i $j]}]
            incr n
        }
    }
    return [expr {$sum / double($n)}]
}

proc _id_hclust {matrix n ncluster cutoff} {
    for {set i 0} {$i < $n} {incr i} {
        set groups($i) [list $i]
        lappend ids $i
    }
    set next $n
    while {[llength $ids] > 1} {
        if {$ncluster ne "" && [llength $ids] <= $ncluster} {
            break
        }
        set best_d 1e300
        set best_a -1
        set best_b -1
        set nids [llength $ids]
        for {set a 0} {$a < $nids} {incr a} {
            for {set b [expr {$a + 1}]} {$b < $nids} {incr b} {
                set ia [lindex $ids $a]
                set ib [lindex $ids $b]
                set d [_id_avg_link $matrix $groups($ia) $groups($ib)]
                if {$d < $best_d} {
                    set best_d $d
                    set best_a $ia
                    set best_b $ib
                }
            }
        }
        if {$cutoff ne "" && $best_d > $cutoff} {
            break
        }
        set groups($next) [concat $groups($best_a) $groups($best_b)]
        unset groups($best_a)
        unset groups($best_b)
        set remaining {}
        foreach id $ids {
            if {$id != $best_a && $id != $best_b} {
                lappend remaining $id
            }
        }
        lappend remaining $next
        set ids $remaining
        incr next
    }
    set result {}
    foreach id $ids {
        lappend result $groups($id)
    }
    return $result
}

proc _id_write_rmsd_matrix {path frames matrix} {
    set fh [open $path w]
    set n [llength $frames]
    puts $fh "# RMSD of sel2 after fitting sel1"
    puts $fh "# rows and columns follow frame_order"
    puts $fh "# frame_order: $frames"
    puts $fh "# pandas: read_csv(path, sep=r'\\s+', comment='#', index_col=0)"
    puts -nonewline $fh "frame"
    foreach frame $frames {
        puts -nonewline $fh " $frame"
    }
    puts $fh ""
    for {set i 0} {$i < $n} {incr i} {
        puts -nonewline $fh [lindex $frames $i]
        for {set j 0} {$j < $n} {incr j} {
            puts -nonewline $fh " [_id_fmt [_id_pair_rmsd $matrix $i $j]]"
        }
        puts $fh ""
    }
    close $fh
}

proc _id_clear_cluster_dirs {outdir} {
    file mkdir $outdir
    foreach old [glob -nocomplain -directory $outdir -types d cluster*] {
        if {[regexp {cluster[0-9]+$} $old]} {
            file delete -force $old
        }
    }
}

# rmsdcluster: fit sel1, then hierarchical average-linkage on sel2 RMSD.
# Groups frames by structural pose of domain 2 relative to domain 1.
proc _id_cmd_rmsdcluster {argv} {
    array set opt [_id_args $argv {molid sel1 sel2 first last stride cutoff ncluster outdir} {
        molid top
        sel1 {}
        sel2 {}
        first 0
        last end
        stride 1
        cutoff {}
        ncluster {}
        outdir rmsd_clusters
    }]
    if {[info exists opt(_help)]} { return }
    _id_need opt sel1 sel2
    if {$opt(ncluster) eq "" && $opt(cutoff) eq ""} {
        error "InterDomain: set -ncluster or -cutoff"
    }
    if {$opt(ncluster) ne "" && $opt(cutoff) ne ""} {
        error "InterDomain: use -ncluster or -cutoff, not both"
    }
    if {$opt(ncluster) ne "" && (![string is integer -strict $opt(ncluster)] || $opt(ncluster) < 2)} {
        error "InterDomain: -ncluster must be an integer >= 2"
    }
    if {$opt(cutoff) ne "" && (![string is double -strict $opt(cutoff)] || $opt(cutoff) <= 0)} {
        error "InterDomain: -cutoff must be a positive RMSD in Angstrom"
    }
    lassign [_id_frames $opt(molid) $opt(first) $opt(last) $opt(stride)] first last stride nframes nanalyzed
    _id_open_pair $opt(molid) $opt(sel1) $opt(sel2) 3
    set ref1 [atomselect $opt(molid) $opt(sel1) frame $first]

    set frames {}
    set poses {}
    set step 0
    puts "InterDomain: fitting sel1 and recording sel2 for $nanalyzed frames"
    for {set f $first} {$f <= $last} {incr f $stride} {
        incr step
        _id_progress $step $nanalyzed $f
        _id_frame_step $sel1 $sel2 $f
        if {[$sel1 num] != [$ref1 num]} {
            $ref1 delete
            $sel1 delete
            $sel2 delete
            error "InterDomain: sel1 atom count changed at frame $f"
        }
        set matrix [measure fit $sel1 $ref1]
        set aligned {}
        foreach xyz [$sel2 get {x y z}] {
            lappend aligned [_id_transform_point $matrix $xyz]
        }
        if {[llength $poses] > 0 && [llength $aligned] != [llength [lindex $poses 0]]} {
            $ref1 delete
            $sel1 delete
            $sel2 delete
            error "InterDomain: sel2 atom count changed at frame $f"
        }
        lappend frames $f
        lappend poses $aligned
    }
    $ref1 delete
    $sel1 delete
    $sel2 delete

    set n [llength $frames]
    if {$n < 2} {
        error "InterDomain: rmsdcluster needs at least 2 frames"
    }
    puts "InterDomain: pairwise RMSD of sel2 for $n frames"
    set matrix {}
    for {set i 0} {$i < $n} {incr i} {
        set row {}
        for {set j [expr {$i + 1}]} {$j < $n} {incr j} {
            lappend row [_id_rmsd_coords [lindex $poses $i] [lindex $poses $j]]
        }
        lappend matrix $row
        if {$i == 0 || [expr {($i + 1) % 25}] == 0 || $i == [expr {$n - 1}]} {
            puts "InterDomain: RMSD row [expr {$i + 1}]/$n"
        }
    }

    set matrix_path [file join $opt(outdir) rmsd_matrix.dat]
    file mkdir $opt(outdir)
    _id_write_rmsd_matrix $matrix_path $frames $matrix
    puts "InterDomain: wrote $matrix_path ($n x $n)"

    set raw_groups [_id_hclust $matrix $n $opt(ncluster) $opt(cutoff)]
    set ranked {}
    foreach members $raw_groups {
        lappend ranked [list [llength $members] $members]
    }
    set ranked [lsort -integer -decreasing -index 0 $ranked]
    set k [llength $ranked]
    puts "InterDomain: hierarchical clustering produced $k clusters"

    _id_clear_cluster_dirs $opt(outdir)
    set summary [open [file join $opt(outdir) summary.dat] w]
    _id_header $summary "InterDomain RMSD clusters" [list \
        method "fit sel1, average-linkage hierarchical clustering on the sel2 RMSD matrix" \
        sel1 $opt(sel1) sel2 $opt(sel2) \
        reference_frame $first ncluster $opt(ncluster) cutoff_A $opt(cutoff) \
        frames_analyzed $n matrix $matrix_path outdir $opt(outdir)]
    puts $summary "cluster n_frames center_frame mean_rmsd_to_center_A"
    set membership [open [file join $opt(outdir) frames.dat] w]
    _id_header $membership "InterDomain RMSD cluster membership" [list \
        ncluster $opt(ncluster) cutoff_A $opt(cutoff)]
    puts $membership "frame cluster rmsd_to_center_A"

    set writer [atomselect $opt(molid) all]
    set id 0
    foreach item $ranked {
        incr id
        set members [lindex $item 1]
        set folder [file join $opt(outdir) cluster$id]
        file mkdir $folder
        set center -1
        set center_mean 1e300
        foreach candidate $members {
            set sum 0.0
            foreach other $members {
                set sum [expr {$sum + [_id_pair_rmsd $matrix $candidate $other]}]
            }
            set mean [expr {$sum / double([llength $members])}]
            if {$mean < $center_mean} {
                set center_mean $mean
                set center $candidate
            }
        }
        set center_frame [lindex $frames $center]
        puts $summary "$id [llength $members] $center_frame [_id_fmt $center_mean]"
        foreach member [lsort -integer $members] {
            set frame [lindex $frames $member]
            puts $membership "$frame $id [_id_fmt [_id_pair_rmsd $matrix $member $center]]"
            $writer frame $frame
            $writer writepdb [file join $folder frame_${frame}.pdb]
        }
        puts "InterDomain: cluster$id has [llength $members] frames, center frame $center_frame"
    }
    $writer delete
    close $summary
    close $membership
    puts "InterDomain: wrote PDBs under $opt(outdir)/cluster1 .. cluster$k"
    puts "InterDomain: wrote [file join $opt(outdir) summary.dat] and [file join $opt(outdir) frames.dat]"
}

foreach _id_old {interdomain interDomainDistance interDomainAnalysis getInterface drawPrincipalAxes} {
    catch {rename $_id_old {}}
}
unset _id_old

puts "InterDomain loaded. Type: InterDomain help"
