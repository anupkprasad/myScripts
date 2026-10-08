# trajClustering.tcl
# Single-selection RMSD clustering for trajectories loaded in VMD.
#
# Loaded from callScripts.tcl, or:
#   source ~/myScripts/simulations/vmdTcl/mdAnalysis/trajClustering.tcl
#   TrajCluster help

proc TrajCluster {args} {
    if {[llength $args] == 0} {
        _tc_help
        return
    }
    set cmd [lindex $args 0]
    set rest [lrange $args 1 end]
    if {$cmd in {help --help -help -h}} {
        _tc_help
        return
    }
    switch -exact -- $cmd {
        rmsd    { _tc_cmd_rmsd $rest }
        matrix  { _tc_cmd_matrix $rest }
        cluster { _tc_cmd_cluster $rest }
        python  { _tc_cmd_python $rest }
        default {
            puts "TrajCluster: unknown command '$cmd'"
            _tc_help
        }
    }
}

proc _tc_help {} {
    puts {
TrajCluster — RMSD clustering of one atom selection on a VMD molecule

Usage:
  TrajCluster <command> [options]
  TrajCluster help

Commands:
  rmsd      RMSD of each frame to the first analyzed frame (after fitting)
  matrix    Write the N x N pairwise RMSD matrix
  cluster   Hierarchical clustering from that matrix; write PDBs
  python    Call the standalone Python module (needs MDAnalysis)
  help      Show this message

Common options:
  -molid <id|top>   molecule id (default: top)
  -sel <text>       VMD selection (default: protein and name CA)
  -first <frame>    first frame, inclusive (default: 0)
  -last <frame|end> last frame, inclusive (default: end)
  -stride <n>       use every n-th frame (default: 1)
  -out <file>       output for rmsd / matrix

cluster also accepts one of:
  -ncluster <k>     cut the tree into k groups
  -cutoff <A>       cut where average-linkage RMSD exceeds this value
  -outdir <dir>     output directory (default: traj_clusters)

How it works:
  For every pair of frames the selection is fitted (i onto j), then RMSD is
  computed. Fitting uses the loaded molecule; coordinates are transformed in
  Tcl so the trajectory itself is not moved. Clustering is average-linkage
  hierarchical clustering on that pairwise RMSD matrix. cluster1 is largest.

Examples:
  TrajCluster help
  TrajCluster rmsd -molid 0 -sel "protein and name CA" -out traj_rmsd.dat
  TrajCluster matrix -sel "backbone" -stride 5 -out matrix.dat
  TrajCluster cluster -sel "protein and name CA" -ncluster 4 -outdir traj_clusters
  TrajCluster cluster -sel "name CA" -cutoff 2.5 -first 0 -last 200

Python module (optional, outside VMD or via TrajCluster python):
  ~/myScripts/simulations/trajectory/trajClustering.py
  Heatmap of a written matrix:
    conda activate ComBio
    python ~/myScripts/simulations/trajectory/trajClustering.py plot \\
        -matrix traj_clusters/rmsd_matrix.dat -out heatmap.png
}
}

proc _tc_args {argv allowed defaults} {
    array set opt $defaults
    set n [llength $argv]
    for {set i 0} {$i < $n} {incr i} {
        set tok [lindex $argv $i]
        if {![string match -* $tok]} {
            error "TrajCluster: unexpected '$tok'. Type: TrajCluster help"
        }
        set name [string range $tok 1 end]
        if {$name in {h help -help}} {
            _tc_help
            set opt(_help) 1
            return [array get opt]
        }
        if {[lsearch -exact $allowed $name] < 0} {
            error "TrajCluster: unknown option -$name. Allowed: -[join $allowed { -}]"
        }
        incr i
        if {$i >= $n} {
            error "TrajCluster: -$name needs a value"
        }
        set opt($name) [lindex $argv $i]
    }
    return [array get opt]
}

proc _tc_need {optName args} {
    upvar $optName opt
    foreach key $args {
        if {![info exists opt($key)] || $opt($key) eq ""} {
            error "TrajCluster: missing -$key. Type: TrajCluster help"
        }
    }
}

proc _tc_frames {molid first last stride} {
    set nframes [molinfo $molid get numframes]
    if {$nframes < 1} {
        error "TrajCluster: molecule $molid has no frames"
    }
    if {$last eq "end" || $last eq ""} {
        set last [expr {$nframes - 1}]
    }
    if {![string is integer -strict $first] || ![string is integer -strict $last] || ![string is integer -strict $stride]} {
        error "TrajCluster: -first, -last, and -stride must be integers (-last may be end)"
    }
    if {$first < 0 || $last >= $nframes || $first > $last || $stride < 1} {
        error "TrajCluster: bad range first=$first last=$last stride=$stride (nframes=$nframes)"
    }
    set nanalyzed [expr {($last - $first) / $stride + 1}]
    return [list $first $last $stride $nframes $nanalyzed]
}

proc _tc_fmt {x} {
    return [format "%.4f" $x]
}

proc _tc_progress {step total frame} {
    if {$step == 1 || $step == $total || $step % 25 == 0} {
        puts "TrajCluster: frame $frame ($step/$total)"
    }
}

proc _tc_header {fh title meta} {
    puts $fh "# $title"
    foreach {key value} $meta {
        puts $fh "# $key: $value"
    }
    puts $fh "# pandas: read_csv(path, sep=r'\\s+', comment='#')"
}

proc _tc_transform_point {matrix xyz} {
    set x [lindex $xyz 0]
    set y [lindex $xyz 1]
    set z [lindex $xyz 2]
    set r0 [lindex $matrix 0]
    set r1 [lindex $matrix 1]
    set r2 [lindex $matrix 2]
    return [list \
        [expr {[lindex $r0 0]*$x + [lindex $r0 1]*$y + [lindex $r0 2]*$z + [lindex $r0 3]}] \
        [expr {[lindex $r1 0]*$x + [lindex $r1 1]*$y + [lindex $r1 2]*$z + [lindex $r1 3]}] \
        [expr {[lindex $r2 0]*$x + [lindex $r2 1]*$y + [lindex $r2 2]*$z + [lindex $r2 3]}]]
}

proc _tc_transform_coords {matrix coords} {
    set out {}
    foreach xyz $coords {
        lappend out [_tc_transform_point $matrix $xyz]
    }
    return $out
}

proc _tc_rmsd_coords {a b} {
    set n [llength $a]
    if {$n < 1 || $n != [llength $b]} {
        error "TrajCluster: RMSD atom counts differ"
    }
    set s 0.0
    foreach p $a q $b {
        set dx [expr {[lindex $p 0] - [lindex $q 0]}]
        set dy [expr {[lindex $p 1] - [lindex $q 1]}]
        set dz [expr {[lindex $p 2] - [lindex $q 2]}]
        set s [expr {$s + $dx*$dx + $dy*$dy + $dz*$dz}]
    }
    return [expr {sqrt($s / double($n))}]
}

proc _tc_collect_frames {molid seltext first last stride} {
    # Store raw selection coordinates for each analyzed frame.
    set sel [atomselect $molid $seltext]
    if {[$sel num] < 3} {
        set n [$sel num]
        $sel delete
        error "TrajCluster: selection has $n atoms; need at least 3"
    }
    set frames {}
    set poses {}
    set step 0
    set nanalyzed [expr {($last - $first) / $stride + 1}]
    puts "TrajCluster: reading selection coordinates for $nanalyzed frames"
    for {set f $first} {$f <= $last} {incr f $stride} {
        incr step
        _tc_progress $step $nanalyzed $f
        $sel frame $f
        $sel update
        lappend frames $f
        lappend poses [$sel get {x y z}]
    }
    set natom [$sel num]
    $sel delete
    return [list $frames $poses $natom]
}

proc _tc_pair_fit_rmsd {molid seltext frame_i frame_j} {
    # Fit frame_i onto frame_j for the selection, then RMSD.
    # Uses live atomselects on different frames; applies the transform in Tcl
    # so the trajectory coordinates are never moved.
    set si [atomselect $molid $seltext frame $frame_i]
    set sj [atomselect $molid $seltext frame $frame_j]
    if {[$si num] != [$sj num] || [$si num] < 3} {
        set ni [$si num]
        set nj [$sj num]
        $si delete
        $sj delete
        error "TrajCluster: bad selection sizes at frames $frame_i ($ni) and $frame_j ($nj)"
    }
    set mat [measure fit $si $sj]
    set aligned [_tc_transform_coords $mat [$si get {x y z}]]
    set value [_tc_rmsd_coords $aligned [$sj get {x y z}]]
    $si delete
    $sj delete
    return $value
}

proc _tc_build_matrix {molid seltext frames} {
    # True pairwise fit+RMSD for every unique pair of frames.
    set n [llength $frames]
    set matrix {}
    for {set i 0} {$i < $n} {incr i} {
        set row {}
        set fi [lindex $frames $i]
        for {set j [expr {$i + 1}]} {$j < $n} {incr j} {
            set fj [lindex $frames $j]
            lappend row [_tc_pair_fit_rmsd $molid $seltext $fi $fj]
        }
        lappend matrix $row
        if {$i == 0 || [expr {($i + 1) % 10}] == 0 || $i == [expr {$n - 1}]} {
            puts "TrajCluster: RMSD row [expr {$i + 1}]/$n"
        }
    }
    return $matrix
}

proc _tc_pair_rmsd {matrix i j} {
    if {$i == $j} { return 0.0 }
    if {$i > $j} {
        set tmp $i
        set i $j
        set j $tmp
    }
    return [lindex [lindex $matrix $i] [expr {$j - $i - 1}]]
}

proc _tc_avg_link {matrix group_a group_b} {
    set sum 0.0
    set n 0
    foreach i $group_a {
        foreach j $group_b {
            set sum [expr {$sum + [_tc_pair_rmsd $matrix $i $j]}]
            incr n
        }
    }
    return [expr {$sum / double($n)}]
}

proc _tc_hclust {matrix n ncluster cutoff} {
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
                set d [_tc_avg_link $matrix $groups($ia) $groups($ib)]
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

proc _tc_write_matrix {path frames matrix} {
    set fh [open $path w]
    set n [llength $frames]
    _tc_header $fh "TrajCluster pairwise RMSD matrix" [list \
        frame_order $frames nframes $n]
    puts -nonewline $fh "frame"
    foreach frame $frames { puts -nonewline $fh " $frame" }
    puts $fh ""
    for {set i 0} {$i < $n} {incr i} {
        puts -nonewline $fh [lindex $frames $i]
        for {set j 0} {$j < $n} {incr j} {
            puts -nonewline $fh " [_tc_fmt [_tc_pair_rmsd $matrix $i $j]]"
        }
        puts $fh ""
    }
    close $fh
}

proc _tc_clear_dirs {outdir} {
    file mkdir $outdir
    foreach old [glob -nocomplain -directory $outdir -types d cluster*] {
        if {[regexp {cluster[0-9]+$} $old]} {
            file delete -force $old
        }
    }
}

proc _tc_cmd_rmsd {argv} {
    array set opt [_tc_args $argv {molid sel first last stride out} {
        molid top
        sel {protein and name CA}
        first 0
        last end
        stride 1
        out traj_rmsd.dat
    }]
    if {[info exists opt(_help)]} { return }
    lassign [_tc_frames $opt(molid) $opt(first) $opt(last) $opt(stride)] first last stride nframes nanalyzed
    lassign [_tc_collect_frames $opt(molid) $opt(sel) $first $last $stride] frames poses natom
    set fh [open $opt(out) w]
    _tc_header $fh "TrajCluster RMSD to reference" [list \
        molid $opt(molid) sel $opt(sel) reference_frame $first \
        first $first last $last stride $stride \
        note "each frame is fitted to the first analyzed frame before RMSD"]
    puts $fh "frame rmsd_A"
    set ref_frame [lindex $frames 0]
    foreach frame $frames {
        puts $fh "$frame [_tc_fmt [_tc_pair_fit_rmsd $opt(molid) $opt(sel) $frame $ref_frame]]"
    }
    close $fh
    puts "TrajCluster: wrote $opt(out)"
}

proc _tc_cmd_matrix {argv} {
    array set opt [_tc_args $argv {molid sel first last stride out} {
        molid top
        sel {protein and name CA}
        first 0
        last end
        stride 1
        out traj_rmsd_matrix.dat
    }]
    if {[info exists opt(_help)]} { return }
    lassign [_tc_frames $opt(molid) $opt(first) $opt(last) $opt(stride)] first last stride nframes nanalyzed
    lassign [_tc_collect_frames $opt(molid) $opt(sel) $first $last $stride] frames poses natom
    puts "TrajCluster: pairwise fit+RMSD for [llength $frames] frames"
    set matrix [_tc_build_matrix $opt(molid) $opt(sel) $frames]
    _tc_write_matrix $opt(out) $frames $matrix
    puts "TrajCluster: wrote $opt(out)"
}

proc _tc_cmd_cluster {argv} {
    array set opt [_tc_args $argv {molid sel first last stride ncluster cutoff outdir} {
        molid top
        sel {protein and name CA}
        first 0
        last end
        stride 1
        ncluster {}
        cutoff {}
        outdir traj_clusters
    }]
    if {[info exists opt(_help)]} { return }
    if {$opt(ncluster) eq "" && $opt(cutoff) eq ""} {
        error "TrajCluster: set -ncluster or -cutoff"
    }
    if {$opt(ncluster) ne "" && $opt(cutoff) ne ""} {
        error "TrajCluster: use -ncluster or -cutoff, not both"
    }
    if {$opt(ncluster) ne "" && (![string is integer -strict $opt(ncluster)] || $opt(ncluster) < 2)} {
        error "TrajCluster: -ncluster must be an integer >= 2"
    }
    if {$opt(cutoff) ne "" && (![string is double -strict $opt(cutoff)] || $opt(cutoff) <= 0)} {
        error "TrajCluster: -cutoff must be a positive RMSD in Angstrom"
    }
    lassign [_tc_frames $opt(molid) $opt(first) $opt(last) $opt(stride)] first last stride nframes nanalyzed
    if {$nanalyzed < 2} {
        error "TrajCluster: need at least 2 frames"
    }
    lassign [_tc_collect_frames $opt(molid) $opt(sel) $first $last $stride] frames poses natom
    set n [llength $frames]
    puts "TrajCluster: pairwise fit+RMSD for $n frames"
    set matrix [_tc_build_matrix $opt(molid) $opt(sel) $frames]
    file mkdir $opt(outdir)
    set matrix_path [file join $opt(outdir) rmsd_matrix.dat]
    _tc_write_matrix $matrix_path $frames $matrix
    puts "TrajCluster: wrote $matrix_path"

    set raw [_tc_hclust $matrix $n $opt(ncluster) $opt(cutoff)]
    set ranked {}
    foreach members $raw {
        lappend ranked [list [llength $members] $members]
    }
    set ranked [lsort -integer -decreasing -index 0 $ranked]
    set k [llength $ranked]
    puts "TrajCluster: hierarchical clustering produced $k clusters"

    _tc_clear_dirs $opt(outdir)
    set summary [open [file join $opt(outdir) summary.dat] w]
    _tc_header $summary "TrajCluster RMSD clusters" [list \
        method "pairwise fit+RMSD of the selection, average-linkage hierarchical clustering" \
        molid $opt(molid) sel $opt(sel) \
        ncluster $opt(ncluster) cutoff_A $opt(cutoff) \
        frames_analyzed $n matrix $matrix_path outdir $opt(outdir)]
    puts $summary "cluster n_frames center_frame mean_rmsd_to_center_A"
    set membership [open [file join $opt(outdir) frames.dat] w]
    _tc_header $membership "TrajCluster membership" [list \
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
                set sum [expr {$sum + [_tc_pair_rmsd $matrix $candidate $other]}]
            }
            set mean [expr {$sum / double([llength $members])}]
            if {$mean < $center_mean} {
                set center_mean $mean
                set center $candidate
            }
        }
        set center_frame [lindex $frames $center]
        puts $summary "$id [llength $members] $center_frame [_tc_fmt $center_mean]"
        foreach member [lsort -integer $members] {
            set frame [lindex $frames $member]
            puts $membership "$frame $id [_tc_fmt [_tc_pair_rmsd $matrix $member $center]]"
            $writer frame $frame
            $writer writepdb [file join $folder frame_${frame}.pdb]
        }
        puts "TrajCluster: cluster$id has [llength $members] frames, center frame $center_frame"
    }
    $writer delete
    close $summary
    close $membership
    puts "TrajCluster: wrote PDBs under $opt(outdir)/cluster1 .. cluster$k"
    puts "TrajCluster: wrote [file join $opt(outdir) summary.dat] and [file join $opt(outdir) frames.dat]"
}

proc _tc_find_python {} {
    foreach candidate {
        /home/anup/miniconda3/envs/ComBio/bin/python
        /home/anup/miniconda3/bin/python
    } {
        if {[file executable $candidate]} {
            return $candidate
        }
    }
    if {![catch {exec which python3} path]} {
        return [string trim $path]
    }
    error "TrajCluster: no Python interpreter found for the python backend"
}

proc _tc_cmd_python {argv} {
    if {[llength $argv] == 0 || [lindex $argv 0] in {help -h --help}} {
        puts "TrajCluster python — forward arguments to trajClustering.py"
        puts "Example:"
        puts {  TrajCluster python cluster -top top.pdb -traj traj.xtc -sel "name CA" -ncluster 4}
        return
    }
    set py [_tc_find_python]
    set script [file normalize ~/myScripts/simulations/trajectory/trajClustering.py]
    if {![file exists $script]} {
        error "TrajCluster: missing $script"
    }
    puts "TrajCluster: running $py $script $argv"
    set status [catch {eval exec $py [list $script] {*}$argv 2>@1} result]
    puts $result
    if {$status} {
        error "TrajCluster: python backend failed"
    }
}

puts "TrajCluster loaded. Type: TrajCluster help"
