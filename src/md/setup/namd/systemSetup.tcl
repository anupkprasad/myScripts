#!/usr/bin/tclsh
# VMD helpers for PSF generation / solvation / ionization and NAMD system setup.
# Anup K. Prasad — anupkprasad121@gmail.com
#
# Loaded via src/md/vmd_autoload.tcl (whitelist). Uses globals:
#   pathSetup    → src/md/setup/namd
#   pathVmdScript → repo root (~/myScripts)

proc get_min {items} {
    set min [lindex $items 0]
    foreach i $items {
        set ai [expr {abs($i)}]
        if {$ai < $min} {
            set min $ai
        }
    }
    puts "min is = $min"
    return $min
}

proc get_max {items} {
    set max [lindex $items 0]
    foreach i $items {
        set ai [expr {abs($i)}]
        if {$ai > $max} {
            set max $ai
        }
    }
    puts "max is = $max"
    return $max
}

proc moveSysOrigin {{molid top} {newcom {0 0 0}}} {
    set pdbfile [molinfo $molid get name]
    set inputname [lindex [split $pdbfile .] 0]
    set a [atomselect $molid all]
    set c [measure center $a]
    $a moveby [vecscale $c -1.0]
    $a moveby $newcom
    $a writepdb ${inputname}Cent.pdb
    $a delete
}

proc protMoveOnLipSurface {{pdbfile prot_lip.pdb}} {
    mol new $pdbfile
    set prot [atomselect top "protein"]
    moveSysOrigin
    $prot moveby {0 0 33.34}  ;# Z of lip surface is 23.34 and +10 pad
    set sel [atomselect top all]
    $sel writepdb prot_lip.pdb
    $prot delete
    $sel delete
}

proc rotateMol {{axis x} {angle 5}} {
    set sel [atomselect top "protein"]
    set com [measure center $sel weight mass]
    set matrix [transaxis $axis $angle]
    $sel moveby [vecscale -1.0 $com]
    $sel move $matrix
    $sel moveby $com
    $sel delete
}

proc calculate_distance {point1 point2} {
    set x1 [lindex $point1 0]
    set y1 [lindex $point1 1]
    set z1 [lindex $point1 2]
    set x2 [lindex $point2 0]
    set y2 [lindex $point2 1]
    set z2 [lindex $point2 2]
    set dx [expr {$x2 - $x1}]
    set dy [expr {$y2 - $y1}]
    set dz [expr {$z2 - $z1}]
    return [expr {sqrt($dx * $dx + $dy * $dy + $dz * $dz)}]
}

proc getBoxSetup {{molid top} {pad 7} {cube_padsize yes}} {
    set a [atomselect $molid all]
    set c [measure center $a]
    set box {{? ? ?} {? ? ?}}
    set basisvec {? ? ?}
    set origin {? ? ?}
    set minmax [measure minmax $a]

    if {$cube_padsize eq "yes"} {
        set j [lindex $minmax 0]
        set k [lindex $minmax 1]
        set distance [calculate_distance $j $k]
        set cube_halfside [expr {0.5 * $distance + $pad}]
        puts "box size = $cube_halfside"
        foreach d {0 1 2} {
            set cd [lindex $c $d]
            lset box 0 $d [expr {$cd - $cube_halfside}]
            lset box 1 $d [expr {$cd + $cube_halfside}]
            set b0 [lindex $box 0 $d]
            set b1 [lindex $box 1 $d]
            lset basisvec $d [expr {$b1 - $b0}]
            lset origin $d [expr {0.5 * ($b1 + $b0)}]
        }
    } else {
        foreach d {0 1 2} {
            set lo [lindex $minmax 0 $d]
            set hi [lindex $minmax 1 $d]
            lset box 0 $d [expr {$lo - $pad}]
            lset box 1 $d [expr {$hi + $pad}]
            set b0 [lindex $box 0 $d]
            set b1 [lindex $box 1 $d]
            lset basisvec $d [expr {$b1 - $b0}]
            lset origin $d [expr {0.5 * ($b1 + $b0)}]
        }
    }
    puts "Box size: $box"
    $a delete
    return $box
}

proc getBoxProtLip {pdbfile {pad 7}} {
    set inputname $pdbfile
    mol new ${inputname}.pdb
    set a [atomselect top "segname XO1"]
    set box [measure minmax $a]
    set zmin [lindex $box 0 2]
    set boxlower [measure minmax $a]
    lset boxlower 0 2 [expr {$zmin - $pad}]
    lset boxlower 1 2 $zmin
    set zmax [lindex $box 1 2]
    set boxupper [measure minmax $a]
    lset boxupper 0 2 $zmax
    lset boxupper 1 2 [expr {$zmax + $pad}]
    puts "boxupper size: $boxupper"
    puts "boxlower size: $boxlower"
    set boxboth [list $boxupper $boxlower]
    $a delete
    return $boxboth
}

proc getpsf {pdbfile {amidated "No"}} {
    ## Generate pdb/psf; optional CTER amidation
    set inputname $pdbfile
    resetpsf
    mol new $inputname.pdb
    set sel [atomselect top "protein"]
    set cter_index [lindex [$sel get resid] end]
    $sel delete
    topology "/usr/local/lib/vmd/plugins/noarch/tcl/readcharmmtop1.2/top_all36_prot.rtf"
    segment P1 [list pdb $inputname.pdb]
    if {$amidated ne "No"} {
        patch CT2 P1:$cter_index
    }
    regenerate angles dihedrals
    coordpdb $inputname.pdb P1
    guesscoord
    writepdb ${inputname}p.pdb
    writepsf ${inputname}p.psf
    mol new ${inputname}p.psf
    mol addfile ${inputname}p.pdb
    puts "\npsf finished\n"
}

proc getsolv {input box} {
    puts "\n############  solvation start...    ##############\n"
    package require solvate
    solvate ${input}.psf ${input}.pdb -minmax $box -o solvate
}

proc getionze {input} {
    puts "\n######### autoionize start... ###########\n"
    package require autoionize
    autoionize -psf ${input}.psf -pdb ${input}.pdb -sc 0.15 -o ionized
}

proc constrinedFixedFiles {{molid top}} {
    set all [atomselect $molid all]
    set p [atomselect $molid "protein and backbone"]
    $all set beta 0
    $p set beta 1
    $all writepdb ./fixed.pdb

    set q [atomselect $molid "protein and name CA"]
    $all set beta 0
    $q set beta 0.5
    $all writepdb ./constrained.pdb
    $all delete
    $p delete
    $q delete
}

proc restraintRefFiles {{molid top}} {
    set lip_head_popc {N C12 C13 C14 C15 C11 P O13 O14 O12 O11}
    set lip_head_popg {C13 OC3 C12 OC2 C11 P O13 O14 O12 O11}

    set all [atomselect $molid all]
    set p [atomselect $molid "protein and backbone"]
    $all set beta 0
    $p set beta 1
    $all writepdb ../restraints/prot_posres.ref
    $p delete
    $all delete

    set all [atomselect $molid "all and not water"]
    set popc_UH [atomselect $molid \
        "segname XO1 and not resid > 100 and resname POPC and name $lip_head_popc"]
    $all set beta 0
    $popc_UH set beta 1
    $all writepdb ../restraints/popc_head_upper.ref
    $popc_UH delete

    set popc_LH [atomselect $molid \
        "segname XO1 and resid > 100 and resname POPC and name $lip_head_popc"]
    $all set beta 0
    $popc_LH set beta 1
    $all writepdb ../restraints/popc_head_lower.ref
    $popc_LH delete
    $all delete

    set all [atomselect $molid "all and not water"]
    set popg_UH [atomselect $molid \
        "segname XO1 and not resid > 100 and resname POPC and name $lip_head_popg"]
    $all set beta 0
    $popg_UH set beta 1
    $all writepdb ../restraints/popg_head_upper.ref
    $popg_UH delete

    set popg_LH [atomselect $molid \
        "segname XO1 and resid > 100 and resname POPC and name $lip_head_popg"]
    $all set beta 0
    $popg_LH set beta 1
    $all writepdb ../restraints/popg_head_lower.ref
    $popg_LH delete
    $all delete
}

proc PBC_info {{molid top}} {
    set all [atomselect $molid all]
    set minmax [measure minmax $all]
    set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]]
    set a "set CBV1 {[lindex $vec 0] 0 0}"
    set b "set CBV2 {0 [lindex $vec 1] 0}"
    set c "set CBV3 {0 0 [lindex $vec 2]}"
    set center "set ori {[measure center $all]}"
    $all delete
    return [list $a\n$b\n$c\n$center]
}

proc insert_lines_in_file {filename line_number new_lines} {
    set in [open $filename r]
    set lines [split [read $in] "\n"]
    close $in
    set lines [linsert $lines $line_number $new_lines]
    set out [open $filename w]
    puts $out [join $lines "\n"]
    close $out
}

########### MD setup functions ###########
proc sysSetup {{molid top} {cter_amidation "No"}} {
    global pathSetup
    set pdbfile [molinfo $molid get name]
    set inputname [lindex [split $pdbfile .] 0]
    set box [getBoxSetup $molid 7]
    getsolv ${inputname} $box
    getionze solvate
    constrinedFixedFiles

    exec cp -r [file join $pathSetup protein equilibration] ./
    set pbcinfo [lindex [getCellDimention] 0]
    insert_lines_in_file ./equilibration/equilibration.inp 50 $pbcinfo

    exec cp -r [file join $pathSetup protein prod1] ./
    exec cp -r [file join $pathSetup protein parameter] ./

    getSysInfo top yes
}

proc sysSetupSDS {{molid top} {cter_amidation "No"}} {
    global pathSetup
    set pdbfile [molinfo $molid get name]
    set inputname [lindex [split $pdbfile .] 0]
    set box [getBoxSetup $molid 7]
    getsolv ${inputname} $box
    getionze solvate
    constrinedFixedFiles

    exec cp -r [file join $pathSetup proteinSDS equilibration] ./
    set pbcinfo [lindex [getCellDimention] 0]
    insert_lines_in_file ./equilibration/equilibration.inp 50 $pbcinfo

    exec cp -r [file join $pathSetup proteinSDS prod1] ./
    exec cp -r [file join $pathSetup proteinSDS parameter] ./
    set sds [atomselect top "segname XO1"]
    set min [expr {[get_min [$sds get index]] + 1}]
    set max [expr {[get_max [$sds get index]] + 1}]
    set line "\t\t\tatomnumbersrange {$min-$max}"
    insert_lines_in_file ./parameter/com_colvar.in 11 $line
    $sds delete

    sys_info top yes
}

proc sysSetupLipbilayr {pdbfile {cter_amidation "No"}} {
    global pathSetup pathVmdScript
    set box [getBoxProtLip ${pdbfile}_cp 20]
    puts "box1: [lindex $box 0]"
    puts "box2: [lindex $box 1]"
    getsolv ${pdbfile}_cp [lindex $box 0]
    getsolv ${pdbfile}_cpw [lindex $box 1]
    getionze ${pdbfile}_cpww

    if {![file isdirectory ../restraints]} { file mkdir ../restraints }
    restraintRefFiles
    file copy -force \
        [file join $pathSetup refsParaInp lipProt lipColvarTemplate.dat] \
        ../restraints

    if {![file isdirectory ../equilbration]} { file mkdir ../equilbration }
    file copy -force [file join $pathSetup refsParaInp namd.sh] ../equilbration
    for {set a 1} {$a < 7} {incr a} {
        file copy -force \
            [file join $pathSetup refsParaInp lipProt step6.${a}_equilibration.inp] \
            ../equilbration
    }

    # CHARMM/NAMD parameters: data/forcefields/namd/toppar
    file copy -force \
        [file join $pathVmdScript data forcefields namd toppar] \
        ../
    sysinfo_inp
    sys_info top yes
}
