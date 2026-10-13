#!/usr/bin/tclsh
#### Get transformation matrix @2023
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy

proc transmatrix {filename} {
    set fp [open $filename r]
    set matrix {}
    set count 0
    while {[gets $fp line]!=-1} {
        incr count
        if {$count==27} {
            set a $line
            lappend matrix $line
            }
        if {$count==28} {
            set b $line
            lappend matrix $line
            }
        if {$count==29} {
            set c $line
            lappend matrix $line
            }
        if {$count==30} {
            set d $line
            lappend matrix $line
                break
            }
    }
    close $fp
    puts "Transformation matrix retrieved"
    return $matrix
}

proc get_matrixenergy {filename} {
    set fp [open $filename r]
    set count 0
    while {[gets $fp line]!= -1} {
        incr count
        if {$count==17} {
            set energy $line
            break
            }
    }
    close $fp
    set energy [lindex $energy 1]
    return $energy
}


proc getPropPdb {pdbfile {num 7} {clash_cutoff 0}} {
    mol new ${pdbfile}.pdb
    set matrix [transmatrix "${pdbfile}_matrix.hex"]
    set sel [atomselect top all]
    set n $num
    for {set i 1} {$i <= $n} {incr i} {
        animate dup frame [expr {$i-1}] [molinfo top get id]
        $sel frame $i
        $sel move $matrix
        $sel set segid Seg$i
        $sel set chain $i
        $sel writepdb ${i}_temp.pdb
    }
    $sel delete
    mol delete top
    set outfile [open ${pdbfile}_docked.pdb w+]
    for {set i 1} {$i <= $n} {incr i} {
        set myfile [open ${i}_temp.pdb]
        set data [read $myfile]
        puts $outfile $data
        close $myfile;
        file delete ${i}_temp.pdb
    }
    close $outfile
#    animate write pdb ${pdbfile}_docked.pdb
    assign_ENDtoTER_pdb ${pdbfile}_docked
    mol new ${pdbfile}_docked.pdb
    set clash [get_clash]
    puts "clash value $clash"
    if {$clash < $clash_cutoff } {
        file copy ${pdbfile}_docked.pdb ./prop_str
    }
    mol delete top
    #mol new ${pdbfile}_docked.pdb
}


proc load_dockedfiles { } {
    set docked_pdbs [glob -directory ./ *docked.pdb]
    set outfile [open docked_info.dat w]
    puts $outfile "Frame Energy Gyration"
    foreach pdb $docked_pdbs {
        puts $pdb
        mol new $pdb
        set clash_atoms [measure contacts 1.5 [atomselect top all]]
        set atom_num [llength [lindex $clash_atoms 0]]
        set frame [lindex [split [lindex [split [lindex [split $pdb .] 1] /] 1] _] 0]
        set enr [get_matrixenergy ../${frame}_matrix.hex]
        set gyra [get_gyra]
        puts $outfile "$frame $enr $gyra"
        #if {$atom_num != 0} {
            #mol delete top
            #puts "filename: $pdb \n Clashed atoms = $atom_num"
        #}
    }
    close $outfile
}
