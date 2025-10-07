#!/usr/bin/tclsh
# @2022, VMD script for automatic psf genration/solvation/ionisation of system from a pdb
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy

proc get_min { items } {
  set min [lindex $items 0]
  foreach i $items {
    if { [expr abs($i)] < $min } {
      set min [expr abs($i)]
    }
  }
  puts "min is = $min"
  return $min
}

proc get_max { items } {
  set max [lindex $items 0]
  foreach i $items {
    if { [expr abs($i)] > $max } {
      set max [expr abs($i)]
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
	$a moveby [vecscale $c -1.0]; #system moves to orgin
	$a moveby $newcom
	$a writepdb ${inputname}Cent.pdb
}

proc protMoveOnLipSurface {{pdbfile prot_lip.pdb}} {
	mol new $pdbfile
	set prot [atomselect top "protein"]
	moveSysOrigin
	$prot moveby {0 0 33.34}  ; # Z of lip surface is 23.34 and +10 pad
	set sel [atomselect top all]
	$sel writepdb prot_lip.pdb
}

proc rotateMol { {axis x} {angle 5}} {
	set sel [atomselect top "protein"]
	set com [measure center $sel weight mass]
	set matrix [transaxis $axis $angle]
	$sel moveby [vecscale -1.0 $com]
	$sel move $matrix
	$sel moveby $com
}


#set sel [atomselect top "protein"]
#set com [measure center $sel weight mass]
#set matrix [transaxis z 90]
#$sel moveby [vecscale -1.0 $com]
#$sel move $matrix
#$sel moveby $com

# Function to calculate the distance between two points
proc calculate_distance {point1 point2} {
    # Extract the coordinates of the first point
    set x1 [lindex $point1 0]
    set y1 [lindex $point1 1]
    set z1 [lindex $point1 2]
    
    # Extract the coordinates of the second point
    set x2 [lindex $point2 0]
    set y2 [lindex $point2 1]
    set z2 [lindex $point2 2]
    
    # Calculate the differences in each coordinate
    set dx [expr {$x2 - $x1}]
    set dy [expr {$y2 - $y1}]
    set dz [expr {$z2 - $z1}]
    
    # Calculate the distance using the Euclidean distance formula
    set distance [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]
    
    return $distance
}



proc getBoxSetup {{molid top} {pad 7} {cube_padsize yes} } {
	set a [atomselect $molid all]
	set c [measure center $a]
	set box { { ? ? ? } { ? ? ? } }
	set basisvec { ? ? ? }
	set origin { ? ? ? }	
	set minmax [measure minmax $a]
	
	if {$cube_padsize == yes} {
		#earlier I used following two lines but now modified with next three lines of it
		#set minmax_sub [vecsub [lindex $minmax 0] [lindex $minmax 1]]
		#set cube_halfside [expr 0.5*([get_max $minmax_sub]) + $pad]
		set j [lindex $minmax 0]
		set k [lindex $minmax 1]
		set distance [calculate_distance $j $k]
		set cube_halfside [expr 0.5*([get_max $distance]) + $pad]
		puts "box size = $cube_halfside"
		foreach d {0 1 2} {
	  	lset box 0 $d [expr [lindex $c $d] - $cube_halfside]
	  	lset box 1 $d [expr [lindex $c $d] + $cube_halfside]
	  	lset basisvec $d [expr [lindex $box 1 $d ] - [lindex $box 0 $d]] 
	  	lset origin $d [expr 0.5*([lindex $box 1 $d ] + [lindex $box 0 $d])] 
		}
	} else {
		foreach d {0 1 2} {
	  	lset box 0 $d [expr [lindex $minmax 0 $d] - $pad]
	  	lset box 1 $d [expr [lindex $minmax 1 $d] + $pad]
	  	lset basisvec $d [expr [lindex $box 1 $d ] - [lindex $box 0 $d]] 
	  	lset origin $d [expr 0.5*([lindex $box 1 $d ] + [lindex $box 0 $d])] 
	  }

	}
	puts "Box size: $box"
	return $box
	mol delete top
}


proc getBoxProtLip {pdbfile {pad 7} } {
	set inputname $pdbfile
	mol new ${inputname}.pdb
	set a [atomselect top "segname XO1"]	
	set box [measure minmax $a]
	set zmin [lindex $box 0 2]
	set boxlower [measure minmax $a]
	lset boxlower 0 2 [expr $zmin - $pad]
	lset boxlower 1 2 $zmin
	set zmax [lindex $box 1 2]
	set boxupper [measure minmax $a]
	lset boxupper 0 2 $zmax
	lset boxupper 1 2 [expr $zmax + $pad]
	puts "boxupper size: $boxupper"
	puts "boxlower size: $boxlower"
	set boxboth [list $boxupper $boxlower]
	return $boxboth
	mol delete top
}


proc getpsf {pdbfile {amidated "No"} } {
	## generate pdb and psf files with CTER amidation of given input pdb 
	set inputname $pdbfile
	resetpsf
	mol new $inputname.pdb
	set sel [atomselect top "protein"]
	set cter_index [lindex [$sel get resid] end]; ### get resid of cter residue
	topology "/usr/local/lib/vmd/plugins/noarch/tcl/readcharmmtop1.2/top_all36_prot.rtf"
	segment P1 [list pdb $inputname.pdb]
	if {$amidated != "No"} {
		patch CT2 P1:$cter_index  ; ##CTER amidations
	}
	regenerate angles dihedrals
	coordpdb $inputname.pdb P1
	guesscoord
	writepdb ${inputname}p.pdb
	writepsf ${inputname}p.psf
	mol new ${inputname}p.psf
	mol addfile ${inputname}p.pdb
	puts "\n \n \n psf finished \n \n \n"
}


proc getsolv {input {box} } {
	#vmdcon -info "Usage: solvate <psffile> <pdbfile> <option1?> <option2?>..."
	#vmdcon -info "Usage: solvate <option1?> <option2?>...  to just create a water box" 
	#vmdcon -info "Options:"
	#vmdcon -info "    -o <output prefix> (data will be written to output.psf/output.pdb)"
	#vmdcon -info "    -s <segid prefix> (should be either one or two letters; default WT)"
	#vmdcon -info "    -b <boundary> (minimum distance between water and solute, default 2.4)"
	#vmdcon -info "    -minmax {{xmin ymin zmin} {xmax ymax zmax}}"
	#vmdcon -info "    -rotate (rotate molecule to minimize water volume)"
	#vmdcon -info "    -rotsel <selection> (selection of atoms to check for rotation)"
	#vmdcon -info "    -rotinc <increment> (degree increment for rotation)"
	#vmdcon -info "    -t <pad in all directions> (override with any of the following)"
	#vmdcon -info "    -x <pad negative x>"
	#vmdcon -info "    -y <pad negative y>"
	#vmdcon -info "    -z <pad negative z>"
	#vmdcon -info "    +x <pad positive x>"
	#vmdcon -info "    +y <pad positive y>"
	#vmdcon -info "    +z <pad positive z>"
	#vmdcon -info "    The following options allow the use of solvent other than water:"
	#vmdcon -info "      -spsf <solventpsf> (PSF file for nonstandard solvent)"
	#vmdcon -info "      -spdb <solventpdb> (PDB file for nonstandard solvent)"
	#vmdcon -info "      -stop <solventtop> (Topology file for nonstandard solvent)"
	#vmdcon -info "      -ws <size> (Box length for nonstandard solvent)"
	#vmdcon -info "      -ks <keyatom> (Atom occuring once per residue for nonstandard solvent)"
	puts "\n \n \n ############  solvation start...    ############## \n \n \n"
	package require solvate
	solvate ${input}.psf ${input}.pdb -minmax $box -o solvate

}

proc getionze {input} {
	#puts "Usage: autoionize -psf file.psf -pdb file.pdb <mode> \[options\]"
	#puts "Ion placement mode (choose one):"
	#puts "  -neutralize              -- only neutralize system"
	#puts "  -sc <salt concentration> -- neutralize and set salt concentration (mol/L)"
	## puts "  -is <ionic strength>     -- neutralize and set ionic strength (mol/L)"
	#puts "  -nions {{ion1 num1} {ion2 num2} ...} -- user defined number of ions"
	#puts "Other options:"
	#puts "  -cation <ion resname>    -- default: SOD"
	#puts "  -anion <ion resname>     -- default: CLA"
	#puts "  -o <prefix>              -- output file prefix (default: ionized)"
	#puts "  -from <distance>         -- min. distance from solute (default: 5A)"
	#puts "  -between <distance>      -- min. distance between ions (default: 5A)"
	#puts "  -seg <segname>           -- specify new segment name (default: ION)"
	#puts "Supported ions (CHARMM force field resnames):#
	puts "\n \n \n ######### autoionize start... ########### \n \n \n"
	package require autoionize
	autoionize -psf ${input}.psf -pdb ${input}.pdb -sc 0.15 -o ionized
}



proc constrinedFixedFiles {{molid top}} {
	set all [atomselect top all] 
	set p [atomselect top "protein and backbone"] 
	$all set beta 0 
	$p set beta 1
	$all writepdb ./fixed.pdb 

	set all [atomselect top all] 
	set q [atomselect top "protein and name CA"] 
	$all set beta 0 
	$q set beta 0.5
	$all writepdb ./constrained.pdb
}


proc restraintRefFiles {{molid top}} {
	### for protein backbone
	set all [atomselect top all] 
	set p [atomselect top "protein and backbone"] 
	$all set beta 0 
	$p set beta 1
	$all writepdb ../restraints/prot_posres.ref
	#### for given POPC:POPG = 4:1 (generated file)
	## for lower head POPC (atoms = N C12 C13 C14 C15 C11 P O13 O14 O12 O11)
	set all [atomselect top "all and not water"]
	set popc_UH [atomselect top "segname XO1 and not resid > 100 and resname POPC and name N C12 C13 C14 C15 C11 P O13 O14 O12 O11"]
	$all set beta 0
	$popc_UH set beta 1
	$all writepdb ../restraints/popc_head_upper.ref
	set popc_LH [atomselect top "segname XO1 and resid > 100 and resname POPC and name N C12 C13 C14 C15 C11 P O13 O14 O12 O11"]
	$all set beta 0
	$popc_LH set beta 1
	$all writepdb ../restraints/popc_head_lower.ref
	## for lower head POPG (atoms = C13 OC3 C12 OC2 C11 P O13 O14 O12 O11)
	set all [atomselect top "all and not water"]
	set popg_UH [atomselect top "segname XO1 and not resid > 100 and resname POPC and name C13 OC3 C12 OC2 C11 P O13 O14 O12 O11"]
	$all set beta 0
	$popg_UH set beta 1
	$all writepdb ../restraints/popg_head_upper.ref
	set popg_LH [atomselect top "segname XO1 and resid > 100 and resname POPC and name C13 OC3 C12 OC2 C11 P O13 O14 O12 O11"]
	$all set beta 0
	$popg_LH set beta 1
	$all writepdb ../restraints/popg_head_lower.ref

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



# Function to insert lines into a file at a specific line number
proc insert_lines_in_file {filename line_number new_lines} {
    # Open the file for reading
    set in [open $filename r]
    
    # Read all lines into a list
    set lines [split [read $in] "\n"]
    close $in
   
    # Insert the new lines into the list of lines
    #set new_lines_list [split $new_lines "\n"]
    set lines [linsert $lines $line_number $new_lines]
    
    # Open the file for writing and write the modified content
    set out [open $filename w]
    puts $out [join $lines "\n"]
    close $out
}


########### MD setup functions ###########
proc sysSetup {{molid top} {cter_amidation "No"}} {
	global pathSetup
	set pdbfile [molinfo $molid get name]
	set inputname [lindex [split $pdbfile .] 0]
	#getpsf ${pdbfile}_c $cter_amidation
	set box [getBoxSetup $molid 7]
	getsolv ${inputname} $box
	getionze solvate
	constrinedFixedFiles

	exec cp -r $pathSetup/protein/equilibration ./
	set pbcinfo [lindex [getCellDimention] 0]
	insert_lines_in_file ./equilibration/equilibration.inp 50 $pbcinfo

	exec cp -r $pathSetup/protein/prod1 ./
	exec cp -r $pathSetup/protein/parameter ./

	getSysInfo top yes
}


proc sysSetupSDS {{molid top} {cter_amidation "No"}} {
	global pathSetup
	set pdbfile [molinfo $molid get name]
	set inputname [lindex [split $pdbfile .] 0]
	#getpsf ${pdbfile}_c $cter_amidation
	set box [getBoxSetup $molid 7]
	getsolv ${inputname} $box
	getionze solvate
	constrinedFixedFiles

	exec cp -r $pathSetup/proteinSDS/equilibration ./
	set pbcinfo [lindex [getCellDimention] 0]
	insert_lines_in_file ./equilibration/equilibration.inp 50 $pbcinfo

	exec cp -r $pathSetup/proteinSDS/prod1 ./
	exec cp -r $pathSetup/proteinSDS/parameter ./
	set sds [atomselect top "segname XO1"]
	set min [expr [get_min [$sds get index]] + 1]
	set max [expr [get_max [$sds get index]] + 1]
	set line "			atomnumbersrange {$min-$max}"
	insert_lines_in_file ./parameter/com_colvar.in 11 $line

	sys_info top yes
}


proc sysSetupLipbilayr {pdbfile {cter_amidation "No"}} {
	global pathSetup
	###pwd == equilbration
	#getpsf ${pdbfile}_c $cter_amidation
	set box [getBoxProtLip ${pdbfile}_cp 20]
	puts "box1: [lindex $box 0]"
	puts "box2: [lindex $box 1]"
	getsolv ${pdbfile}_cp [lindex $box 0]
	getsolv ${pdbfile}_cpw [lindex $box 1]
	getionze ${pdbfile}_cpww

	file mkdir -p ../restraints
	restraintRefFiles
	file copy -force $pathSetup/refsParaInp/lipProt/lipColvarTemplate.dat ../restraints

	file mkdir -p ../equilbration
	file copy -force $pathSetup/refsParaInp/namd.sh ../equilbration
	for { set a 1}  {$a < 7} {incr a} {
		file copy -force $pathSetup/refsParaInp/lipProt/step6.${a}_equilibration.inp ../equilbration
	}

	file copy -force $pathSetup/refsParaInp/toppar ../
	sysinfo_inp
	sys_info top yes

}


