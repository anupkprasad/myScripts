#!/usr/bin/tclsh
# @2022, VMD script for psf genration, solvation and ionisation of system from a pdb
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy

proc getmax { items } {
  set max 1
  foreach i $items {
    if { [expr abs($i)] > $max } {
      set max [expr abs($i)]
    }
  }
  puts "max is = $max"
  return $max
}


proc getBoxofCenteredSystem {{molid top} {pad 10} {cube_padsize "No"} } {
	set pdbfile [molinfo $molid get name]
	set inputname [lindex [split $pdbfile .] 0]
	set a [atomselect $molid all]
	set c [measure center $a]
	
	$a moveby [vecscale $c -1.0]; #system moves to orgin
	set new_c [measure center $a]
	set box { { ? ? ? } { ? ? ? } }
	set basisvec { ? ? ? }
	set origin { ? ? ? }	
	set minmax [measure minmax $a]
	
	if {$cube_padsize != "No"} {
		set minmax_sub [vecsub [lindex $minmax 0] [lindex $minmax 1]]
		set cube_halfside [expr 0.5*([getmax $minmax_sub]) + $pad]
		puts "box size = $cube_halfside"
		foreach d {0 1 2} {
	  	lset box 0 $d [expr [lindex $new_c $d] - $cube_halfside]
	  	lset box 1 $d [expr [lindex $new_c $d] + $cube_halfside]
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

	$a writepdb ${inputname}_c.pdb
	puts "recenter finish"
	return $box
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
	solvate ${input}.psf ${input}.pdb -minmax $box -o ${input}w

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
	autoionize -psf ${input}.psf -pdb ${input}.pdb -sc 0.15 -o ${input}i
}



proc constrined_fixed_files {{molid top}} {
	set all [atomselect top all] 
	set p [atomselect top "protein and backbone"] 
	$all set beta 0 
	$p set beta 1
	$all writepdb fixed.pdb 

	set all [atomselect top all] 
	set q [atomselect top "protein and name CA"] 
	$all set beta 0 
	$q set beta 0.5
	$all writepdb constrained.pdb
}


########### code initializations  ###########


proc md_setup {pdbfile {cter_amidation "No"}} {
	set box [getBoxofCenteredSystem top 7]
	#getpsf ${pdbfile}_c $cter_amidation
	getsolv ${pdbfile}_cp $box
	getionze ${pdbfile}_cpw
	constrined_fixed_files
	sys_info top yes
}
