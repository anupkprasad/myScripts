#!/usr/bin/tclsh
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy
# @2023, VMD script for system information, amino seqs

source "$pathVmdScript/clashCheck.tcl"

proc get_cell {{molid top}} {
 set all [atomselect $molid all]
 set minmax [measure minmax $all]
 set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]]
 set a "cellBasisVector1 [lindex $vec 0] 0 0"
 set b "cellBasisVector2 0 [lindex $vec 1] 0"
 set c "cellBasisVector3 0 0 [lindex $vec 2]"
 set center "cellOrigin [measure center $all]"
 $all delete
 return [list $a\n$b\n$c\n$center]
}

proc get_aminoseq {{molid top}} {
	set aa_triple { ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR}
	set aa_single { A C D E F G H I K L M N P Q R S T V W Y }
	set sel [atomselect $molid "protein"]
	set uresid [lsort -unique [$sel get resid]]
	## sorting the prod_list according ASCII
	set uresid [lsort -dictionary $uresid]
	set amino_acids ""
	foreach rid $uresid {
		set sel_rid [atomselect $molid "resid $rid and protein"]
	    set rsnm [lsort -unique [$sel_rid get resname]]
	    #puts $rsnm
	    set aa_in [lsearch -exact $aa_triple $rsnm]
	    append amino_acids [lindex $aa_single $aa_in]
	    $sel_rid delete
	}
	set cter [atomselect $molid "resid [lindex $uresid end]"]
	set cter_atms [$cter get name]
	set is_NH2 "No"
	if {[lsearch -exact $cter_atms "NT"] >= 0} {
    append amino_acids "-NH2"
    set is_NH2 "Yes"
	}
	$sel delete
	return [list $amino_acids $is_NH2]
}



proc sys_info {{molid top} {write_file "No"}} {
	set aa_seqs "Amino acid Sequence: [lindex [get_aminoseq $molid] 0]"
	if {[lindex [get_aminoseq $molid] 1] == "Yes"} {
		set aa_num "Peptide length: [expr [llength [split [lindex [get_aminoseq $molid] 0] ""]] - 4]"
	} else {
		set aa_num "Peptide length: [llength [split [lindex [get_aminoseq $molid] 0] ""]]"
	}
	set cell_dimention "system dimention: [get_cell $molid]"
	set frame "Frame_num: [molinfo top get numframes]"
	set Schar "System charge: [eval "vecadd [[atomselect $molid all] get charge]"]"
	set Pchar "Peptide charge: [eval "vecadd [[atomselect $molid "protein"] get charge]"]"
	set all "all_atoms_num: [[atomselect top all] num]"
	set wat "Water_num: [expr double ([[atomselect top "water"] num])/3]"
	#set TFE "TFE: [expr double ([[atomselect top "resname TFE"] num])/9]"
	#set Wat_TFE "Water:TFE= [expr double ([expr [[atomselect top "water"] num]/3])/[expr [[atomselect top "resname TFE"] num]/9]]"
	set ion "Ions: [[atomselect top "ions"] num]"
 	set all_info [list $aa_seqs\n$aa_num\n$Pchar\n$cell_dimention\n$Schar\n$frame\n$all\n$wat\n$ion]
 	if {$write_file != "No"} {
		set outfile [open system_info.dat w]
		puts $outfile $all_info
		puts $outfile "Total clash_atoms =  [get_clash]"
		close $outfile
 	}
 	
 	return $all_info
}


