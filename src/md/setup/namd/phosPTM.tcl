#!/usr/bin/tclsh
# @2025, VMD script to apply constraints
# Anup K. Prasad
# anupkprasad121@gmail.com

### list all resids of phophorylated aminos
proc getPhosResidues {{molid top}} {
	### THP1, TPO == phosphothreonine and SP2 and SEP == Phosphoserine
	set selP [atomselect $molid "protein and resname THP1 TPO SP2 SEP PTR"]
	set residsP [lsort -unique [$selP get resid]]
	foreach res $residsP {
		set selresid [atomselect $molid "protein and resid $res"]
		puts "$res : [lsort -unique [$selresid get resname]]" }
	return
}

proc convertPhosphoToNative {{molid top} resids {outfile convertPhosphoToNative.pdb}} {
	foreach res $resids {
		set selres [atomselect $molid "protein and resid $res"]
		set resnameP [lsort -unique [$selres get resname]]
		if {$resnameP == "THP1"} {
			$selres set resname "THR"  ;# Change THP1 to THR
		} elseif {$resnameP == "SP2"} {
			$selres set resname "SER"  ;
		}
		set phos [atomselect $molid "resid $res and name P O1P O2P OT"]
		$phos delete
		## it does not delete actually so save by the follwoiing selection
		puts "modified $res $resnameP"
	}
	set all [atomselect $molid "protein and not {resid $resids and name P O1P O2P OT}"]
	$all writepdb $outfile
}

proc renamePhosphoResidues {{molid top} resids {outfile phosResnameModify.pdb}} {

	foreach res $resids {
		set selres [atomselect $molid "protein and resid $res"]
		set resnameP [lsort -unique [$selres get resname]]
		if {$resnameP == "TPO"} {
			$selres set resname "THP1"
		} elseif {$resnameP == "SEP"} {
			$selres set resname "SP2"  ;
		}

		set phos [atomselect $molid "resid $res and name O3P"]
		$phos set name "OT"
		puts "modified $res $resnameP"
	}
	set all [atomselect $molid all]
	$all writepdb $outfile
}
