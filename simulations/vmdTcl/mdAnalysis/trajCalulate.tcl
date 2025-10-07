#!/usr/bin/tclsh
#### Trajectory analysis functions for VMD
# Author: Anup K. Prasad
# Email: anupkprasad121@gmail.com
# Description: Collection of functions for RMSD, RMSF, and secondary structure analysis

proc calcSecondaryStructure {{molid top} {selection "protein and name CA"} {outfile "ss_traj.dat"}} {
    # Calculate secondary structure of protein atoms for each frame
    # molid: molecule ID (default: top)
    # selection: atom selection (default: "protein and name CA")
    # outfile: output filename (default: "ss_traj.dat")
    
    set sel [atomselect $molid $selection]
    set num_frames [molinfo $molid get numframes]
    set output_handle [open $outfile w]
    
    puts "Calculating secondary structure for $num_frames frames..."
    puts $output_handle "# Frame-by-frame secondary structure data"
    puts $output_handle "# Each line represents one frame"
    
    display update off
    
    for {set frame 0} {$frame < $num_frames} {incr frame} {
        animate goto $frame
        $sel frame $frame
        mol ssrecalc $molid
        set secstruct [$sel get structure]
        puts $output_handle "$secstruct"
        
        # Progress update every 100 frames
        if {[expr $frame % 100] == 0 || $frame == [expr $num_frames-1]} {
            puts "Frame $frame processed"
        }
    }
    
    display update on
    close $output_handle
    $sel delete
    
    puts "Secondary structure calculation completed. Results saved to: $outfile"
}


proc calcRmsdReference {{molid_ref top} {sel_ref "protein and name CA"} {molid_traj top} {sel_traj "protein and name CA"} {outfile "rmsd.dat"}} {
    # Calculate RMSD of trajectory against reference structure
    # molid_ref: molecule ID of reference structure
    # sel_ref: selection string for reference structure
    # molid_traj: molecule ID of trajectory
    # sel_traj: selection string for trajectory
    # outfile: output filename (default: rmsd.dat)
    
    # Create reference selection (frame 0 of reference molecule)
    set reference [atomselect $molid_ref $sel_ref frame 0]
    
    # Create trajectory selection
    set compare [atomselect $molid_traj $sel_traj]
    
    # Check atom count compatibility
    set ref_atoms [$reference num]
    set traj_atoms [$compare num]
    puts "Reference atoms: $ref_atoms, Trajectory atoms: $traj_atoms"
    
    if {$ref_atoms != $traj_atoms} {
        puts "Error: Atom count mismatch between reference and trajectory selections!"
        $reference delete
        $compare delete
        return
    }
    
    # Open output file
    set outfile_handle [open $outfile w]
    puts $outfile_handle "# Frame RMSD(Angstrom)"
    puts $outfile_handle "# Reference: molid $molid_ref, selection: $sel_ref"
    puts $outfile_handle "# Trajectory: molid $molid_traj, selection: $sel_traj"
    
    # Get number of frames in trajectory
    set num_frames [molinfo $molid_traj get numframes]
    puts "Calculating RMSD for $num_frames frames..."
    
    # Loop through all trajectory frames
    for {set frame 0} {$frame < $num_frames} {incr frame} {
        # Set trajectory selection to current frame
        $compare frame $frame
        
        # Compute the transformation matrix (align trajectory to reference)
        set trans_mat [measure fit $compare $reference]
        
        # Apply alignment to trajectory selection
        $compare move $trans_mat
        
        # Compute the RMSD after alignment
        set rmsd [measure rmsd $compare $reference]
        
        # Write to output file
        puts $outfile_handle [format "%d %.6f" $frame $rmsd]
        
        # Print progress every 100 frames
        if {[expr $frame % 100] == 0 || $frame == [expr $num_frames-1]} {
            puts "Frame $frame: RMSD = [format "%.3f" $rmsd] Å"
        }
    }
    
    # Clean up
    close $outfile_handle
    $reference delete
    $compare delete
    
    puts "RMSD calculation completed. Results saved to: $outfile"
}


proc calcRmsfReference {{molid_ref top} {sel_ref "protein and name CA"} {molid_traj top} {sel_traj "protein and name CA"} {outfile "rmsf.dat"} {pdbfile ""}} {
    # Calculate Root Mean Square Fluctuation (RMSF) against reference structure
    # molid_ref: molecule ID of reference structure (default: top)
    # sel_ref: selection string for reference structure (default: "protein and name CA")
    # molid_traj: molecule ID of trajectory (default: top)
    # sel_traj: selection string for trajectory (default: "protein and name CA")
    # outfile: output filename (default: "rmsf.dat")
    # pdbfile: optional PDB output filename with RMSF as B-factor and RMSF std dev as occupancy (default: "")
    # Returns: list of {rmsf_values rmsf_std_deviations}
    
    # Create atom selections
    set refsel [atomselect $molid_ref $sel_ref frame 0]
    set trajsel [atomselect $molid_traj $sel_traj]

    # Check selection compatibility
    set num_atoms [$refsel num]
    set traj_atoms [$trajsel num]
    
    if {$num_atoms != $traj_atoms} {
        puts "Error: Selection mismatch - Reference: $num_atoms atoms, Trajectory: $traj_atoms atoms"
        $refsel delete
        $trajsel delete
        return
    }

    # Get reference coordinates and initialize arrays
    set ref_coords [$refsel get {x y z}]
    set num_frames [molinfo $molid_traj get numframes]
    set sqsum [lrepeat $num_atoms 0.0]
    set sq_sqsum [lrepeat $num_atoms 0.0]  ;# For standard deviation calculation

    puts "Calculating RMSF for $num_atoms atoms over $num_frames frames..."

    # Calculate squared deviations for each frame
    for {set frame 0} {$frame < $num_frames} {incr frame} {
        $trajsel frame $frame
        set traj_coords [$trajsel get {x y z}]
        
        # Calculate per-atom squared deviations
        for {set atom 0} {$atom < $num_atoms} {incr atom} {
            set ref_atom [lindex $ref_coords $atom]
            set traj_atom [lindex $traj_coords $atom]
            
            set dx [expr {[lindex $traj_atom 0] - [lindex $ref_atom 0]}]
            set dy [expr {[lindex $traj_atom 1] - [lindex $ref_atom 1]}]
            set dz [expr {[lindex $traj_atom 2] - [lindex $ref_atom 2]}]
            set sqdist [expr {$dx*$dx + $dy*$dy + $dz*$dz}]
            
            # Accumulate for mean and variance calculations
            lset sqsum $atom [expr {[lindex $sqsum $atom] + $sqdist}]
            lset sq_sqsum $atom [expr {[lindex $sq_sqsum $atom] + $sqdist*$sqdist}]
        }
        
        # Progress update every 100 frames
        if {[expr $frame % 100] == 0 || $frame == [expr $num_frames-1]} {
            puts "Processed frame $frame"
        }
    }

    # Calculate RMSF and standard deviation
    set rmsf_values {}
    set rmsf_std_devs {}
    set output_file ""
    
    if {$outfile ne ""} {
        set output_file [open $outfile w]
        puts $output_file "# Atom_Index RMSF(Angstrom) RMSF_StdDev(Angstrom)"
        puts $output_file "# Reference: molid $molid_ref, selection: $sel_ref"
        puts $output_file "# Trajectory: molid $molid_traj, selection: $sel_traj"
        puts $output_file "# Frames analyzed: $num_frames"
    }

    for {set atom 0} {$atom < $num_atoms} {incr atom} {
        # Calculate RMSF (root mean square fluctuation)
        set mean_sq [expr {[lindex $sqsum $atom] / double($num_frames)}]
        set rmsf_val [expr {sqrt($mean_sq)}]
        
        # Calculate standard deviation of squared distances
        set mean_sq_sq [expr {[lindex $sq_sqsum $atom] / double($num_frames)}]
        set variance [expr {$mean_sq_sq - $mean_sq*$mean_sq}]
        set std_dev [expr {$variance > 0 ? sqrt($variance) : 0.0}]
        
        lappend rmsf_values $rmsf_val
        lappend rmsf_std_devs $std_dev
        
        if {$output_file ne ""} {
            puts $output_file [format "%d %.6f %.6f" $atom $rmsf_val $std_dev]
        }
    }

    if {$output_file ne ""} {
        close $output_file
        puts "RMSF data written to: $outfile"
    }

    # Write PDB file with RMSF as B-factor and RMSF standard deviation as occupancy if requested
    if {$pdbfile ne ""} {
        # Create a new selection for the reference structure to modify B-factors and occupancy
        set pdb_sel [atomselect $molid_ref $sel_ref frame 0]
        
        # Set RMSF values as B-factors and RMSF standard deviations as occupancy
        $pdb_sel set beta $rmsf_values
        $pdb_sel set occupancy $rmsf_std_devs
        
        # Write PDB file
        $pdb_sel writepdb $pdbfile
        
        puts "PDB file with RMSF data written to: $pdbfile"
        puts "B-factors contain per-atom RMSF values"
        puts "Occupancy contains per-atom RMSF standard deviations"
        
        # Clean up PDB selection
        $pdb_sel delete
    }

    # Clean up
    $refsel delete
    $trajsel delete
    
    puts "RMSF calculation completed."
    return [list $rmsf_values $rmsf_std_devs]
}



# Procedure to analyze contacts between two domains across trajectory frames
# Arguments:
#   molid   : molecule id
#   selA    : atom selection for domain A
#   selB    : atom selection for domain B
#   cutoff  : contact cutoff in Å (default = 4.0 Å)
#   outfile : output filename (optional)

proc contactAnalysis {molid selA selB {cutoff 4.0} {outfile ""}} {
    set nframes [molinfo $molid get numframes]
    puts "Analyzing contacts for $nframes frames"

    # Open output file if specified
    set output_file ""
    if {$outfile ne ""} {
        set output_file [open $outfile w]
        puts $output_file "# Frame, Contacts, Dom_A_Residues, Dom_B_Residues"
    }

    # Store total contact counts for frequency
    set total_contacts 0
    set total_possible 0

    for {set f 0} {$f < $nframes} {incr f} {
        # Select domains for this frame
        set selAframe [atomselect $molid $selA frame $f]
        set selBframe [atomselect $molid $selB frame $f]

        # --- Contact calculation ---
        set contacts [measure contacts $cutoff $selAframe $selBframe]
        set num_contacts [llength [lindex $contacts 0]]
        
        # Get unique residue IDs that make contacts
        set contact_atoms_A [lindex $contacts 0]
        set contact_atoms_B [lindex $contacts 1]
        set unique_resids_A {}
        set unique_resids_B {}
        set resinfo_A {}
        set resinfo_B {}
        
        if {$num_contacts > 0} {
            # Get residue IDs and names for contacting atoms in domain A
            foreach atom_idx $contact_atoms_A {
                set atom_sel [atomselect $molid "index $atom_idx" frame $f]
                set resid [$atom_sel get resid]
                set resname [$atom_sel get resname]
                $atom_sel delete
                if {[lsearch -exact $unique_resids_A $resid] == -1} {
                    lappend unique_resids_A $resid
                    lappend resinfo_A "$resid $resname"
                }
            }
            
            # Get residue IDs and names for contacting atoms in domain B
            foreach atom_idx $contact_atoms_B {
                set atom_sel [atomselect $molid "index $atom_idx" frame $f]
                set resid [$atom_sel get resid]
                set resname [$atom_sel get resname]
                $atom_sel delete
                if {[lsearch -exact $unique_resids_B $resid] == -1} {
                    lappend unique_resids_B $resid
                    lappend resinfo_B "$resid $resname"
                }
            }
            
            # Sort the residue lists for consistent output
            set unique_resids_A [lsort -integer $unique_resids_A]
            set unique_resids_B [lsort -integer $unique_resids_B]
            
            # Sort resinfo lists by residue number (simple approach using the sorted unique_resids)
            set sorted_resinfo_A {}
            set sorted_resinfo_B {}
            
            foreach resid $unique_resids_A {
                foreach resinfo $resinfo_A {
                    if {[lindex [split $resinfo] 0] == $resid} {
                        lappend sorted_resinfo_A $resinfo
                        break
                    }
                }
            }
            
            foreach resid $unique_resids_B {
                foreach resinfo $resinfo_B {
                    if {[lindex [split $resinfo] 0] == $resid} {
                        lappend sorted_resinfo_B $resinfo
                        break
                    }
                }
            }
            
            set resinfo_A $sorted_resinfo_A
            set resinfo_B $sorted_resinfo_B
        }

        # Update totals
        incr total_contacts $num_contacts
        incr total_possible [$selAframe num]

        # Print contact information for frames with contacts
        if {$num_contacts > 0} {
            puts "Frame $f: Contacts = $num_contacts"
            puts "  Contact residues - SelA ([join $resinfo_A ", "]), SelB ([join $resinfo_B ", "])"
        }
        
        # Write to output file if specified
        if {$output_file ne ""} {
            if {$num_contacts > 0} {
                puts $output_file [format "%d, %d, SelA (%s), SelB (%s)" $f $num_contacts [join $resinfo_A ", "] [join $resinfo_B ", "]]
            } else {
                puts $output_file [format "%d, %d," $f $num_contacts]
            }
        }

        $selAframe delete
        $selBframe delete
    }

    # Average contact frequency
    set freq [expr {double($total_contacts) / $total_possible}]
    puts "Average contact frequency = $freq"
    
    # Close output file and add summary
    if {$output_file ne ""} {
        puts $output_file "# Average contact frequency = $freq"
        close $output_file
        puts "Contact analysis results written to: $outfile"
    }
}


proc calcCenterOfMassMultiple {molid selection_list {outfile "com_trajectory.dat"}} {
    # Calculate center of mass for multiple selections across all frames
    # molid: molecule ID (default: top)
    # selection_list: list of selection strings (e.g., {"resid 5 7 9" "resname ADE"})
    #                 OR single selection string (e.g., "resid 45 and name CA")
    # outfile: output filename (default: "com_trajectory.dat")
    
    set num_frames [molinfo $molid get numframes]
    
    # Handle both single string and list input
    # Check if it's already a proper list (has curly braces or is a variable containing a list)
    if {[catch {llength $selection_list} list_length]} {
        # If llength fails, treat as single string
        set selection_list [list $selection_list]
        set num_selections 1
    } else {
        # Check if it's a single string that got split by spaces
        if {$list_length > 1} {
            # Try to create atomselect with first element to see if it's valid
            if {[catch {atomselect $molid [lindex $selection_list 0]} test_sel]} {
                # First element is not a valid selection, so treat whole thing as single string
                set selection_list [list $selection_list]
                set num_selections 1
            } else {
                # First element is valid, so it's likely a proper list
                $test_sel delete
                set num_selections $list_length
            }
        } else {
            set num_selections $list_length
        }
    }
    
    if {$num_selections == 0} {
        puts "Error: No selections provided"
        return
    }
    
    puts "Calculating center of mass for $num_selections selections over $num_frames frames..."
    
    # Create atom selections for all groups
    set sel_objects {}
    set sel_names {}
    
    for {set i 0} {$i < $num_selections} {incr i} {
        set sel_string [lindex $selection_list $i]
        set sel_obj [atomselect $molid $sel_string]
        
        # Check if selection is valid and has atoms
        if {[$sel_obj num] == 0} {
            puts "Warning: Selection '$sel_string' contains no atoms"
            $sel_obj delete
            continue
        }
        
        lappend sel_objects $sel_obj
        lappend sel_names $sel_string
        puts "Selection [expr $i + 1]: '$sel_string' - [$sel_obj num] atoms"
    }
    
    set valid_selections [llength $sel_objects]
    if {$valid_selections == 0} {
        puts "Error: No valid selections found"
        return
    }
    
    # Open output file
    set output_file [open $outfile w]
    
    # Write header
    puts $output_file "# Center of mass trajectory data"
    puts $output_file "# Molecule ID: $molid"
    puts $output_file "# Frames: $num_frames"
    puts $output_file "# Selections:"
    for {set i 0} {$i < $valid_selections} {incr i} {
        puts $output_file "#   Selection [expr $i + 1]: [lindex $sel_names $i]"
    }
    
    # Create column headers (comma-separated)
    set header "Frame"
    for {set i 0} {$i < $valid_selections} {incr i} {
        set header "$header, Sel[expr $i + 1]_X, Sel[expr $i + 1]_Y, Sel[expr $i + 1]_Z"
    }
    puts $output_file $header
    
    # Calculate center of mass for each frame
    for {set frame 0} {$frame < $num_frames} {incr frame} {
        set line_data "$frame"
        
        # Update all selections to current frame and calculate COM
        for {set i 0} {$i < $valid_selections} {incr i} {
            set sel_obj [lindex $sel_objects $i]
            $sel_obj frame $frame
            
            # Calculate center of mass
            set com [measure center $sel_obj weight mass]
            set x [format "%.4f" [lindex $com 0]]
            set y [format "%.4f" [lindex $com 1]]
            set z [format "%.4f" [lindex $com 2]]
            
            set line_data "$line_data, $x, $y, $z"
        }
        
        puts $output_file $line_data
        
        # Progress update every 100 frames
        if {[expr $frame % 100] == 0 || $frame == [expr $num_frames-1]} {
            puts "Processed frame $frame"
        }
    }
    
    # Clean up
    foreach sel_obj $sel_objects {
        $sel_obj delete
    }
    
    close $output_file
    puts "Center of mass calculation completed. Results saved to: $outfile"
    
    return $valid_selections
}

# ========================================================================
# USAGE EXAMPLES
# ========================================================================

# SECONDARY STRUCTURE ANALYSIS:
# Calculate secondary structure for all frames
# calcSecondaryStructure top "protein and name CA" "secondary_structure.dat"

# RMSD CALCULATIONS:
# Calculate RMSD of trajectory against reference structure
# calcRmsdReference 0 "protein and backbone" 1 "protein and backbone" "rmsd_backbone.dat"
#
# Calculate RMSD of CA atoms
# calcRmsdReference 0 "protein and name CA" 1 "protein and name CA" "rmsd_ca.dat"
#
# Calculate RMSD for specific residue range
# calcRmsdReference 0 "resid 100 to 200 and backbone" 1 "resid 100 to 200 and backbone" "rmsd_domain.dat"

# RMSF CALCULATIONS:
# Calculate RMSF of CA atoms comparing trajectory to reference structure
# set result [calcRmsfReference 0 "protein and name CA" 1 "protein and name CA" "rmsf_ca.dat"]
# set rmsf_values [lindex $result 0]
# set rmsf_std_devs [lindex $result 1]
#
# Calculate RMSF for backbone atoms of specific residue range
# set result [calcRmsfReference 0 "protein and backbone and resid 100 to 200" 1 "protein and backbone and resid 100 to 200" "rmsf_backbone.dat"]
#
# Calculate RMSF without saving to file (returns data only)
# set result [calcRmsfReference 0 "protein and name CA" 1 "protein and name CA" ""]


# CONTACT ANALYSIS:
# Example usage for contact analysis:
# contactAnalysis top "resid 54 to 311" "resid 437 to 623" 4.5 "contact_analysis.dat"

# CENTER OF MASS CALCULATIONS:
# Calculate center of mass for multiple selections (using list)
# set selections {"resid 5 7 9" "resname ADE" "protein and name CA and resid 100 to 200"}
# calcCenterOfMassMultiple top $selections "com_multiple.dat"
#
# Calculate COM for specific residues and ligands (using list)
# set my_selections {"resid 25 50 75" "resname ATP" "resname MG"}
# calcCenterOfMassMultiple 0 $my_selections "ligand_protein_com.dat"
#
# Single selection examples:
# calcCenterOfMassMultiple top "resid 45 and name CA" "single_com.dat"
# calcCenterOfMassMultiple top {"protein and name CA"} "protein_com.dat"
#
# Multiple selections without using a variable:
# calcCenterOfMassMultiple top {"resid 45" "resname ATP" "resname MG"} "multi_com.dat"


puts "Trajectory analysis functions loaded successfully!"
puts "Available functions:"
puts "  - calcSecondaryStructure: Calculate secondary structure for each frame"
puts "  - calcRmsdReference: Calculate RMSD against reference structure"
puts "  - calcRmsfReference: Calculate RMSF against reference structure"
puts "  - calcCenterOfMassMultiple: Calculate center of mass for multiple selections"
puts "  - contactAnalysis: Analyze contacts between two domains"
puts "  - calcCenterOfMassMultiple: Calculate center of mass for multiple selections"
puts ""
puts "Use 'info proc calc*' to see all available functions"

