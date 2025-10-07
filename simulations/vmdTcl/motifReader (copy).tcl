# -*- coding: utf-8 -*-
# motifReader.tcl
# Created on 2025-09-18
# @author: Anup K. Prasad
# email: anupkprasad121@gmail.com
# 
# TCL version of motif_reader.py functionality for VMD
# Parses MMA format files and extracts motif information
# Returns data as arrays/lists without saving files
# Compatible with VMD Tcl environment

# Parse MMA format file and return sequence data as nested lists
proc parseMMAFile {mma_file} {
    set sequences {}
    
    # Read file content
    set fp [open $mma_file r]
    set content [read $fp]
    close $fp
    
    # Split by sequence entries (each starts with $)
    set entries [split $content "\$"]
    
    # Process each entry (skip first empty part)
    foreach entry [lrange $entries 1 end] {
        set entry [string trim $entry]
        if {$entry eq ""} continue
        
        set lines [split $entry "\n"]
        if {[llength $lines] < 2} continue
        
        # Parse sequence name and data
        set seq_name_line ""
        set seq_data_line ""
        
        foreach line [lrange $lines 1 end] {
            set line [string trim $line]
            if {[string match ">*" $line]} {
                set seq_name_line [string range $line 1 end]
            } elseif {[string match "\{*\}*" $line]} {
                # Remove { at start and }* at end
                set seq_data_line [string range $line 1 end-2]
            }
        }
        
        if {$seq_name_line ne "" && $seq_data_line ne ""} {
            lappend sequences [list $seq_name_line $seq_data_line]
        }
    }
    
    return $sequences
}

# Extract motifs from MMA format alignment
# Returns nested list with motif data for each sequence
proc extractMotifsFromMMA {mma_file column_ranges} {
    set sequences [parseMMAFile $mma_file]
    
    if {[llength $sequences] == 0} {
        puts "Error: No sequences found in MMA file"
        return {}
    }
    
    # Process sequences and separate flanking regions from alignment
    set processed_sequences {}
    set reference_clean_alignment ""
    
    for {set i 0} {$i < [llength $sequences]} {incr i} {
        set seq_data [lindex $sequences $i]
        set name [lindex $seq_data 0]
        set full_seq [lindex $seq_data 1]
        
        # Extract flanking regions
        set flanking_start ""
        set flanking_end ""
        set alignment_seq $full_seq
        
        # Remove start flanking region in parentheses
        if {[regexp {^\(([^)]*)\)} $full_seq match flanking_start]} {
            set alignment_seq [string range $full_seq [string length $match] end]
        }
        
        # Remove end flanking region in parentheses
        if {[regexp {\(([^)]*)\)$} $alignment_seq match flanking_end]} {
            set end_pos [string first $match $alignment_seq]
            set alignment_seq [string range $alignment_seq 0 [expr $end_pos - 1]]
        }
        
        # Create clean alignment (uppercase + dashes only)
        set clean_alignment ""
        for {set j 0} {$j < [string length $alignment_seq]} {incr j} {
            set char [string index $alignment_seq $j]
            if {[string is upper $char] || $char eq "-"} {
                append clean_alignment $char
            }
        }
        
        lappend processed_sequences [list $name $flanking_start $flanking_end $alignment_seq $clean_alignment]
        
        # Use first sequence clean alignment as reference
        if {$i == 0} {
            set reference_clean_alignment $clean_alignment
        }
    }
    
    # Validate column ranges against first sequence clean alignment
    set ref_length [string length $reference_clean_alignment]
    foreach range $column_ranges {
        set start_col [lindex $range 0]
        set end_col [lindex $range 1]
        if {$start_col < 1 || $end_col > $ref_length} {
            puts "Error: Column range $start_col-$end_col exceeds reference clean alignment length ($ref_length)"
            return {}
        }
    }
    
    # Extract motifs for each sequence
    set all_results {}
    
    foreach seq_data $processed_sequences {
        set seq_name [lindex $seq_data 0]
        set flanking_start [lindex $seq_data 1]
        set flanking_end [lindex $seq_data 2]
        set full_alignment [lindex $seq_data 3]
        set clean_alignment [lindex $seq_data 4]
        
        set seq_motifs {}
        
        foreach range $column_ranges {
            set start_col [lindex $range 0]
            set end_col [lindex $range 1]
            
            # Handle sequences shorter than reference
            if {$end_col > [string length $clean_alignment]} {
                lappend seq_motifs [list "---" "N/A" "TRUNCATED"]
            } else {
                # Extract motif from clean alignment (1-based indexing)
                set motif_region_clean [string range $clean_alignment [expr $start_col - 1] [expr $end_col - 1]]
                set ungapped_motif [string map {"-" ""} $motif_region_clean]
                
                if {$ungapped_motif ne ""} {
                    # Calculate original residue positions in the full protein sequence
                    # Map clean alignment start position to full alignment position
                    set clean_pos 0
                    set full_pos 0
                    set start_full_pos 0
                    
                    while {$clean_pos < [expr $start_col - 1] && $full_pos < [string length $full_alignment]} {
                        set char [string index $full_alignment $full_pos]
                        if {[string is upper $char] || $char eq "-"} {
                            if {$clean_pos == [expr $start_col - 1]} {
                                set start_full_pos $full_pos
                                break
                            }
                            incr clean_pos
                        }
                        incr full_pos
                    }
                    
                    if {$clean_pos == [expr $start_col - 1]} {
                        set start_full_pos $full_pos
                    }
                    
                    # Count residues before the motif start in original protein
                    set residues_before_motif [string range $full_alignment 0 [expr $start_full_pos - 1]]
                    set residue_count_before [string length $flanking_start]
                    
                    for {set k 0} {$k < [string length $residues_before_motif]} {incr k} {
                        set char [string index $residues_before_motif $k]
                        if {[string is alpha $char]} {
                            incr residue_count_before
                        }
                    }
                    
                    set original_start [expr $residue_count_before + 1]
                    set original_end [expr $original_start + [string length $ungapped_motif] - 1]
                    
                    # Get the alignment region with insertions for display
                    set end_full_pos $start_full_pos
                    set clean_chars_found 0
                    while {$end_full_pos < [string length $full_alignment] && $clean_chars_found < [expr $end_col - $start_col + 1]} {
                        set char [string index $full_alignment $end_full_pos]
                        if {[string is upper $char] || $char eq "-"} {
                            incr clean_chars_found
                        }
                        incr end_full_pos
                    }
                    
                    set alignment_region_with_insertions [string range $full_alignment $start_full_pos [expr $end_full_pos - 1]]
                    
                    lappend seq_motifs [list $ungapped_motif "$original_start-$original_end" $alignment_region_with_insertions]
                } else {
                    # All gaps case
                    set clean_pos 0
                    set full_pos 0
                    set start_full_pos 0
                    
                    while {$clean_pos < [expr $start_col - 1] && $full_pos < [string length $full_alignment]} {
                        set char [string index $full_alignment $full_pos]
                        if {[string is upper $char] || $char eq "-"} {
                            if {$clean_pos == [expr $start_col - 1]} {
                                set start_full_pos $full_pos
                                break
                            }
                            incr clean_pos
                        }
                        incr full_pos
                    }
                    
                    set end_full_pos $start_full_pos
                    set clean_chars_found 0
                    while {$end_full_pos < [string length $full_alignment] && $clean_chars_found < [expr $end_col - $start_col + 1]} {
                        set char [string index $full_alignment $end_full_pos]
                        if {[string is upper $char] || $char eq "-"} {
                            incr clean_chars_found
                        }
                        incr end_full_pos
                    }
                    
                    set alignment_region_with_insertions [string range $full_alignment $start_full_pos [expr $end_full_pos - 1]]
                    
                    lappend seq_motifs [list "---" "N/A" $alignment_region_with_insertions]
                }
            }
        }
        
        lappend all_results [list $seq_name $seq_motifs]
    }
    
    return $all_results
}

# Get alignment residue IDs for each sequence
# Returns nested list mapping alignment positions to residue IDs
proc getAlignmentResidueIds {mma_file} {
    set sequences [parseMMAFile $mma_file]
    set result {}
    
    foreach seq $sequences {
        set seq_name [lindex $seq 0]
        set full_seq [lindex $seq 1]
        
        # Extract flanking regions
        set flanking_start ""
        set flanking_end ""
        set alignment_seq $full_seq
        
        # Remove start flanking region in parentheses
        if {[regexp {^\(([^)]*)\)} $full_seq match flanking_start]} {
            set alignment_seq [string range $full_seq [string length $match] end]
        }
        
        # Remove end flanking region in parentheses
        if {[regexp {\(([^)]*)\)$} $alignment_seq match flanking_end]} {
            set end_pos [string first $match $alignment_seq]
            set alignment_seq [string range $alignment_seq 0 [expr $end_pos - 1]]
        }
        
        # Calculate residue IDs for alignment positions only (uppercase + dashes)
        # Two-pass algorithm to handle insertions correctly
        
        set residue_array {}
        set protein_resnum [string length $flanking_start]  ;# Start counting from flanking region
        
        # Create a mapping of alignment positions to protein residue IDs
        set alignment_to_protein {}
        
        for {set i 0} {$i < [string length $alignment_seq]} {incr i} {
            set char [string index $alignment_seq $i]
            
            if {[string is upper $char]} {
                # Uppercase: aligned residue - increment protein ID and record mapping
                incr protein_resnum
                lappend alignment_to_protein $protein_resnum
            } elseif {$char eq "-"} {
                # Dash: gap - no protein residue at this alignment position
                lappend alignment_to_protein "NaN"
            } elseif {[string is lower $char]} {
                # Lowercase: insertion - increment protein ID but don't record for alignment
                incr protein_resnum
                # Don't add to alignment_to_protein - insertions are not alignment positions
            }
        }
        
        # The alignment_to_protein list now contains the correct mapping
        set residue_array $alignment_to_protein
        
        lappend result [list $seq_name $residue_array]
    }
    
    return $result
}

# Print motif extraction results in formatted output
proc printMotifResults {results column_ranges} {
    puts "MOTIF EXTRACTION FROM MMA (Column ranges: $column_ranges)"
    puts [string repeat "=" 120]
    
    # Debug: Check if we received the wrong type of results
    if {[llength $results] > 0} {
        set first_seq [lindex $results 0]
        set first_seq_data [lindex $first_seq 1]
        
        # Check if this looks like residue ID data (list of numbers/NaN) instead of motif data
        if {[llength $first_seq_data] > 50} {
            # This looks like residue ID data, not motif data
            puts "ERROR: printMotifResults received residue ID data instead of motif data!"
            puts "Expected: nested list with motif information"
            puts "Received: list with [llength $first_seq_data] elements (likely residue IDs)"
            puts "HINT: Use printMotifResults with motif_results, not residue_results"
            return
        }
    }
    
    # Print header
    set header "Sequence_Name"
    for {set i 0} {$i < [llength $column_ranges]} {incr i} {
        append header [format ", Motif_%d, Alignment_Region_%d, Original_Range_%d" [expr $i + 1] [expr $i + 1] [expr $i + 1]]
    }
    puts $header
    puts [string repeat "-" [string length $header]]
    
    # Print data rows
    foreach seq_result $results {
        set seq_name [lindex $seq_result 0]
        set motifs [lindex $seq_result 1]
        
        set line $seq_name
        foreach motif_data $motifs {
            set motif [lindex $motif_data 0]
            set range [lindex $motif_data 1]
            set alignment_region [lindex $motif_data 2]
            append line [format ", %s, %s, %s" $motif $alignment_region $range]
        }
        puts $line
    }
}

# Print residue ID mapping results
proc printResidueIds {results {max_positions 50}} {
    puts "\nRESIDUE ID MAPPING"
    puts [string repeat "=" 80]
    
    # Print header
    set header "Sequence"
    for {set i 1} {$i <= $max_positions} {incr i} {
        append header [format "\tPos_%d" $i]
    }
    puts $header
    puts [string repeat "-" 120]
    
    # Print data rows
    foreach seq_result $results {
        set seq_name [lindex $seq_result 0]
        set residues [lindex $seq_result 1]
        
        set line $seq_name
        set pos_count 0
        foreach residue $residues {
            if {$pos_count >= $max_positions} break
            append line "\t$residue"
            incr pos_count
        }
        puts $line
    }
}

# Helper function to get motifs for a specific sequence
proc getSequenceMotifs {results sequence_name} {
    foreach seq_result $results {
        if {[lindex $seq_result 0] eq $sequence_name} {
            return [lindex $seq_result 1]
        }
    }
    return {}
}

# Helper function to get residue IDs for a specific sequence
proc getSequenceResidueIds {results sequence_name} {
    foreach seq_result $results {
        if {[lindex $seq_result 0] eq $sequence_name} {
            return [lindex $seq_result 1]
        }
    }
    return {}
}

# Helper function to list all sequence names in results
proc getSequenceNames {results} {
    set names {}
    foreach seq_result $results {
        lappend names [lindex $seq_result 0]
    }
    return $names
}

# Main analysis function
proc analyzeMMAlignment {mma_file column_ranges} {
    puts "Analyzing MMA alignment file: $mma_file"
    puts "Column ranges: $column_ranges"
    puts ""
    
    # Extract motifs
    set motif_results [extractMotifsFromMMA $mma_file $column_ranges]
    if {[llength $motif_results] == 0} {
        puts "Error: Failed to extract motifs"
        return {}
    }
    
    # Get residue ID mapping
    set residue_results [getAlignmentResidueIds $mma_file]
    
    # Print results
    printMotifResults $motif_results $column_ranges
    puts ""
    printResidueIds $residue_results
    
    # Return both results as a nested list
    return [list $motif_results $residue_results]
}

# Helper function to get aligned residue pairs between two sequences (excluding NaN positions)
proc getAlignedResiduePairs {residue_results seq1_name seq2_name} {
    # Get residue arrays for both sequences
    set seq1_residues [getSequenceResidueIds $residue_results $seq1_name]
    set seq2_residues [getSequenceResidueIds $residue_results $seq2_name]
    
    if {[llength $seq1_residues] == 0} {
        puts "Error: Sequence '$seq1_name' not found in results"
        return {}
    }
    
    if {[llength $seq2_residues] == 0} {
        puts "Error: Sequence '$seq2_name' not found in results"
        return {}
    }
    
    # Find aligned residue pairs (both must be non-NaN)
    set aligned_pairs {}
    set min_length [expr {[llength $seq1_residues] < [llength $seq2_residues] ? [llength $seq1_residues] : [llength $seq2_residues]}]
    
    for {set i 0} {$i < $min_length} {incr i} {
        set res1 [lindex $seq1_residues $i]
        set res2 [lindex $seq2_residues $i]
        
        # Only include if both residues are not NaN
        if {$res1 ne "NaN" && $res2 ne "NaN"} {
            lappend aligned_pairs [list $res1 $res2 [expr $i + 1]]
        }
    }
    
    return $aligned_pairs
}

# Helper function to print aligned residue pairs in formatted output
proc printAlignedResiduePairs {aligned_pairs seq1_name seq2_name} {
    puts "\nALIGNED RESIDUE PAIRS (Non-gap positions only)"
    puts "Sequence 1: $seq1_name"
    puts "Sequence 2: $seq2_name"
    puts [string repeat "=" 60]
    puts [format "%-15s %-15s %-15s" "Seq1_ResId" "Seq2_ResId" "Alignment_Pos"]
    puts [string repeat "-" 60]
    
    foreach pair $aligned_pairs {
        set res1 [lindex $pair 0]
        set res2 [lindex $pair 1]
        set align_pos [lindex $pair 2]
        puts [format "%-15s %-15s %-15s" $res1 $res2 $align_pos]
    }
    
    puts "\nTotal aligned residue pairs: [llength $aligned_pairs]"
}

# Get aligned residue lists for VMD structural alignment
# This function processes all alignment positions row-wise and returns only those positions 
# where BOTH sequences have actual residues (non-NaN/non-gap). The returned lists will
# have equal lengths and can be used directly for VMD structural alignment.
proc getAlignedResidueListsForVMD {residue_results seq1_name seq2_name} {
    # Get residue arrays for both sequences (all 240 alignment positions)
    set seq1_residues [getSequenceResidueIds $residue_results $seq1_name]
    set seq2_residues [getSequenceResidueIds $residue_results $seq2_name]
    
    if {[llength $seq1_residues] == 0} {
        puts "Error: Sequence '$seq1_name' not found in results"
        return [list {} {}]
    }
    
    if {[llength $seq2_residues] == 0} {
        puts "Error: Sequence '$seq2_name' not found in results"
        return [list {} {}]
    }
    
    # Create lists for residues that exist in both sequences at the same alignment position
    set seq1_list {}
    set seq2_list {}
    
    # Go through all alignment positions
    set min_length [expr {[llength $seq1_residues] < [llength $seq2_residues] ? [llength $seq1_residues] : [llength $seq2_residues]}]
    
    for {set i 0} {$i < $min_length} {incr i} {
        set res1 [lindex $seq1_residues $i]
        set res2 [lindex $seq2_residues $i]
        
        # Only include if both sequences have actual residues (non-NaN) at this alignment position
        if {$res1 ne "NaN" && $res2 ne "NaN"} {
            lappend seq1_list $res1
            lappend seq2_list $res2
        }
    }
    
    puts "  Alignment row-wise filtering results:"
    puts "    Total alignment positions: $min_length"
    puts "    Positions where both sequences have residues: [llength $seq1_list]"
    puts "    Seq1 residues for alignment: [llength $seq1_list]"
    puts "    Seq2 residues for alignment: [llength $seq2_list]"
    
    return [list $seq1_list $seq2_list]
}

# Helper function to filter aligned residue lists further based on actual residues present in VMD structures
# This function takes the output from getAlignedResidueListsForVMD and checks if residues actually exist in structures
proc getAlignedResidueListsForVMDFiltered {residue_results seq1_name seq2_name molid1 molid2} {
    # Get the aligned residue lists (only positions where both sequences have non-NaN residues)
    set aligned_lists [getAlignedResidueListsForVMD $residue_results $seq1_name $seq2_name]
    set seq1_list [lindex $aligned_lists 0]
    set seq2_list [lindex $aligned_lists 1]
    
    # Check which of these aligned residues actually exist in the VMD structures
    set filtered_seq1_list {}
    set filtered_seq2_list {}
    
    for {set i 0} {$i < [llength $seq1_list]} {incr i} {
        set res1 [lindex $seq1_list $i]
        set res2 [lindex $seq2_list $i]
        
        # Check if residues exist in both structures
        set sel1 [atomselect $molid1 "resid $res1 and name CA"]
        set sel2 [atomselect $molid2 "resid $res2 and name CA"]
        
        set count1 [$sel1 num]
        set count2 [$sel2 num]
        
        $sel1 delete
        $sel2 delete
        
        # Only include if both residues exist in the actual VMD structures
        if {$count1 > 0 && $count2 > 0} {
            lappend filtered_seq1_list $res1
            lappend filtered_seq2_list $res2
        }
    }
    
    puts "  VMD structure filtering results:"
    puts "    Aligned residue pairs from sequences: [llength $seq1_list]"
    puts "    Pairs with residues present in both VMD structures: [llength $filtered_seq1_list]"
    puts "    Final seq1 residues for VMD alignment: [llength $filtered_seq1_list]"
    puts "    Final seq2 residues for VMD alignment: [llength $filtered_seq2_list]"
    
    return [list $filtered_seq1_list $filtered_seq2_list]
}

# Helper function to check residue coverage in VMD structures
proc checkResidueCoverage {residue_list molid sequence_name} {
    set present_count 0
    set missing_residues {}
    
    foreach resid $residue_list {
        set sel [atomselect $molid "resid $resid and name CA"]
        set count [$sel num]
        $sel delete
        
        if {$count > 0} {
            incr present_count
        } else {
            lappend missing_residues $resid
        }
    }
    
    puts "  Residue coverage for $sequence_name (molid $molid):"
    puts "    Total expected residues: [llength $residue_list]"
    puts "    Present residues: $present_count"
    puts "    Missing residues: [llength $missing_residues]"
    if {[llength $missing_residues] > 0 && [llength $missing_residues] <= 10} {
        puts "    Missing residue IDs: $missing_residues"
    } elseif {[llength $missing_residues] > 10} {
        puts "    Missing residue IDs (first 10): [lrange $missing_residues 0 9]..."
    }
    
    return [list $present_count [llength $missing_residues]]
}

# Debug function to understand alignment structure
proc debugAlignmentStructure {mma_file {max_seqs 3}} {
    set sequences [parseMMAFile $mma_file]
    puts "DEBUG: Alignment Structure Analysis"
    puts [string repeat "=" 80]
    
    set seq_count 0
    foreach seq $sequences {
        if {$seq_count >= $max_seqs} break
        
        set seq_name [lindex $seq 0]
        set full_seq [lindex $seq 1]
        
        puts "\nSequence: $seq_name"
        puts "Full sequence: $full_seq"
        
        # Extract components
        set flanking_start ""
        set flanking_end ""
        set alignment_seq $full_seq
        
        # Remove start flanking region
        if {[regexp {^\(([^)]*)\)} $full_seq match flanking_start]} {
            set alignment_seq [string range $full_seq [string length $match] end]
        }
        
        # Remove end flanking region
        if {[regexp {\(([^)]*)\)$} $alignment_seq match flanking_end]} {
            set end_pos [string first $match $alignment_seq]
            set alignment_seq [string range $alignment_seq 0 [expr $end_pos - 1]]
        }
        
        puts "Flanking start: ($flanking_start)"
        puts "Flanking end: ($flanking_end)"
        puts "Alignment sequence: $alignment_seq"
        
        # Create clean alignment (only uppercase and dashes)
        set clean_alignment ""
        set uppercase_count 0
        set dash_count 0
        set lowercase_count 0
        
        for {set i 0} {$i < [string length $alignment_seq]} {incr i} {
            set char [string index $alignment_seq $i]
            if {[string is upper $char]} {
                append clean_alignment $char
                incr uppercase_count
            } elseif {$char eq "-"} {
                append clean_alignment $char
                incr dash_count
            } elseif {[string is lower $char]} {
                incr lowercase_count
            }
        }
        
        puts "Clean alignment (uppercase + dashes): $clean_alignment"
        puts "Stats - Uppercase: $uppercase_count, Dashes: $dash_count, Insertions: $lowercase_count"
        puts "Clean alignment length: [string length $clean_alignment]"
        
        incr seq_count
    }
}

# Debug function to check alignment lengths and identify discrepancies
proc debugAlignmentLengths {residue_results seq1_name seq2_name} {
    set seq1_residues [getSequenceResidueIds $residue_results $seq1_name]
    set seq2_residues [getSequenceResidueIds $residue_results $seq2_name]
    
    puts "DEBUG: Alignment Length Analysis"
    puts "================================="
    puts "Sequence 1: $seq1_name"
    puts "Sequence 2: $seq2_name"
    puts "Seq1 total positions: [llength $seq1_residues]"
    puts "Seq2 total positions: [llength $seq2_residues]"
    
    # Count non-NaN positions for each sequence
    set seq1_non_nan 0
    set seq2_non_nan 0
    set common_non_nan 0
    
    set min_length [expr {[llength $seq1_residues] < [llength $seq2_residues] ? [llength $seq1_residues] : [llength $seq2_residues]}]
    
    for {set i 0} {$i < $min_length} {incr i} {
        set res1 [lindex $seq1_residues $i]
        set res2 [lindex $seq2_residues $i]
        
        if {$res1 ne "NaN"} {
            incr seq1_non_nan
        }
        if {$res2 ne "NaN"} {
            incr seq2_non_nan
        }
        if {$res1 ne "NaN" && $res2 ne "NaN"} {
            incr common_non_nan
        }
    }
    
    puts "Seq1 non-NaN positions: $seq1_non_nan"
    puts "Seq2 non-NaN positions: $seq2_non_nan" 
    puts "Common non-NaN positions: $common_non_nan"
    
    # Show first 20 positions to see the pattern
    puts "\nFirst 20 alignment positions:"
    puts "Pos\tSeq1\tSeq2\tBoth_Valid"
    for {set i 0} {$i < 20 && $i < $min_length} {incr i} {
        set res1 [lindex $seq1_residues $i]
        set res2 [lindex $seq2_residues $i]
        set both_valid [expr {$res1 ne "NaN" && $res2 ne "NaN" ? "YES" : "NO"}]
        puts "[expr $i+1]\t$res1\t$res2\t$both_valid"
    }
    
    return [list $seq1_non_nan $seq2_non_nan $common_non_nan]
}

# Example usage and configuration
# Uncomment the following lines to run the analysis

# CONFIGURATION
# set path "/mnt/mydrive/pseudokinase/James_murphy/seq_analysis/"
# set mma_file "${path}pseudokinases_kinase_domains_msa.mma"
# set column_ranges {{27 30} {117 119} {137 139}}

# if {[file exists $mma_file]} {
#     puts "Processing MMA format file..."
#     set results [analyzeMMAlignment $mma_file $column_ranges]
#     
#     # Access specific results
#     set motif_results [lindex $results 0]
#     set residue_results [lindex $results 1]
#     
#     # Example: Get motifs for a specific sequence
#     set seq_names [getSequenceNames $motif_results]
#     if {[llength $seq_names] > 0} {
#         set first_seq [lindex $seq_names 0]
#         set motifs [getSequenceMotifs $motif_results $first_seq]
#         puts "\nMotifs for $first_seq: $motifs"
#     }
# } else {
#     puts "MMA file not found: $mma_file"
# }




# Check if running in VMD environment
if {[info exists ::env(VMDDIR)]} {
    puts "Info) Motif Reader loaded in VMD environment"
} else {
    puts "Motif Reader loaded successfully!"
}

puts "Available motif analysis functions:"
puts "  Core Functions:"
puts "    - parseMMAFile {mma_file}"
puts "    - extractMotifsFromMMA {mma_file column_ranges}"
puts "    - getAlignmentResidueIds {mma_file}"
puts "    - analyzeMMAlignment {mma_file column_ranges}"
puts "  Helper Functions:"
puts "    - getSequenceMotifs {results sequence_name}"
puts "    - getSequenceResidueIds {results sequence_name}"
puts "    - getSequenceNames {results}"
puts "  Alignment Functions:"
puts "    - getAlignedResiduePairs {residue_results seq1_name seq2_name}"
puts "    - printAlignedResiduePairs {aligned_pairs seq1_name seq2_name}"
puts "    - getAlignedResidueListsForVMD {residue_results seq1_name seq2_name}"
puts "    - getAlignedResidueListsForVMDFiltered {residue_results seq1_name seq2_name molid1 molid2}"
puts "    - checkResidueCoverage {residue_list molid sequence_name}"
puts "  Debug Functions:"
puts "    - debugAlignmentStructure {mma_file {max_seqs 3}}"
puts "    - debugAlignmentLengths {residue_results seq1_name seq2_name}"

proc listMotifReaderProcs {} {
    # List all available procedures in motifReader.tcl with descriptions
    # Use this in VMD console to see what functions are available
    
    puts "\n=== Motif Reader Procedures ==="
    puts "Available procedures for MMA alignment analysis:\n"
    
    puts "CORE PARSING FUNCTIONS:"
    puts "  parseMMAFile mma_file"
    puts "    -> Parse MMA format file and return sequence data as nested lists"
    puts "  extractMotifsFromMMA mma_file column_ranges"
    puts "    -> Extract motifs from MMA alignment for specified column ranges"
    puts "  getAlignmentResidueIds mma_file"
    puts "    -> Get alignment residue ID mapping for each sequence"
    puts "  analyzeMMAlignment mma_file column_ranges"
    puts "    -> Main analysis function - extract motifs and residue IDs"
    puts ""
    
    puts "SEQUENCE HELPER FUNCTIONS:"
    puts "  getSequenceMotifs results sequence_name"
    puts "    -> Get motifs for a specific sequence from results"
    puts "  getSequenceResidueIds results sequence_name"
    puts "    -> Get residue IDs for a specific sequence from results"
    puts "  getSequenceNames results"
    puts "    -> List all sequence names in results"
    puts ""
    
    puts "ALIGNMENT FUNCTIONS:"
    puts "  getAlignedResiduePairs residue_results seq1_name seq2_name"
    puts "    -> Find aligned residue pairs between two sequences (excluding gaps)"
    puts "  printAlignedResiduePairs aligned_pairs seq1_name seq2_name"
    puts "    -> Print aligned residue pairs in formatted output"
    puts "  getAlignedResidueListsForVMD residue_results seq1_name seq2_name"
    puts "    -> Get aligned residue lists for VMD structural alignment"
    puts "  getAlignedResidueListsForVMDFiltered residue_results seq1_name seq2_name molid1 molid2"
    puts "    -> Get filtered aligned residue lists checking VMD structure availability"
    puts ""
    
    puts "STRUCTURE VALIDATION FUNCTIONS:"
    puts "  checkResidueCoverage residue_list molid sequence_name"
    puts "    -> Check residue coverage in VMD structures"
    puts ""
    
    puts "OUTPUT FUNCTIONS:"
    puts "  printMotifResults results column_ranges"
    puts "    -> Print motif extraction results in formatted output"
    puts "  printResidueIds results ?max_positions?"
    puts "    -> Print residue ID mapping results (default max_positions=50)"
    puts ""
    
    puts "DEBUG FUNCTIONS:"
    puts "  debugAlignmentStructure mma_file ?max_seqs?"
    puts "    -> Debug function to understand alignment structure (default max_seqs=3)"
    puts "  debugAlignmentLengths residue_results seq1_name seq2_name"
    puts "    -> Debug function to check alignment lengths and discrepancies"
    puts "  identifyResultType results"
    puts "    -> Identify if results are motif results or residue ID results"
    puts ""
    
    puts "EXAMPLE FUNCTIONS:"
    puts "  examplePrintMotifResults"
    puts "    -> Show detailed example of how to use printMotifResults correctly"
    puts "  quickExampleMotifExtraction"
    puts "    -> Show quick example of complete motif extraction workflow"
    puts ""
    
    puts "HELP:"
    puts "  listMotifReaderProcs"
    puts "    -> Display this help message"
    puts "  listMotifReaderProcsCompact"
    puts "    -> Display compact list of procedures"
    puts ""
    
    puts "EXAMPLE USAGE:"
    puts "  set mma_file \"alignment.mma\""
    puts "  set column_ranges {{27 30} {117 119} {137 139}}"
    puts "  set results [analyzeMMAlignment \$mma_file \$column_ranges]"
    puts "  set motif_results [lindex \$results 0]"
    puts "  set residue_results [lindex \$results 1]"
    puts "  set aligned_lists [getAlignedResidueListsForVMD \$residue_results \"seq1\" \"seq2\"]"
    puts ""
    
    puts "=== End of Motif Reader Procedures ==="
}

proc helpMotifReader {} {
    # Alias for listMotifReaderProcs - shorter command
    listMotifReaderProcs
}

proc motifReaderHelp {} {
    # Another alias for listMotifReaderProcs
    listMotifReaderProcs
}

proc listMotifReaderProcsCompact {} {
    # Show a compact list of all available procedures
    
    puts "\n=== Motif Reader Procedures (Compact) ==="
    puts "Core:     parseMMAFile, extractMotifsFromMMA, getAlignmentResidueIds"
    puts "          analyzeMMAlignment"
    puts "Helper:   getSequenceMotifs, getSequenceResidueIds, getSequenceNames"
    puts "Align:    getAlignedResiduePairs, printAlignedResiduePairs"
    puts "          getAlignedResidueListsForVMD, getAlignedResidueListsForVMDFiltered"
    puts "Valid:    checkResidueCoverage"
    puts "Output:   printMotifResults, printResidueIds"
    puts "Debug:    debugAlignmentStructure, debugAlignmentLengths, identifyResultType"
    puts "Examples: examplePrintMotifResults, quickExampleMotifExtraction"
puts "Help:     listMotifReaderProcs, helpMotifReader, listMotifReaderProcsCompact"
    puts "=========================================="
}

proc identifyResultType {results} {
    # Helper function to identify if results are motif results or residue ID results
    # This helps users understand which print function to use
    
    if {[llength $results] == 0} {
        puts "Empty results"
        return "empty"
    }
    
    set first_seq [lindex $results 0]
    set seq_name [lindex $first_seq 0]
    set seq_data [lindex $first_seq 1]
    
    puts "Result type analysis for first sequence: $seq_name"
    puts "Data structure: [llength $seq_data] elements"
    
    if {[llength $seq_data] > 50} {
        puts "This appears to be RESIDUE ID results"
        puts "Use: printResidueIds \$results"
        puts "First few elements: [lrange $seq_data 0 4]"
        return "residue_ids"
    } else {
        puts "This appears to be MOTIF results"
        puts "Use: printMotifResults \$results \$column_ranges"
        puts "Number of motifs: [llength $seq_data]"
        if {[llength $seq_data] > 0} {
            set first_motif [lindex $seq_data 0]
            puts "First motif structure: $first_motif"
        }
        return "motifs"
    }
}

# EXAMPLE USAGE FOR printMotifResults
# ==================================

# Example: How to use printMotifResults correctly
proc examplePrintMotifResults {} {
    puts "\n=== Example Usage of printMotifResults ==="
    puts "This example shows the correct workflow for using printMotifResults\n"
    
    # Step 1: Define your MMA file and column ranges
    set mma_file "/path/to/your/alignment.mma"
    set column_ranges {{27 30} {117 119} {137 139}}
    
    puts "Step 1: Set up file and column ranges"
    puts "  MMA file: $mma_file"
    puts "  Column ranges: $column_ranges"
    puts "  (These represent motif regions you want to extract)\n"

    
    # Step 2: Extract motifs from MMA file
    puts "Step 2: Extract motifs from MMA file"
    set motif_results [extractMotifsFromMMA $mma_file $column_ranges]
    
    if {[llength $motif_results] == 0} {
        puts "Error: No motifs extracted"
        return
    }
    
    puts "  Successfully extracted motifs for [llength $motif_results] sequences\n"
    
    # Step 3: Print the motif results
    puts "Step 3: Print motif results using printMotifResults"
    printMotifResults $motif_results $column_ranges
    
    puts "\nStep 4: Alternative - print motifs for specific sequence"
    set seq_names [getSequenceNames $motif_results]
    if {[llength $seq_names] > 0} {
        set first_seq [lindex $seq_names 0]
        set motifs [getSequenceMotifs $motif_results $first_seq]
        puts "Motifs for $first_seq:"
        foreach motif $motifs {
            puts "  Motif: [lindex $motif 0], Range: [lindex $motif 1], Alignment: [lindex $motif 2]"
        }
    }
}

