# Motif Reader Tcl Port - Complete Implementation

## Overview
The `motif_reader.tcl` file contains a complete Tcl port of the Python `motif_reader.py` functionality. All functions return results as arrays/lists instead of saving to files, making it suitable for integration with VMD and other Tcl-based workflows.

## Key Features
- **Complete Python Functionality Port**: All major functions from the Python version
- **Array-based Results**: No file I/O for results - everything returned as Tcl lists
- **MMA Format Support**: Full parsing of MMA (Multiple Alignment) format files
- **Motif Extraction**: Extract specific sequence motifs based on column ranges
- **Residue Mapping**: Map alignment positions to original protein residue numbers
- **Helper Functions**: Utilities for accessing specific sequence data

## Available Functions

### Core Functions

#### `parseMMAFile {mma_file}`
- **Purpose**: Parse MMA format file and return sequence data
- **Returns**: List of sequences, each containing [name, sequence_data]
- **Example**: 
  ```tcl
  set sequences [parseMMAFile "alignment.mma"]
  # Returns: {{seq1 "ACGT..."} {seq2 "ATCG..."} ...}
  ```

#### `extractMotifsFromMMA {mma_file column_ranges}`
- **Purpose**: Extract motifs from specific column ranges in the alignment
- **Parameters**:
  - `mma_file`: Path to MMA format file
  - `column_ranges`: List of {start end} pairs (1-based indexing)
- **Returns**: Nested list with motif data for each sequence
- **Example**:
  ```tcl
  set motifs [extractMotifsFromMMA "file.mma" {{27 30} {117 119}}]
  # Returns: {{seq1 {{motif1 "27-30" "align_region1"} {motif2 "117-119" "align_region2"}}} ...}
  ```

#### `getAlignmentResidueNumbers {mma_file}`
- **Purpose**: Map alignment positions to original protein residue numbers
- **Returns**: List mapping each alignment position to residue number
- **Example**:
  ```tcl
  set residue_map [getAlignmentResidueNumbers "file.mma"]
  # Returns: {{seq1 {1 2 3 NaN 4 5...}} {seq2 {1 NaN 2 3 4...}} ...}
  ```

#### `analyzeMMAlignment {mma_file column_ranges}`
- **Purpose**: Complete analysis with formatted output
- **Returns**: Combined results [motif_results residue_results]
- **Example**:
  ```tcl
  set results [analyzeMMAlignment "file.mma" {{27 30} {117 119}}]
  set motif_results [lindex $results 0]
  set residue_results [lindex $results 1]
  ```

### Helper Functions

#### `getSequenceMotifs {results sequence_name}`
- **Purpose**: Get motifs for a specific sequence from results
- **Returns**: List of motifs for the specified sequence

#### `getSequenceResidueNumbers {results sequence_name}`
- **Purpose**: Get residue number mapping for a specific sequence
- **Returns**: List of residue numbers for the specified sequence

#### `getSequenceNames {results}`
- **Purpose**: Extract all sequence names from results
- **Returns**: List of sequence names

#### `printMotifResults {results column_ranges}`
- **Purpose**: Print formatted motif extraction results
- **Output**: Table format showing motifs, alignment regions, and ranges

#### `printResidueNumbers {results {max_positions 50}}`
- **Purpose**: Print residue number mapping in table format
- **Output**: Table showing alignment positions mapped to residue numbers

## Data Structure Details

### Motif Results Format
```tcl
{
  {sequence_name {
    {motif_sequence "start-end" "alignment_region_with_insertions"}
    {motif_sequence "start-end" "alignment_region_with_insertions"}
    ...
  }}
  ...
}
```

### Residue Number Results Format
```tcl
{
  {sequence_name {residue_num1 residue_num2 NaN residue_num3 ...}}
  ...
}
```

## Usage Examples

### Basic Motif Extraction
```tcl
# Source the script
source motif_reader.tcl

# Define parameters
set mma_file "pseudokinases_alignment.mma"
set column_ranges {{27 30} {117 119} {137 139}}

# Extract motifs
set motif_results [extractMotifsFromMMA $mma_file $column_ranges]

# Access results for specific sequence
set sequence_names [getSequenceNames $motif_results]
set first_seq [lindex $sequence_names 0]
set motifs [getSequenceMotifs $motif_results $first_seq]

puts "Motifs for $first_seq:"
foreach motif $motifs {
    set seq [lindex $motif 0]
    set range [lindex $motif 1]
    set align_region [lindex $motif 2]
    puts "  $seq ($range): $align_region"
}
```

### Complete Analysis with Output
```tcl
# Run complete analysis with formatted output
set results [analyzeMMAlignment $mma_file $column_ranges]

# Extract components
set motif_results [lindex $results 0]
set residue_results [lindex $results 1]

# Use results for further analysis
foreach seq_result $motif_results {
    set seq_name [lindex $seq_result 0]
    set motifs [lindex $seq_result 1]
    
    # Process each motif for this sequence
    foreach motif_data $motifs {
        set motif_seq [lindex $motif_data 0]
        set orig_range [lindex $motif_data 1]
        # ... custom processing ...
    }
}
```

### Integration with VMD
```tcl
# Example: Use motif results to select residues in VMD
source motif_reader.tcl

set motif_results [extractMotifsFromMMA $mma_file $column_ranges]
set target_sequence "your_sequence_name"
set motifs [getSequenceMotifs $motif_results $target_sequence]

# Process each motif for VMD selection
foreach motif_data $motifs {
    set range_str [lindex $motif_data 1]
    if {$range_str ne "N/A"} {
        # Parse range (e.g., "27-30" -> start=27, end=30)
        set range_parts [split $range_str "-"]
        set start [lindex $range_parts 0]
        set end [lindex $range_parts 1]
        
        # Create VMD selection
        set sel [atomselect top "resid $start to $end"]
        # ... use selection ...
        $sel delete
    }
}
```

## Testing
The functionality has been tested with the provided test script (`test_motif_tcl.tcl`) which demonstrates:
- Parsing of MMA files
- Motif extraction for multiple column ranges
- Residue number mapping
- Helper function usage
- Array-based result handling

## Comparison with Python Version
- **Full Feature Parity**: All Python functionality replicated
- **Performance**: Similar performance for typical alignment sizes
- **Integration**: Better suited for VMD/Tcl workflows
- **Output Format**: Arrays/lists instead of CSV files
- **Dependencies**: No external dependencies (unlike Python Bio module)

## Files
- **Main Script**: `/mnt/mydrive/pseudokinase/motif_reader.tcl`
- **Test Script**: `/mnt/mydrive/pseudokinase/test_motif_tcl.tcl`
- **Original Python**: `/mnt/mydrive/pseudokinase/motif_reader.py`

The Tcl port is complete and fully functional, providing all the capabilities of the original Python script while being optimized for array-based data handling suitable for VMD and other Tcl-based analysis workflows.
