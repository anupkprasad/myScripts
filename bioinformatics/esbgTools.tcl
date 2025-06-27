# esbgTools.tcl - A collection of scripts for ESBG tools

proc getFilelistFromSSCR {{input_file "provide_input.sscr" } {output_file "filelist_scorebased.txt"}} {
    # Open and read the file, skipping header line
    set fp [open $input_file r]
    set data {}

    # Skip first header line
    gets $fp line

    set count 0
    while {[gets $fp line] >= 0 && $count < 501} {
        set fields [regexp -inline -all {\S+} $line]
        lappend data $fields
        incr count
    }
    close $fp

    # Sort rows by score (second column, index 1)
    set sorted_data [lsort -real -index 1 $data]

    # Extract only filenames (filename from path)
    set filelist {}
    foreach row $sorted_data {
        set filepath [lindex $row 2]
        set filename [file tail $filepath]
        lappend filelist $filename
    }

    # Write the output to file
    set out [open $output_file w]
    foreach f $filelist {
        puts $out $f
    }
    close $out

    # Return the final list
    return $filelist
}