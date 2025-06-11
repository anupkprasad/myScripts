

#draw delete all
proc label_atom {selection_string label_string pad} {
    set sel [atomselect top $selection_string]
    if {[$sel num] != 1} {
        error "label_atom: '$selection_string' must select 1 atom"
    }
    # get the coordinates of the atom
    lassign [$sel get {x y z}] coord
    # and draw the text
    draw color red
    draw text $coord $label_string size 3 thickness 5
}