#!/usr/bin/tclsh
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy

proc get_repre { {molecule_number 0}} {
    mol modselect 0 $molecule_number all ; # modselect rep_number molecule_number select_method
    mol modstyle 0 $molecule_number Newcartoon
    mol modcolor 0 $molecule_number Chain
    #color Chain A cyan
    #color Chain B lime

    mol addrep $molecule_number
    mol modselect 1 $molecule_number not protein
    mol modstyle 1 $molecule_number licorice

    #You can do something like this to get both the repid and the repname.
    #Note that the repname is guaranteed to be unique for a given molecule,
    #but the repid may change if reps are deleted.
    mol addrep $molid
    set repid [expr [molinfo $molid get numreps] - 1]
    # repname is guaranteed to be unique within a molecule
    set repname [mol repname $molid $repid]
    #Later on, if you want to, e.g., delete a rep given the repname, you
    #can do something like:
    set repid [mol repindex $molid $repname]
    mol delrep $repid $molid
}