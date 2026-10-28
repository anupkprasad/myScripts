### VMD entry point for myScripts custom Tcl
# Add to ~/.vmdrc (or VMD startup):
#   source /path/to/myScripts/callScripts.tcl
# or, if the repo is at ~/myScripts:
#   source ~/myScripts/callScripts.tcl
#
# Repo root is the directory of this file — works for any username / clone path.
# Modules under src/md/**/*.tcl are auto-loaded (see src/md/vmd_autoload.tcl).

set pathVmdScript [file normalize [file dirname [info script]]]
if {![file isdirectory [file join $pathVmdScript src md]]} {
    # Fallback if this file was copied elsewhere: ~/myScripts
    set pathVmdScript [file normalize [file join $::env(HOME) myScripts]]
}
source [file join $pathVmdScript src md vmd_autoload.tcl]
puts "Custom myScripts functions are integrated"
