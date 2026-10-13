### VMD entry point for myScripts custom Tcl
# Add to ~/.vmdrc (or VMD startup):
#   source ~/myScripts/callScripts.tcl
#
# Modules under src/md/**/*.tcl are auto-loaded (see src/md/vmd_autoload.tcl).
# Drop a new analysis .tcl into src/md/analysis/ — no edit needed here.

set pathVmdScript [file normalize ~/myScripts]
source [file join $pathVmdScript src md vmd_autoload.tcl]
puts "Custom myScripts functions are integrated"
