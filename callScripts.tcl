set pathVmdScript "/home/anup/myScripts/simulations/vmdTcl"
set pathAnalysis "$pathVmdScript/mdAnalysis"
set pathGraphic "$pathVmdScript/graphic"
set pathSetup "$pathVmdScript/mdSetup"
set pathCentroid "$pathVmdScript/centroid"
set refStructures "/home/anup/scripts/reference"




source "$pathVmdScript/basicFunctions.tcl"
source "$pathVmdScript/trajFunctions.tcl"
source "$pathVmdScript/parseFiles.tcl"
source "$pathVmdScript/graphic/playAllMol.tcl"

source "$pathVmdScript/mdSetup/systemSetup.tcl"
source "$pathVmdScript/mdSetup/phosPTM.tcl"


source "$pathVmdScript/centroid/propgate2.tcl"
source "$pathVmdScript/mdAnalysis/ssTraj.tcl"
source "$pathVmdScript/mdAnalysis/gyra.tcl"






puts "Custom functions are integrated, thank you"