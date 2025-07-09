set pathVmdScript [file normalize ~/myScripts]
set pathAnalysis "$pathVmdScript/simulations/vmdTcl/mdAnalysis"
set pathGraphic "$pathVmdScript/simulations/vmdTcl/graphic"
set pathSetup "$pathVmdScript/simulations/vmdTcl/mdSetup"
set pathCentroid "$pathVmdScript/simulations/vmdTcl/centroid"
set refStructures "$pathVmdScript/referencefiles/KinaseStructures"
set pathBioinfo "$pathVmdScript/bioinformatics"



source "$pathVmdScript/simulations/vmdTcl/basicFunctions.tcl"
source "$pathVmdScript/simulations/vmdTcl/trajFunctions.tcl"
source "$pathVmdScript/simulations/vmdTcl/parseFiles.tcl"
source "$pathVmdScript/simulations/vmdTcl/graphic/playAllMol.tcl"

source "$pathVmdScript/simulations/vmdTcl/mdSetup/systemSetup.tcl"
source "$pathVmdScript/simulations/vmdTcl/mdSetup/phosPTM.tcl"


source "$pathVmdScript/simulations/vmdTcl/centroid/propgate2.tcl"
source "$pathVmdScript/simulations/vmdTcl/mdAnalysis/ssTraj.tcl"
source "$pathVmdScript/simulations/vmdTcl/mdAnalysis/gyra.tcl"


############ for esbg tools  ########
source "$pathBioinfo/esbgTools.tcl"



puts "Custom functions are integrated, thank you"