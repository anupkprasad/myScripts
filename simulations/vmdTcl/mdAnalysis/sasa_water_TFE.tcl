#Anup_kumar (Ph.D. Student, IITB-Monash research academy)
#Last update: March 30 2023
#sel : selection of all molecules inbetween you dont want surface area  
#sel1/sel2: restrict: caculate sasa particular to this of sel

set sel1 [atomselect top "protein or water or ions"]
set sel2 [atomselect top "protein or resname TFE"]
set sel [atomselect top "protein"]
set n [molinfo top get numframes]
set output [open sasa_water_TFE.dat w]
puts $output "Data information: sasa of system $sel1 and $sel2 in which caculated to restrict1=$sel"
puts $output "frame system1=$sel1 system2=$sel2"
for {set i 0} {$i < $n} {incr i} {
        molinfo top set frame $i
        set sasa1 [measure sasa 1.4 $sel1 -restrict $sel]
        set sasa2 [measure sasa 1.4 $sel2 -restrict $sel]
        puts "\t \t progress: $i/$n"
        puts $output "$i $sasa1 $sasa2"
}
puts "\t \t progress: $n/$n"
puts "Done"
close $output
