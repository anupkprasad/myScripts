# Universidad de Talca
# Centro de Bioinformatica y Simulacion Molecular
# Gonzalo Riadi, Ph.D. student.

# Script that counts the average number of (oxigens actually) of waters in a casquet
# of radius min to max from frame_start to frame_end using atom_ref as reference atom

# Variables to change
set min ### ; # lower limit of the hydration shell
set max ### ; # higher limit of the hydration shell
set atom_ref ###; # VMD index of the reference atom

# Don't need to change from here on
set frame_start 0;
set num_steps [molinfo top get numframes];
set frame_end $num_steps;
set stride 1;
set suma 0;

for {set frame $frame_start} {$frame < $frame_end} {} {
	set oxigens [atomselect top "name OH2 and (within $max of index $atom_ref) and (not within $min of index $atom_ref)" frame $frame];
	set num_oxigenos [$oxigens num];
	set suma [expr $suma + $num_oxigenos];
	set frame [expr $frame + $stride];
#	puts "frame $frame\t$suma\t$num_oxigenos\n";
}
set num_steps2 [expr (($frame_end - $frame_start) / $stride) + 1];
set inv_num_steps [expr 1.0 / $num_steps2];
set resul [vecscale $inv_num_steps $suma];
puts "The number of waters is $resul\n";
display update;

