#!/usr/bin/tclsh
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy
# get image of each frames for movie formations
# get the frames in the movie
proc snapshot {{dim 2000} {filename ""}} {
	display depthcue off
	if {$filename eq ""} {
		set frm [molinfo top get frame]
		set f_name [lindex [split [molinfo top get name] "."] 0]
		set filename ${f_name}_${frm}
	}
	puts $filename
	#render POV3 $filename  povray +W$dim +H$dim -I%s -O%s.tga +X +A +FT +UA +WT8
	render POV3 $filename \
    povray +W$dim +H$dim -I%s -O%s.tga \
    +A0.3 +R2 +UA +FT +WT12
	exec convert $filename.tga $filename.png
	exec rm -f $filename $filename.tga
	puts "rendering finished"
}


# take a snapshot of the trajectory at given frames
proc trajSnapshot {start step stop} {
    for {set i $start} {$i <= $stop} {incr i $step} {
        puts "Processing frame: $i"
        animate goto $i
        display update ui
        snapshot 1000
    }
}


proc snap360degree {{angleIncrement 3}} {
	set frm [molinfo top get frame]
	set f_name [lindex [split [molinfo top get name] "."] 0]

	for {set i 0} {$i < 360} {incr i $angleIncrement} {
		set angle $i
		set filename ${f_name}_${frm}_${angle}
		puts "Taking snapshot at angle $angle"
		snapshot 1000 $filename
		rotate y by $angleIncrement
	}
}
