

proc snapshot {{filename } {dim 2000}} {
	display depthcue off
	set frm [molinfo top get frame]
	if {$filename eq ""} {
		set f_name [lindex [split [molinfo top get name] "."] 0]
		set filename ${f_name}_${frm}.pov
	}

	puts $filename
	render POV3 $filename  povray +W$dim +H$dim -I%s -O%s.tga +X +A +FT +UA +WT8
	exec convert $filename.tga $filename.png
	puts "rendering finished: $filename.png"
}