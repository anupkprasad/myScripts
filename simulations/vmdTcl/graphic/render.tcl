proc snapshot {{molid top} {dim 2000}} {
	display depthcue off
	set frm [molinfo $molid get frame]
	set f_name [lindex [split [molinfo $molid get name] "."] 0]
	set filename ${f_name}_${frm}.pov
	puts $filename
	render POV3 $filename  povray +W$dim +H$dim -I%s -O%s.tga +X +A +FT +UA +WT8
	exec convert $filename.tga $filename.png
	puts "redering finished"
}