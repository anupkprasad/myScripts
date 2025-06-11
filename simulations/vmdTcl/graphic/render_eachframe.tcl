#!/usr/bin/tclsh
# Anup K. Prasad
# anupkprasad121@gmail.com
# Ph.D @IITB-Monash Research Academy
# get image of each frames for movie formations

set num [molinfo top get numframes]
# loop through the frames
for {set i 0} {$i < 150} {incr i} {

        # go to the given frame 
        animate goto $i
        # take the picture
        set filename snap.[format "%04d" $i].dat
        #pbc box_draw -style dashed -width 0.5
                               render Tachyon $filename "/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -format TARGA -res 1500 1500 -o %s.tga
                               #draw delete all 
    set i [expr {$i+4}]
}



set num [molinfo top get numframes]
# loop through the frames
for {set i 0} {$i < 150} {incr i} {

        # go to the given frame 
        animate goto $i
        # take the picture
        set filename snap.[format "%04d" $i].dat
        #pbc box_draw -style dashed -width 0.5
                               render Tachyon $filename "/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -format TARGA -res 1500 1500 -o %s.tga
                               #draw delete all 
    set i [expr {$i+4}]
}
