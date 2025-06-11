# get the number of frames in the movie

        set num [molinfo top get numframes]
        # loop through the frames
        for {set i 0} {$i < $num} {incr i} {
        
                # go to the given frame 
                animate goto $i
                # take the picture
                set filename snap.[format "%04d" $i].dat
                pbc box_draw -style dashed -width 0.5
                                       render Tachyon $filename "/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -format TARGA -res 2500 1500 -o %s.tga
                                       draw delete all 
            set i [expr {$i+4}]
        }
