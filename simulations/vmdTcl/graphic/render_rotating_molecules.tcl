###For transparent image: Depth cueing should be off 

 set frame 0
        for {set i 0} {$i < 360} {incr i 3} {
                set filename snap.[format "%04d" $frame].pov
                render POV3 $filename  povray +W2000 +H2000 -I%s -O%s.tga +X +A +FT +UA
                incr frame
                rotate y by -3
        }


 #set frame 0
        #for {set i 0} {$i < 360} {incr i 5} {
                #set filename snap.[format "%04d" $frame].dat
                #render Tachyon $filename "/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -format TARGA -res 3000 3000 -o %s.tga
                #incr frame
                #rotate y by -5
        #}
