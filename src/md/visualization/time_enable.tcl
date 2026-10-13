################------------------------#########
proc enabletrace {} {
  global vmd_frame;
  trace variable vmd_frame([molinfo top]) w drawcounter

}

proc disabletrace {} {
  global vmd_frame;
  trace vdelete vmd_frame([molinfo top]) w drawcounter
  
}
proc drawcounter { name element op } {
  global vmd_frame;
  draw delete all
# puts "callback!"
 draw color red
set nsperframe 1
set nsoffset 0
set time [format "%d ns" [expr ($vmd_frame([molinfo top]) * $nsperframe) + $nsoffset]]
draw text {0 0 0 } "Time = $time" size 0.7
#pbc box_draw -center com -centersel "segname P9 P11 P14  P9 P6 P13 P19 P10 P3 P4" -style dashed -width 0.1
}
####

