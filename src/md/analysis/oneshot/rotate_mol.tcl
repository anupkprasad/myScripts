set sel [atomselect top "protein"]
set com [measure center $sel weight mass]
set matrix [transaxis z 90]
$sel moveby [vecscale -1.0 $com]
$sel move $matrix
$sel moveby $com
