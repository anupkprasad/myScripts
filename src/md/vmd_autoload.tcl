# VMD auto-loader for myScripts MD Tcl modules.
# Sourced from ~/myScripts/callScripts.tcl
#
# Drop new *.tcl files under src/md/ and they load on the next VMD start
# (see excludes). Core utilities load first for dependency order.

namespace eval ::myscripts {
    variable ROOT
    variable DATA
    variable MD
    variable ANALYSIS
    variable VISUALIZATION
    variable SETUP
    variable BIOINFO
    variable UTILS
    variable REF
}

set ::myscripts::ROOT [file normalize [file join $::env(HOME) myScripts]]
if {![file isdirectory [file join $::myscripts::ROOT src md]]} {
    set ::myscripts::ROOT [file normalize [file join [file dirname [info script]] ../..]]
}

set ::myscripts::DATA          [file join $::myscripts::ROOT data]
set ::myscripts::MD            [file join $::myscripts::ROOT src md]
set ::myscripts::ANALYSIS      [file join $::myscripts::MD analysis]
set ::myscripts::VISUALIZATION [file join $::myscripts::MD visualization]
set ::myscripts::SETUP         [file join $::myscripts::MD setup namd]
set ::myscripts::BIOINFO       [file join $::myscripts::ROOT src bioinformatics]
set ::myscripts::UTILS         [file join $::myscripts::ROOT utils]
set ::myscripts::REF           [file join $::myscripts::DATA reference KinaseStructures]

set pathVmdScript  $::myscripts::ROOT
set pathAnalysis   $::myscripts::ANALYSIS
set pathGraphic    $::myscripts::VISUALIZATION
set pathSetup      $::myscripts::SETUP
set pathCentroid   $::myscripts::ANALYSIS
set refStructures  $::myscripts::REF
set pathBioinfo    $::myscripts::BIOINFO
set pathUtils      $::myscripts::UTILS

proc ::myscripts::_should_skip_tcl {path} {
    set base [file tail $path]
    set norm [file normalize $path]
    if {$base eq "vmd_autoload.tcl"} { return 1 }
    if {[string match "*copy*" $base]} { return 1 }
    if {[string match "test_*" $base]} { return 1 }
    if {[string match "*_test.tcl" $base]} { return 1 }
    if {[string match "*simple_test*" $base]} { return 1 }
    if {[string match "*/jobTemplate/*" $norm]} { return 1 }
    if {[string match "*/proteinSDS/*" $norm]} { return 1 }
    if {[string match "*/refsParaInp/*" $norm]} { return 1 }
    if {[string match "*/setup/namd/protein/*" $norm]} { return 1 }
    if {[string match "*/setup/namd/*" $norm]} { return 1 }
    if {[string match "*/analysis/oneshot/*" $norm]} { return 1 }
    return 0
}

proc ::myscripts::_collect_tcl {root} {
    set files {}
    if {![file isdirectory $root]} { return $files }
    set stack [list $root]
    while {[llength $stack]} {
        set dir [lindex $stack 0]
        set stack [lrange $stack 1 end]
        foreach f [lsort [glob -nocomplain -directory $dir *]] {
            if {[file isdirectory $f]} {
                lappend stack $f
            } elseif {[string match *.tcl $f]} {
                if {![::myscripts::_should_skip_tcl $f]} {
                    lappend files [file normalize $f]
                }
            }
        }
    }
    return [lsort -unique $files]
}

proc ::myscripts::source_file {path} {
    if {[catch {uplevel #0 [list source $path]} err]} {
        puts "WARNING: failed to source $path: $err"
        return 0
    }
    return 1
}

set ::myscripts::_loaded 0

# Core MD helpers first; structure loaders live in utils/parseFiles.tcl
set core [list \
    [file join $::myscripts::ANALYSIS basicFunctions.tcl] \
    [file join $::myscripts::ANALYSIS trajFunctions.tcl] \
    [file join $::myscripts::UTILS parseFiles.tcl] \
]

set core_norm {}
foreach f $core {
    if {[file exists $f]} {
        if {[::myscripts::source_file $f]} { incr ::myscripts::_loaded }
        lappend core_norm [file normalize $f]
    }
}

foreach f [::myscripts::_collect_tcl $::myscripts::MD] {
    if {[lsearch -exact $core_norm $f] >= 0} { continue }
    if {[::myscripts::source_file $f]} { incr ::myscripts::_loaded }
}

foreach f [list \
    [file join $::myscripts::SETUP systemSetup.tcl] \
    [file join $::myscripts::SETUP phosPTM.tcl] \
] {
    if {[file exists $f]} {
        if {[::myscripts::source_file $f]} { incr ::myscripts::_loaded }
    }
}

foreach f [::myscripts::_collect_tcl $::myscripts::BIOINFO] {
    if {[::myscripts::source_file $f]} { incr ::myscripts::_loaded }
}

# Any additional top-level utils/*.tcl (parseFiles already in core)
foreach f [lsort [glob -nocomplain -directory $::myscripts::UTILS *.tcl]] {
    set norm [file normalize $f]
    if {[lsearch -exact $core_norm $norm] >= 0} { continue }
    if {[::myscripts::source_file $norm]} { incr ::myscripts::_loaded }
}

puts "myScripts VMD: loaded $::myscripts::_loaded Tcl modules (md + bioinformatics + utils)"
puts "  ROOT=$::myscripts::ROOT"
puts "  DATA=$::myscripts::DATA  REF=$::myscripts::REF"
