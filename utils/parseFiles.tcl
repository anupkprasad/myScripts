# parseFiles.tcl — all structure / trajectory loaders for VMD (utils/)
#
# CapitalCase entry: LoadStruct
#   LoadStruct help
#
# Autoloaded from ~/myScripts/utils via src/md/vmd_autoload.tcl

proc LoadStruct {args} {
    if {[llength $args] == 0} {
        _ls_help
        return
    }
    set cmd [lindex $args 0]
    set rest [lrange $args 1 end]
    if {$cmd in {help --help -help -h}} {
        _ls_help
        return
    }
    switch -exact -- $cmd {
        pdb            { _ls_cmd_pdb $rest }
        cif            { _ls_cmd_cif $rest }
        dcd            { _ls_cmd_dcd $rest }
        xtc            { _ls_cmd_xtc $rest }
        af3            { _ls_cmd_af3 $rest }
        af3-add        { _ls_cmd_af3_add $rest }
        af3-multimol   { _ls_cmd_af3_multimol $rest }
        af3-multimol-lig { _ls_cmd_af3_multimol_lig $rest }
        default {
            puts "LoadStruct: unknown command '$cmd'"
            _ls_help
        }
    }
}

# Short alias for AF3-only workflows (same file — not a separate module)
proc LoadAf3 {args} {
    if {[llength $args] == 0 || [lindex $args 0] in {help --help -help -h}} {
        _ls_help_af3
        return
    }
    set cmd [lindex $args 0]
    set rest [lrange $args 1 end]
    switch -exact -- $cmd {
        cif         { LoadStruct af3 {*}$rest }
        add         { LoadStruct af3-add {*}$rest }
        multimol    { LoadStruct af3-multimol {*}$rest }
        multimol-lig { LoadStruct af3-multimol-lig {*}$rest }
        default {
            # allow: LoadAf3 -pattern seed  →  LoadStruct af3 ...
            LoadStruct af3 {*}$args
        }
    }
}

proc _ls_help {} {
    puts {
LoadStruct — load structures / trajectories into VMD

Usage:
  LoadStruct <command> [options]
  LoadStruct help
  LoadAf3 <command> [options]     ;# AF3 shortcuts (same module)

Commands:
  pdb              Load PDB file(s) as frames of one molecule
  cif              Load CIF file(s) as frames of one molecule
  dcd              Load ionized.psf + prod*/.dcd under a data path
  xtc              Load GRO topology + XTC trajectory
  af3              Load AF3 seed-*_sample-* CIFs as one multi-frame mol
  af3-add          Append AF3 CIFs onto the current top molecule
  af3-multimol     One molecule per AF3 job directory
  af3-multimol-lig multimol + ligand rep / name heuristics
  help             Show this message

Common options (pdb / cif):
  -files <list>    explicit file list (default: glob *.pdb / *.cif in -dir)
  -dir <path>      directory to glob (default: pwd)
  -multimol 0|1    1 = one molecule per file (default: 0, frames on one mol)

dcd options:
  -dir <path>      data path containing ionized.psf and prod* dirs
  -pattern <text>  optional subdirectory filter (default: prod/Prod)
  -stride <n>      frame stride (default: 10)

xtc options:
  -top <file>      topology GRO/PDB (default: ../equilibration/min.gro)
  -traj <file>     XTC path (default: mdWrap1ns.xtc)
  -stride <n>      frame stride (default: 1)

AF3 options: see LoadAf3 help  (or LoadStruct help-af3 via LoadAf3 help)

Python (same stem, MDAnalysis):
  from utils.parseFiles import LoadUniverse, LoadAf3
  u = LoadUniverse("top.pdb", "traj.xtc")
  u = LoadAf3(".")   # seed-*_sample-* → multi-frame Universe

Examples:
  LoadStruct pdb
  LoadStruct cif -dir ./models
  LoadStruct dcd -dir /path/to/sim -stride 5
  LoadStruct xtc -top ../equilibration/min.gro -traj md.xtc
  LoadStruct af3 -rename af3_ensemble
  LoadAf3 cif -pattern seed
}
}

proc _ls_help_af3 {} {
    puts {
LoadAf3 — AF3 shortcuts into LoadStruct (utils/parseFiles.tcl)

AF3 folders look like:  seed-<seed>_sample-<1..5>/*.cif
Loading them on one mol → multi-frame trajectory.

Usage:
  LoadAf3 cif|add|multimol|multimol-lig [options]
  LoadAf3 help

Options:
  -dir <path>       (default: pwd)
  -pattern <glob>   directory prefix (default: seed)
  -rename <name>    molecule name (cif / add)
  -prefix <text>    job-folder prefix (multimol*; required)
  -ref <pdb>        reference under data/reference/KinaseStructures
  -sel <text>       KD / align selection
  -align 0|1        align to -ref (default: 1 if -ref set)

Examples:
  LoadAf3 cif -rename af3_ensemble
  LoadAf3 multimol -prefix DCLK3_ -ref DCLK3_KD.pdb
}
}

proc _ls_args {argv allowed defaults} {
    array set opt $defaults
    set n [llength $argv]
    for {set i 0} {$i < $n} {incr i} {
        set tok [lindex $argv $i]
        if {![string match -* $tok]} {
            error "LoadStruct: unexpected '$tok'. Type: LoadStruct help"
        }
        set name [string range $tok 1 end]
        if {$name in {h help -help}} {
            _ls_help
            set opt(_help) 1
            return [array get opt]
        }
        if {[lsearch -exact $allowed $name] < 0} {
            error "LoadStruct: unknown option -$name"
        }
        incr i
        if {$i >= $n} {
            error "LoadStruct: -$name needs a value"
        }
        set opt($name) [lindex $argv $i]
    }
    return [array get opt]
}

# ---------- PDB / CIF ----------
proc _ls_cmd_pdb {argv} {
    array set opt [_ls_args $argv {dir files multimol} {
        dir .
        files {}
        multimol 0
    }]
    if {[info exists opt(_help)]} { return }
    _ls_load_coords pdb $opt(dir) $opt(files) $opt(multimol)
}

proc _ls_cmd_cif {argv} {
    array set opt [_ls_args $argv {dir files multimol} {
        dir .
        files {}
        multimol 0
    }]
    if {[info exists opt(_help)]} { return }
    _ls_load_coords cif $opt(dir) $opt(files) $opt(multimol)
}

proc _ls_load_coords {ext dir files multimol} {
    set dir [file normalize $dir]
    if {$files eq ""} {
        set files [lsort -dictionary [glob -nocomplain -directory $dir *.$ext]]
    }
    if {[llength $files] == 0} {
        error "LoadStruct: no .$ext files in $dir"
    }
    set i 0
    foreach f $files {
        puts "LoadStruct: $f"
        if {$multimol} {
            mol new $f waitfor all
        } elseif {$i == 0} {
            mol new $f waitfor all
        } else {
            mol addfile $f waitfor all
        }
        incr i
    }
    if {!$multimol} {
        puts "LoadStruct: molecule [molinfo top] has [molinfo top get numframes] frames"
    }
}

# ---------- DCD / XTC ----------
proc _ls_cmd_dcd {argv} {
    array set opt [_ls_args $argv {dir pattern stride} {
        dir .
        pattern None
        stride 10
    }]
    if {[info exists opt(_help)]} { return }
    set datapath [file normalize $opt(dir)]
    set psf_list [glob -nocomplain -directory $datapath ionized.psf]
    if {[llength $psf_list] == 0} {
        error "LoadStruct dcd: ionized.psf not found under $datapath"
    }
    set psf_name [lindex $psf_list 0]
    puts "LoadStruct: $psf_name"
    mol new $psf_name waitfor all
    if {$opt(pattern) eq "None"} {
        set prod_list [glob -nocomplain -directory $datapath -types d *{prod,Prod}*]
    } else {
        set prod_list [glob -nocomplain -directory $datapath -types d *$opt(pattern)*]
    }
    set prod_list [lsort -dictionary $prod_list]
    if {[llength $prod_list] == 0} {
        error "LoadStruct dcd: no prod* directories under $datapath"
    }
    foreach prod $prod_list {
        set dcds [glob -nocomplain -directory $prod *.dcd]
        if {[llength $dcds] == 0} { continue }
        set filename [lindex [lsort -dictionary $dcds] 0]
        puts "LoadStruct: $filename"
        mol addfile $filename first 0 last -1 step $opt(stride) waitfor all
    }
    puts "LoadStruct: DCD trajectory loaded ([molinfo top get numframes] frames)"
}

proc _ls_cmd_xtc {argv} {
    array set opt [_ls_args $argv {top traj stride} {
        top ../equilibration/min.gro
        traj mdWrap1ns.xtc
        stride 1
    }]
    if {[info exists opt(_help)]} { return }
    if {![file exists $opt(top)]} {
        error "LoadStruct xtc: topology not found: $opt(top)"
    }
    if {![file exists $opt(traj)]} {
        error "LoadStruct xtc: trajectory not found: $opt(traj)"
    }
    mol new $opt(top) waitfor all
    puts "LoadStruct: topology $opt(top)"
    mol addfile $opt(traj) first 0 last -1 step $opt(stride) waitfor all
    puts "LoadStruct: XTC $opt(traj) ([molinfo top get numframes] frames)"
}

# ---------- AF3 ----------
proc _ls_seed_dirs {dir pattern} {
    set dir [file normalize $dir]
    if {![file isdirectory $dir]} {
        error "LoadStruct: not a directory: $dir"
    }
    set matches [glob -nocomplain -directory $dir -types d ${pattern}*]
    if {[llength $matches] == 0} {
        error "LoadStruct: no directories matching '${pattern}*' under $dir"
    }
    return [lsort -dictionary $matches]
}

proc _ls_cif_in_dir {seed_dir} {
    set cifs [lsort -dictionary [glob -nocomplain -directory $seed_dir *.cif]]
    if {[llength $cifs] == 0} {
        error "LoadStruct: no .cif in $seed_dir"
    }
    return [lindex $cifs 0]
}

proc _ls_af3_frames {dir pattern {append 0} {rename ""}} {
    set seed_dirs [_ls_seed_dirs $dir $pattern]
    set n [llength $seed_dirs]
    set first 1
    set i 0
    foreach d $seed_dirs {
        set cif [_ls_cif_in_dir $d]
        incr i
        puts "LoadStruct af3: \[$i/$n\] [file tail $d] <- [file tail $cif]"
        if {$first && !$append} {
            mol new $cif waitfor all
            set first 0
        } else {
            mol addfile $cif waitfor all
            set first 0
        }
    }
    set molid [molinfo top]
    if {$rename ne ""} {
        mol rename $molid $rename
    }
    puts "LoadStruct af3: molecule $molid has [molinfo $molid get numframes] frames"
    return $molid
}

proc _ls_cmd_af3 {argv} {
    array set opt [_ls_args $argv {dir pattern rename} {
        dir .
        pattern seed
        rename {}
    }]
    if {[info exists opt(_help)]} { return }
    _ls_af3_frames $opt(dir) $opt(pattern) 0 $opt(rename)
}

proc _ls_cmd_af3_add {argv} {
    array set opt [_ls_args $argv {dir pattern rename} {
        dir .
        pattern seed
        rename {}
    }]
    if {[info exists opt(_help)]} { return }
    if {[molinfo num] == 0} {
        error "LoadStruct af3-add: no molecule loaded; use LoadStruct af3 first"
    }
    _ls_af3_frames $opt(dir) $opt(pattern) 1 $opt(rename)
}

proc _ls_job_dirs {parent prefix} {
    if {$prefix eq ""} {
        error "LoadStruct: -prefix is required"
    }
    set parent [file normalize $parent]
    set entries [glob -nocomplain -directory $parent -types d *]
    set directories {}
    foreach entry $entries {
        if {[string match ${prefix}* [file tail $entry]]} {
            lappend directories $entry
        }
    }
    if {[llength $directories] == 0} {
        error "LoadStruct: no directories matching '${prefix}*' under $parent"
    }
    return [lsort -dictionary $directories]
}

proc _ls_maybe_load_ref {refStr sel} {
    if {$refStr eq ""} { return "" }
    global refStructures
    set path [file join $refStructures $refStr]
    if {![file exists $path]} {
        error "LoadStruct: reference not found: $path"
    }
    mol new $path waitfor all
    set refid [molinfo top]
    mol modselect 0 $refid "protein and $sel"
    mol modcolor 0 $refid ColorID 8
    mol selection "protein and not ($sel)"
    mol representation NewCartoon
    mol color ColorID 2
    mol addrep $refid
    puts "LoadStruct: reference mol $refid <- $path"
    return $refid
}

proc _ls_cmd_af3_multimol {argv} {
    array set opt [_ls_args $argv {dir prefix ref sel align pattern} {
        dir .
        prefix {}
        ref {}
        sel {resid 390 to 648 and name CA}
        align {}
        pattern seed
    }]
    if {[info exists opt(_help)]} { return }
    if {$opt(align) eq ""} {
        set opt(align) [expr {$opt(ref) ne "" ? 1 : 0}]
    }
    set refid [_ls_maybe_load_ref $opt(ref) $opt(sel)]
    set job_dirs [_ls_job_dirs $opt(dir) $opt(prefix)]
    set n [llength $job_dirs]
    set i 0
    foreach job $job_dirs {
        incr i
        puts "LoadStruct af3-multimol: \[$i/$n\] [file tail $job]"
        _ls_af3_frames $job $opt(pattern) 0 [file tail $job]
        set molid [molinfo top]
        if {$refid ne "" && $opt(align)} {
            trajAlignTwoMol $refid $molid "backbone and $opt(sel)" "backbone and $opt(sel)"
        }
        mol modselect 0 $molid "protein and $opt(sel)"
        mol selection "protein and not ($opt(sel))"
        mol representation NewCartoon
        mol color ColorID 4
        mol addrep $molid
        mol selection "protein and chain B"
        mol representation NewCartoon
        mol color ColorID 10
        mol addrep $molid
        mol off $molid
    }
    puts "LoadStruct af3-multimol: loaded $n molecules"
}

proc _ls_cmd_af3_multimol_lig {argv} {
    array set opt [_ls_args $argv {dir prefix ref sel align pattern} {
        dir .
        prefix {}
        ref {}
        sel {}
        align {}
        pattern seed
    }]
    if {[info exists opt(_help)]} { return }
    if {$opt(align) eq ""} {
        set opt(align) [expr {$opt(ref) ne "" ? 1 : 0}]
    }
    set refid [_ls_maybe_load_ref $opt(ref) $opt(sel)]
    set job_dirs [_ls_job_dirs $opt(dir) $opt(prefix)]
    set n [llength $job_dirs]
    set i 0
    foreach job $job_dirs {
        incr i
        set dirname [file tail $job]
        puts "LoadStruct af3-multimol-lig: \[$i/$n\] $dirname"
        set parts [split $dirname _]
        set rename [lindex $parts 1]_[lindex $parts 2]
        if {$rename eq "_"} { set rename $dirname }
        _ls_af3_frames $job $opt(pattern) 0 $rename
        set molid [molinfo top]
        if {$refid ne "" && $opt(align) && $opt(sel) ne ""} {
            trajAlignTwoMol $refid $molid "backbone and $opt(sel)" "backbone and $opt(sel)"
        }
        if {[llength [info procs trajAlign]]} {
            trajAlign
        }
        mol selection "not protein"
        mol representation Licorice
        mol color ColorID 4
        mol addrep $molid
        mol off $molid
    }
    puts "LoadStruct af3-multimol-lig: loaded $n molecules"
}

# ---------- legacy aliases ----------
proc trajLoadDCD {{datapath .} {pattern "None"} {step 10}} {
    puts "NOTE: trajLoadDCD deprecated; use: LoadStruct dcd -dir …"
    LoadStruct dcd -dir $datapath -pattern $pattern -stride $step
}

proc trajLoadXtc {{filename "mdWrap1ns.xtc"} {step 1}} {
    puts "NOTE: trajLoadXtc deprecated; use: LoadStruct xtc -traj …"
    LoadStruct xtc -traj $filename -stride $step
}

proc loadPdbs {{pdb_list ""}} {
    puts "NOTE: loadPdbs deprecated; use: LoadStruct pdb"
    if {$pdb_list eq ""} {
        LoadStruct pdb
    } else {
        LoadStruct pdb -files $pdb_list
    }
}

proc loadPdbsMultimol {{pdb_list ""}} {
    puts "NOTE: loadPdbsMultimol deprecated; use: LoadStruct pdb -multimol 1"
    if {$pdb_list eq ""} {
        LoadStruct pdb -multimol 1
    } else {
        LoadStruct pdb -files $pdb_list -multimol 1
    }
}

proc loadCifs {{cif_list ""}} {
    puts "NOTE: loadCifs deprecated; use: LoadStruct cif"
    if {$cif_list eq ""} {
        LoadStruct cif
    } else {
        LoadStruct cif -files $cif_list
    }
}

proc loadCifAF3 {{pattern seed}} {
    puts "NOTE: loadCifAF3 deprecated; use: LoadStruct af3 / LoadAf3 cif"
    LoadStruct af3 -pattern $pattern
}

proc addCifAF3 {{rename ""}} {
    puts "NOTE: addCifAF3 deprecated; use: LoadStruct af3-add"
    if {$rename ne ""} {
        LoadStruct af3-add -rename $rename
    } else {
        LoadStruct af3-add
    }
}

proc loadCifAF3MultiMolecules {{dirPrefix ""} {refStr ""} {KDResids "resid 390 to 648 and name CA"}} {
    puts "NOTE: loadCifAF3MultiMolecules deprecated; use: LoadStruct af3-multimol"
    set cmd [list LoadStruct af3-multimol -prefix $dirPrefix -sel $KDResids]
    if {$refStr ne ""} { lappend cmd -ref $refStr }
    {*}$cmd
}

proc loadCifAF3MultiMoleculesLigand {{dirPrefix ""} {refStr ""} {KDResids ""}} {
    puts "NOTE: loadCifAF3MultiMoleculesLigand deprecated; use: LoadStruct af3-multimol-lig"
    set cmd [list LoadStruct af3-multimol-lig -prefix $dirPrefix]
    if {$refStr ne ""} { lappend cmd -ref $refStr }
    if {$KDResids ne ""} { lappend cmd -sel $KDResids }
    {*}$cmd
}

proc readFileAsList {filename} {
    set fp [open $filename r]
    set lines {}
    while {[gets $fp line] >= 0} {
        lappend lines $line
    }
    close $fp
    return $lines
}

puts "LoadStruct loaded. Type: LoadStruct help   (AF3: LoadAf3 help)"
