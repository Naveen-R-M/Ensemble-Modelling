# align_trajectory.tcl — minimal, quiet, uses "protein and backbone" with autobonds ON
# Inputs come from environment variables (robust, no -args needed):
#   IN_PDB        : input PDB path  (required)
#   OUT_PDB       : output PDB path (required)
#   ALIGN_SEL_ENV : optional atom selection (defaults to "protein and backbone")

proc die {msg} { puts stderr "FATAL: $msg" ; flush stderr ; exit 1 }

# ---- read inputs (env only) ----
if {![info exists ::env(IN_PDB)] || ![info exists ::env(OUT_PDB)]} {
    die "Missing IN_PDB or OUT_PDB environment variables."
}
set inFile  $::env(IN_PDB)
set outFile $::env(OUT_PDB)
set seltext "protein and backbone"
if {[info exists ::env(ALIGN_SEL_ENV)] && $::env(ALIGN_SEL_ENV) ne ""} {
    set seltext $::env(ALIGN_SEL_ENV)
}

if {![file exists $inFile]} { die "Input PDB does not exist: $inFile" }

# ---- load and align ----
# Autobonds is ON by default; we rely on it so 'protein and backbone' works.
mol new $inFile type pdb first 0 last -1 waitfor all

set ref_sel   [atomselect top $seltext frame 0]
set check_sel [atomselect top $seltext]
set mov_sel   [atomselect top all]

if {[$ref_sel num] == 0} { die "Selection '$seltext' matched 0 atoms." }

set nf [molinfo top get numframes]
if {$nf <= 0} { die "No frames found in input: $inFile" }

for {set i 0} {$i < $nf} {incr i} {
    $check_sel frame $i
    set trans_mat [measure fit $check_sel $ref_sel]
    $mov_sel frame $i
    $mov_sel move $trans_mat
}

animate write pdb $outFile
exit
