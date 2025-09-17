# align_trajectory.tcl
# Usage: vmd -dispdev text -e align_trajectory.tcl -args INPUT_PDB OUTPUT_PDB [SEL]

proc die {msg} { puts stderr $msg ; exit 1 }

if {[llength $argv] < 2} {
    die "Usage: vmd -dispdev text -e align_trajectory.tcl -args INPUT_PDB OUTPUT_PDB [SEL]"
}

set inFile  [lindex $argv 0]
set outFile [lindex $argv 1]
set seltext [expr {[llength $argv] >= 3 ? [lindex $argv 2] : "protein and backbone"}]

mol new $inFile type pdb first 0 last -1 waitfor all
set ref_frame 0
set ref_sel   [atomselect top "$seltext" frame $ref_frame]
set check_sel [atomselect top "$seltext"]
set mov_sel   [atomselect top all]

set num_frames [molinfo top get numframes]
for {set i 0} {$i < $num_frames} {incr i} {
    $check_sel frame $i
    set trans_mat [measure fit $check_sel $ref_sel]
    $mov_sel frame $i
    $mov_sel move $trans_mat
}

animate write pdb $outFile
exit