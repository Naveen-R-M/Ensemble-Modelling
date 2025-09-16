# Usage from shell:
#   vmd -dispdev text -e align_traj.tcl -args INPUT_PDB OUTPUT_PDB [SEL]
# Example:
#   vmd -dispdev text -e align_traj.tcl -args out.pdb aligned.pdb "protein and backbone"

# --- parse args ---
set inFile  [lindex $argv 0]
set outFile [lindex $argv 1]
set seltext [lindex $argv 2]
if {$seltext eq ""} { set seltext "protein and backbone" }

# --- load the multi-model PDB as a trajectory ---
mol new $inFile type pdb first 0 last -1 waitfor all

# --- align all frames to frame 0 using the selection ---
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

# --- write aligned trajectory ---
animate write pdb $outFile
exit