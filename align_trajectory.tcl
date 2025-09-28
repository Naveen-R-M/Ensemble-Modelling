# align_trajectory.tcl
# Usage:
#   vmd -dispdev text -e align_trajectory.tcl -args IN.pdb OUT.pdb

# ---- args ----
set in_pdb  [lindex $argv 0]
set out_pdb [lindex $argv 1]
if {[llength $argv] < 2} {
    puts stderr "Usage: vmd -dispdev text -e align_trajectory.tcl -args IN.pdb OUT.pdb"
    exit 1
}

# ---- load multi-model PDB ----
mol new $in_pdb type pdb first 0 last -1 waitfor all
set nframes [molinfo top get numframes]
puts "numframes=$nframes"

# ---- selections ----
set ref_frame 0
set seltext "protein and backbone"
set ref_sel   [atomselect top $seltext frame $ref_frame]
set check_sel [atomselect top $seltext]
set mov_sel   [atomselect top all]

# ---- align all frames to reference ----
for {set i 0} {$i < $nframes} {incr i} {
    $check_sel frame $i
    set M [measure fit $check_sel $ref_sel]
    $mov_sel frame $i
    $mov_sel move $M
}

# ---- write ALL frames with MODEL/ENDMDL ----
set fh [open $out_pdb "w"]
for {set i 0} {$i < $nframes} {incr i} {
    puts $fh [format "MODEL     %d" [expr {$i+1}]]

    # write current frame to a temp file
    set tmpfile [format "%s.frame%05d.pdb" $out_pdb $i]
    $mov_sel frame $i
    $mov_sel writepdb $tmpfile

    # append only atom records (strip any MODEL/END/TER/REMARK lines)
    set pf [open $tmpfile "r"]
    while {[gets $pf line] >= 0} {
        if {[regexp {^(CRYST1|REMARK|TER|END|ENDMDL|MODEL)} $line]} {
            continue
        }
        puts $fh $line
    }
    close $pf
    file delete -force $tmpfile

    puts $fh "ENDMDL"
}
close $fh
exit
