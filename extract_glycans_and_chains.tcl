# extract_glycans_and_chains.tcl
# Runs inside VMD with glycan parameters

# -------------------------------------------------------------------
# Usage (from bash):
# vmd -dispdev text -e extract_glycans_and_chains.tcl -args input.pdb "{g1_start g1_end}" "{g2_start g2_end}" "{g3_start g3_end}"
# -------------------------------------------------------------------

# Get arguments from bash
set filename [lindex $argv 0]
set g1_range [lindex $argv 1]
set g2_range [lindex $argv 2]
set g3_range [lindex $argv 3]

puts ">>> Processing file: $filename"
puts ">>> Glycan ranges: G1=$g1_range  G2=$g2_range  G3=$g3_range"

# Load structure
mol new $filename type pdb waitfor all

# Select glycans based on ranges
set sel_g1 [atomselect top "not protein and (resid [lindex $g1_range 0] to [lindex $g1_range 1])"]
set sel_g2 [atomselect top "not protein and (resid [lindex $g2_range 0] to [lindex $g2_range 1])"]
set sel_g3 [atomselect top "not protein and (resid [lindex $g3_range 0] to [lindex $g3_range 1])"]

# Write out selections for downstream analysis
$sel_g1 writepdb CAR1.pdb
$sel_g2 writepdb CAR2.pdb
$sel_g3 writepdb CAR3.pdb

# Example: also write out protein chains if needed

set sel_chA [atomselect top "protein and chain A"]
$sel_chA writepdb CHA.pdb

set sel_chB [atomselect top "protein and chain B"]
$sel_chB writepdb CHB.pdb

set sel_chC [atomselect top "protein and chain C"]
$sel_chC writepdb CHC.pdb

set sel_chD [atomselect top "protein and chain D"]
$sel_chD writepdb CHD.pdb

set sel_chE [atomselect top "protein and chain E"]
$sel_chE writepdb CHE.pdb

set sel_chF [atomselect top "protein and chain F"]
$sel_chF writepdb CHF.pdb

# Clean up
$sel_g1 delete
$sel_g2 delete
$sel_g3 delete
$sel_chA delete
$sel_chB delete
$sel_chC delete
$sel_chD delete
$sel_chE delete
$sel_chF delete

mol delete top
exit