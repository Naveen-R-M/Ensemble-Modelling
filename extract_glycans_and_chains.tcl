# Usage:
# vmd -dispdev text -e extract_glycans_and_chains.tcl -args input.pdb g1s g1e g2s g2e g3s g3e

proc die {msg} {
    puts stderr $msg
    exit 1
}

if {[llength $argv] < 7} {
    die "Usage: vmd -dispdev text -e extract_glycans_and_chains.tcl -args input.pdb g1s g1e g2s g2e g3s g3e"
}

set filename  [lindex $argv 0]
set g1_start  [lindex $argv 1]
set g1_end    [lindex $argv 2]
set g2_start  [lindex $argv 3]
set g2_end    [lindex $argv 4]
set g3_start  [lindex $argv 5]
set g3_end    [lindex $argv 6]

# Validate integers
foreach v {g1_start g1_end g2_start g2_end g3_start g3_end} {
    if {![string is integer -strict [set $v]]} {
        die "Bad integer for $v='[set $v]'"
    }
}

puts ">>> Processing file: $filename"
puts ">>> Glycan ranges: G1={$g1_start $g1_end}  G2={$g2_start $g2_end}  G3={$g3_start $g3_end}"

mol new $filename type pdb waitfor all

set sel_all [atomselect top all]
$sel_all set beta 0

# Glycan selections
set sel_g1 [atomselect top "not protein and (resid $g1_start to $g1_end)"]
set sel_g2 [atomselect top "not protein and (resid $g2_start to $g2_end)"]
set sel_g3 [atomselect top "not protein and (resid $g3_start to $g3_end)"]

# Set segname and write glycans
$sel_g1 set segname CAR1
$sel_g1 writepdb CAR1.pdb

$sel_g2 set segname CAR2
$sel_g2 writepdb CAR2.pdb

$sel_g3 set segname CAR3
$sel_g3 writepdb CAR3.pdb

# Define standard amino acids once
set std_amino_acids "resname ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR"

# Protein chain selections
set sel_chA [atomselect top "chain A and protein"]
$sel_chA set segname CHA
$sel_chA writepdb CHA.pdb

set sel_chB [atomselect top "chain B and protein"]
$sel_chB set segname CHB
$sel_chB writepdb CHB.pdb

set sel_chC [atomselect top "chain C and protein"]
$sel_chC set segname CHC
$sel_chC writepdb CHC.pdb

set sel_chD [atomselect top "chain D and protein"]
$sel_chD set segname CHD
$sel_chD writepdb CHD.pdb

set sel_chE [atomselect top "chain E and protein"]
$sel_chE set segname CHE
$sel_chE writepdb CHE.pdb

set sel_chF [atomselect top "chain F and protein"]
$sel_chF set segname CHF
$sel_chF writepdb CHF.pdb

# Clean up selections
foreach s {sel_g1 sel_g2 sel_g3 sel_chA sel_chB sel_chC sel_chD sel_chE sel_chF} {
    [set $s] delete
}

mol delete top
exit
