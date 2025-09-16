# extract_glycans_and_chains.tcl

proc die {msg} { puts stderr $msg ; exit 1 }

if {[llength $argv] < 4} {
    die "Usage: vmd -dispdev text -e extract_glycans_and_chains.tcl -args input.pdb {s1 e1} {s2 e2} {s3 e3}"
}
set filename [lindex $argv 0]
set g1_range [lindex $argv 1]
set g2_range [lindex $argv 2]
set g3_range [lindex $argv 3]

# Validate: each is a list of two integers
foreach rname {g1_range g2_range g3_range} {
    set r [set $rname]
    if {[catch {set n [llength $r]}] || $n != 2} {
        die "Bad $rname='$r' (expect {start end})"
    }
    foreach v $r {
        if {![string is integer -strict $v]} {
            die "Bad value in $rname='$r' (non-integer '$v')"
        }
    }
}
lassign $g1_range g1_start g1_end
lassign $g2_range g2_start g2_end
lassign $g3_range g3_start g3_end

puts ">>> Processing file: $filename"
puts ">>> Glycan ranges: G1={$g1_start $g1_end}  G2={$g2_start $g2_end}  G3={$g3_start $g3_end}"

mol new $filename type pdb waitfor all

set sel_g1 [atomselect top "not protein and (resid $g1_start to $g1_end)"]
set sel_g2 [atomselect top "not protein and (resid $g2_start to $g2_end)"]
set sel_g3 [atomselect top "not protein and (resid $g3_start to $g3_end)"]

$sel_g1 writepdb CAR1.pdb
$sel_g2 writepdb CAR2.pdb
$sel_g3 writepdb CAR3.pdb

set sel_chA [atomselect top "protein and chain A"]; $sel_chA writepdb CHA.pdb
set sel_chB [atomselect top "protein and chain B"]; $sel_chB writepdb CHB.pdb
set sel_chC [atomselect top "protein and chain C"]; $sel_chC writepdb CHC.pdb
set sel_chD [atomselect top "protein and chain D"]; $sel_chD writepdb CHD.pdb
set sel_chE [atomselect top "protein and chain E"]; $sel_chE writepdb CHE.pdb
set sel_chF [atomselect top "protein and chain F"]; $sel_chF writepdb CHF.pdb

foreach s {sel_g1 sel_g2 sel_g3 sel_chA sel_chB sel_chC sel_chD sel_chE sel_chF} { [set $s] delete }
mol delete top
exit
