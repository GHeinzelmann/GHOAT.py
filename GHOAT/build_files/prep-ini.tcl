mol load pdb hhhh-gggg.pdb
mol load pdb dum.pdb

package require Orient
namespace import Orient::orient

set sel [atomselect 0 "all"]
set hos [atomselect 0 "resid FIRST to LAST and noh"]
set I [draw principalaxes $hos]
set A [orient $hos [lindex $I 2] {0 0 1}]
$sel move $A
set I [draw principalaxes $hos]
set A [orient $hos [lindex $I 1] {0 1 0}]
$sel move $A
set I [draw principalaxes $hos]
$sel move [transaxis y -90]

set hos [atomselect 0 "resid FIRST to LAST and noh"]
$sel moveby [vecinvert [measure center $hos weight mass]]

set a [atomselect 0 "resname MMM and noh"]
set tot [$a get name]
set n 0
set leng 100
set lengmax 100
set ran RANG
set dmax DMAX
set dmin DMIN
set mat {}
set amx 90


foreach i $tot {
set t [atomselect 0 "resname MMM and name $i"]
set m [measure center $t weight mass]
foreach {x1 y1 z1} $m {break}
set xl $x1
set yl $y1
set zl $z1
if {[expr abs([expr $z1]) < $ran] && [expr sqrt([expr pow($x1,2) + pow($y1,2)]) < $ran ]} {
lappend mat $i
}
}

foreach i $mat {
set t [atomselect 0 "resname MMM and name $i"]
set d1 [measure center $t]
foreach {x1 y1 z1} $d1 {break}
set xl $x1
set yl $y1
set zl $z1
set leng1 [expr sqrt([expr pow($x1,2) + pow($y1,2)])]
puts $i
puts $leng1
if [expr $leng1 < $leng] {
set leng $leng1
set aa1 $i }
}


set exist [info exists aa1]
if {[expr $exist == 0]} {
set data ""
set filename "anchors.txt"
set fileId [open $filename "w"]
puts -nonewline $fileId $data
close $fileId
puts "Ligand first anchor not found"
exit
}

puts "" 
puts "anchor 1 is" 
puts $aa1
puts $leng
puts "" 

set amat {}
foreach i $tot {
set alis {}
set angle1 {}
set angle2 {}
set angle3 {}
set angle {}
set t [atomselect 0 "resname MMM and name $i"]
set p [atomselect 0 "resname MMM and name $aa1"]
set d [atomselect 0 "resid H1A and name NN1"]
if {$i ne $aa1} { set a1 [$d get index]
set d1 [measure center $t]
set d2 [measure center $p]
set leng [veclength [vecsub $d1 $d2]]
lappend angle1 $a1
lappend angle1 "0"
lappend angle $angle1
lappend alis [$d get name]
set a2 [$p get index]
lappend angle2 $a2
lappend angle2 "0"
lappend angle $angle2
lappend alis [$p get name]
set a3 [$t get index]
lappend angle3 $a3
lappend angle3 "0"
lappend angle $angle3
lappend alis [$t get name]
set ang [measure angle $angle]
puts $i
puts $ang
puts $leng
if {[expr $leng > $dmin] && [expr $leng < $dmax]} {
lappend amat $i}
}
}

puts "" 
foreach i $amat {
set t [atomselect 0 "resname MMM and name $i"]
set m [atomselect 0 "resname MMM and name $aa1"]
set d1 [measure center $t weight mass]
set d2 [measure center $m weight mass]
set diff [vecsub $d1 $d2]
foreach {x1 y1 z1} $diff {break}
set xl $x1
set yl $y1
set zl $z1
set leng1 [expr sqrt([expr pow($x1,2) + pow($y1,2)])]
puts $i
puts $leng1
if {[expr $leng1 < $leng]} {
set leng $leng1
set aa2 $i }
}

set exist [info exists aa2]
if {[expr $exist == 0]} {
set data "$aa1\n"
set filename "anchors.txt"
set fileId [open $filename "w"]
puts -nonewline $fileId $data
close $fileId
puts "Ligand second anchor not found"
exit
}

puts ""
puts "anchor 2 is" 
puts $aa2
puts $leng
puts ""

set amat {}
foreach i $tot {
set alis {}
set angle1 {}
set angle2 {}
set angle3 {}
set angle {}
set t [atomselect 0 "resname MMM and name $i"]
set p [atomselect 0 "resname MMM and name $aa2"]
set d [atomselect 0 "resname MMM and name $aa1"]
if {$i ne $aa1 && $i ne $aa2} { set a1 [$d get index]
set d1 [measure center $t]
set d2 [measure center $p]
set leng [veclength [vecsub $d1 $d2]]
lappend angle1 $a1
lappend angle1 "0"
lappend angle $angle1
lappend alis [$d get name]
set a2 [$p get index]
lappend angle2 $a2
lappend angle2 "0"
lappend angle $angle2
lappend alis [$p get name]
set a3 [$t get index]
lappend angle3 $a3
lappend angle3 "0"
lappend angle $angle3
lappend alis [$t get name]
set ang [measure angle $angle]
puts $i
puts $ang
puts $leng
if {[expr $leng > $dmin] && [expr $leng < $dmax]} {
lappend amat $i}
}
}

puts "" 
set amx 90
foreach i $amat {
puts $i
set angle1 {}
set angle2 {}
set angle3 {}
set angle {}
set t [atomselect 0 "resname MMM and name $i"]
set p [atomselect 0 "resname MMM and name $aa2"]
set d [atomselect 0 "resname MMM and name $aa1"]
set d1 [measure center $t weight mass]
set d2 [measure center $p weight mass]
set a1 [$d get index]
lappend angle1 $a1
lappend angle1 "0"
lappend angle $angle1
lappend alis [$d get name]
set a2 [$p get index]
lappend angle2 $a2
lappend angle2 "0"
lappend angle $angle2
lappend alis [$p get name]
set a3 [$t get index]
lappend angle3 $a3
lappend angle3 "0"
lappend angle $angle3
lappend alis [$t get name]
set ang [measure angle $angle]
if {[expr abs([expr $ang - 90.0])] < $amx} {
set amx [expr abs([expr $ang - 90.0])]
set angl $ang
set aa3 $i
set leng [veclength [vecsub $d1 $d2]]
}
}

set exist [info exists aa3]
if {[expr $exist == 0]} {
set data "$aa1 $aa2\n"
set filename "anchors.txt"
set fileId [open $filename "w"]
puts -nonewline $fileId $data
close $fileId
puts "Ligand third anchor not found"
exit
}

puts ""
puts "anchor 3 is" 
puts $aa3
puts $angl 
puts $leng

puts "" 
puts "The three anchors are" 
puts "$aa1 $aa2 $aa3"
puts "" 

set data "$aa1 $aa2 $aa3\n"
set filename "anchors.txt"
set fileId [open $filename "w"]
puts -nonewline $fileId $data
close $fileId

set a [atomselect 0 "resid FIRST to LAST and noh"]
set b [atomselect 0 "resname MMM and noh"]
set c [atomselect 1 all]
$c moveby [vecsub [measure center $a weight mass] [measure center $c weight mass]]
$c writepdb dum1.pdb
$b moveby {0 0 30}
$c moveby [vecsub [measure center $b weight mass] [measure center $c weight mass]]
$c set resid 2
$c writepdb dum2.pdb
$b moveby {0 0 -30}

set all [atomselect 0 all]
$all writepdb hhhh-gggg-aligned.pdb

exit
