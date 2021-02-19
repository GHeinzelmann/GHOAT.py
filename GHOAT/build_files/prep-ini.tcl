mol load pdb hhhh-gggg.pdb
mol load pdb dum.pdb

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
set m [atomselect 0 "(resid H1A and name NN1) or (resid H2A and name NN2) or (resid H3A and name NN3)"]
set d1 [measure center $t]
set d2 [measure center $m]
set leng1 [veclength [vecsub $d1 $d2]]
if {[expr $leng1 < $ran]} {
lappend mat $i
}
}


foreach i $mat {
set t [atomselect 0 "resname MMM and name $i"]
set m [atomselect 0 "(resid H1A and name NN1) or (resid H2A and name NN2) or (resid H3A and name NN3)"]
set d1 [measure center $t]
set d2 [measure center $m]
set leng1 [veclength [vecsub $d1 $d2]]
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
puts $i
set angle1 {}
set angle2 {}
set angle3 {}
set angle {}
set t [atomselect 0 "resname MMM and name $i"]
set p [atomselect 0 "resname MMM and name $aa1"]
set d [atomselect 0 "resid H1A and name NN1"]
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
set aa2 $i
set leng [veclength [vecsub $d1 $d2]]
}
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
puts $angl 
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

set a [atomselect 0 "resid H1A and name NN1"]
set b [atomselect 0 "resname MMM and name $aa1"]
set c [atomselect 0 all]
set d1 [measure center $a weight mass]
set d2 [measure center $b weight mass]
set diff [vecsub $d1 $d2]
set d [veclength $diff]
foreach {x1 y1 z1} $diff {break}
set x $x1
set y $y1
set z $z1
$c move [transvecinv "$x $y $z"]
$c move [transaxis y 90]
$c moveby [vecinvert [measure center $c weight mass]]

set a [atomselect 0 "resid FIRST to LAST and not water and not resname MMM and noh"]
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
