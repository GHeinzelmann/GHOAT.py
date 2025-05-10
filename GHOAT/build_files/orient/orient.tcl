package provide Orient 1.0
package require La

namespace eval ::Orient:: {
    namespace export orient
}

# package require Orient
# namespace import Orient::orient
# ... load your molecules and make a selection ...
#
# set I [draw principalaxes $sel]           <--- show/calc the principal axes
# set A [orient $sel [lindex $I 2] {0 0 1}] <--- rotate axis 2 to match Z
# $sel move $A
# set I [draw principalaxes $sel]           <--- recalc principal axes to check
# set A [orient $sel [lindex $I 1] {0 1 0}] <--- rotate axis 1 to match Y
# $sel move $A
# set I [draw principalaxes $sel]           <--- recalc principal axes to check#
proc Orient::sel_com { sel weights } {
    set x [ $sel get x ]
    set y [ $sel get y ]
    set z [ $sel get z ]
    set m $weights
    
    set comx 0
    set comy 0
    set comz 0
    set totalm 0
    foreach xx $x yy $y zz $z mm $m {
        # use the abs of the weights
        set mm [expr abs($mm)]
	set comx [ expr "$comx + $xx*$mm" ]
	set comy [ expr "$comy + $yy*$mm" ]
	set comz [ expr "$comz + $zz*$mm" ]
	set totalm [ expr "$totalm + $mm" ]
    }
    set comx [ expr "$comx / $totalm" ]
    set comy [ expr "$comy / $totalm" ]
    set comz [ expr "$comz / $totalm" ]
    puts "Total weight: $totalm"
    return [list $comx $comy $comz]
}

proc Orient::sel_it { sel COM weights} {
    set x [ $sel get x ]
    set y [ $sel get y ]
    set z [ $sel get z ]
    set m $weights

    # compute I
    set Ixx 0
    set Ixy 0
    set Ixz 0
    set Iyy 0
    set Iyz 0
    set Izz 0
    foreach xx $x yy $y zz $z mm $m {
        # use the abs of the weights
        set mm [expr abs($mm)]
        
        # subtract the COM
        set xx [expr $xx - [lindex $COM 0]]
        set yy [expr $yy - [lindex $COM 1]]
        set zz [expr $zz - [lindex $COM 2]]

        set rr [expr $xx + $yy + $zz]

        set Ixx [expr $Ixx + $mm*($yy*$yy+$zz*$zz)]
        set Ixy [expr $Ixy - $mm*($xx*$yy)]
        set Ixz [expr $Ixz - $mm*($xx*$zz)]
        set Iyy [expr $Iyy + $mm*($xx*$xx+$zz*$zz)]
        set Iyz [expr $Iyz - $mm*($yy*$zz)]
        set Izz [expr $Izz + $mm*($xx*$xx+$yy*$yy)]

    }
    
    return [list 2 3 3 $Ixx $Ixy $Ixz $Ixy $Iyy $Iyz $Ixz $Iyz $Izz]
}

proc vmd_draw_arrow {mol start end} {
    set scaling [expr [veclength [vecsub $end $start]]/100]
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.8 [vecsub $end $start]]]
    graphics $mol cylinder $start $middle radius [expr 2*$scaling]
    puts [list cone $middle $end radius [expr 5*$scaling]]
    graphics $mol cone $middle $end radius [expr 5*$scaling]
}

proc vmd_draw_vector { mol pos val } {
    set end   [ vecadd $pos [ vecscale +1 $val ] ]
    vmd_draw_arrow $mol $pos $end
}    

# find the max of some numbers
proc Orient::max { args } {
    set maxval [lindex $args 0]
    foreach arg $args {
        if { $arg > $maxval } {
            set maxval $arg
        }
    }
    return $maxval
}

# draws the three principal axes
proc vmd_draw_principalaxes { mol sel {weights domass} } {
    if { $weights == "domass" } {
        set weights [ $sel get mass ]
    }

    set I [Orient::calc_principalaxes $sel $weights]
    set a1 [lindex $I 0]
    set a2 [lindex $I 1]
    set a3 [lindex $I 2]

    # find the size of the system
    set minmax [measure minmax $sel]
    set ranges [vecsub [lindex $minmax 1] [lindex $minmax 0]]
    set scale [expr .7*[Orient::max [lindex $ranges 0] \
                             [lindex $ranges 1] \
                             [lindex $ranges 2]]]
    set scale2 [expr 1.02 * $scale]

    # draw some nice vectors
    graphics $mol delete all
    graphics $mol color yellow
    set COM [Orient::sel_com $sel $weights]
    vmd_draw_vector $mol $COM [vecscale $scale $a1]
    vmd_draw_vector $mol $COM [vecscale $scale $a2]
    vmd_draw_vector $mol $COM [vecscale $scale $a3]

    graphics $mol color white
    graphics $mol text [vecadd $COM [vecscale $scale2 $a1]] "1"
    graphics $mol text [vecadd $COM [vecscale $scale2 $a2]] "2"
    graphics $mol text [vecadd $COM [vecscale $scale2 $a3]] "3"
    
    return [list $a1 $a2 $a3]
}

# returns the three principal axes
proc Orient::calc_principalaxes { sel {weights domass} } {
    puts "Calculating principal axes."
    if { $weights == "domass" } {
        set weights [ $sel get mass ]
    }

    puts "Getting the center-of-mass..."
    # get the COM
    set COM [Orient::sel_com $sel $weights]
    puts "Computing the inertia tensor..."
    # get the I
    set I [Orient::sel_it $sel $COM $weights]
    puts "Drawing the principal components..."
    La::mevsvd_br I evals
    # now $I holds in its columns the principal axes
    set a1 "[lindex $I 3] [lindex $I 6] [lindex $I 9]"
    set a2 "[lindex $I 4] [lindex $I 7] [lindex $I 10]"
    set a3 "[lindex $I 5] [lindex $I 8] [lindex $I 11]"

    return [list $a1 $a2 $a3]
}

# rotate a selection about its COM, taking <vector1> to <vector2>
# e.g.: orient $sel [lindex $I 2] {0 0 1}
# (this aligns the third principal axis with z)
proc Orient::orient { sel vector1 vector2 {weights domass}} {
    if { $weights == "domass" } {
        set weights [ $sel get mass ]
    }

    set COM [Orient::sel_com $sel $weights]

    set vec1 [vecnorm $vector1]

    set vec2 [vecnorm $vector2]

    # compute the angle and axis of rotation
    set rotvec [veccross $vec1 $vec2]
    set sine   [veclength $rotvec]
    set cosine [vecdot $vec1 $vec2]
    set angle [expr atan2($sine,$cosine)]
    
    # return the rotation matrix
    return [trans center $COM axis $rotvec $angle rad]
}
