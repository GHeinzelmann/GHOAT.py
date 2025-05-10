# $Header: /usr/cvsroot/la/la.tcl,v 1.3 2001/11/14 13:27:35 hume Exp $
#
# La == Linear Algebra and similar math
# 
# Part of the La Package
# (C)Copyright 2001, Hume Integration Software
#
# $Log: la.tcl,v $
# Revision 1.3  2001/11/14 13:27:35  hume
# Changed Services to Software.
#
# Revision 1.2  2001/07/12 23:59:51  hume
# Added K proc and boosted performance bigtime!
#
# Revision 1.1.1.1  2001/07/12 15:50:16  hume
# Initial check-in.
# 
#
#

package provide La 1.0



#
# Package scope and addressed requirements:
#

# Procedures cover manipulation of vectors and matrices such as
# scaling, concatenation by rows or columns, subsetting by
# rows or columns, formatted printing, transpose, dot product,
# matrix multiplication, solution of linear equation sets, 
# matrix inversion, and eigenvalue/eigenvector solutions.
# The user can mix vectors and arrays in linear algebra operations.  
# The logic does reasonable conversion of types.
# Sophisticated operations such as evaluating a custom procedure
# against each element of a matrix are easily possible.
#
# There are typically two calls for a procedure.  A plainly 
# named procedure such as "transpose", and another procedure
# with the suffix "_br" such as "transpose_br".  The plainly
# named procedures expect their data arguments to be passed
# by value, which is the usual argument passing convention with
# Tcl programming.  The plain calls are designed for ease of
# interactive use, and in general have been coded to perform 
# more conversion of arguments, more error checking, and to
# trade off efficiency for convenience.  The "_br" procedures
# are intended for efficent use, with the "_br" indicating that
# data arguments are passed "By Reference" to avoid copying the
# data.  In Tcl, to pass by reference means that the caller 
# has the data in a named variable, and the caller passes the
# name of the variable instead of a copy of the data.  You can
# see that passing by reference becomes important for larger
# vectors and arrays.  The _br procedures in general assume 
# that the data arguments have the correct structure, so the
# caller may need to use the "promote", "demote", and "transpose" 
# calls to prepare arguments for a _br call.



#
# La operand formats - scalars, vectors, and arrays
# 
# Tcl list formats are used for better performance than arrays
# and better compatibility with C/C++ code.
# We add dimension information in a simple way at the front of the list.
#   
#  La operands are structured Tcl lists.
#  If the list length is 1 the operand is a scalar, eg., "6.02e23".
#  If the first element in the list is 2, the two following
#  elements provide the number of rows and the number of columns.
#  If the first element is 2:
#     If the number of columns is 0, the operand is a one dimensional
#     vector.  So {2 4 0 1 2 3 4} is the vector {1 2 3 4}.  Vectors
#     have the same structure as matrices so they can be efficiently
#     modified in place.
#     So a vector of length N has the representation:
#        {2 N 0 v0 v1 v2 ... vN-1}
#     The index into the underlying Tcl list for v[i] is
#        set index [expr 3 + $i]
#     A vector of length N is promoted to an N rows by 1 column
#     matrix if used in an operation where a matrix argument is
#     expected.
#     The transpose of a vector of length N is a 1 row by N columns
#     matrix.
#     If the number of rows and the number of columns are both positive,
#     the operand is a matrix.  The data elements follow the number of
#     rows and columns.  The data is concatenated as row0, row1, ...
#     so a 2-D matrix, a, has the representation:
#       2 R C a[0,0] a[0,1] a[0,2] ... a[0,C-1] a[1,0] a[1,1] ... a[1,C-1] ... a[R-1,C-1]  
#
#     The index into the underlying Tcl list for a[i,j] is
#           set index [expr 3 + $i * $C + $j]
#   
# If the first element in the list is > 2, the package assumes the list
# represents a higher dimensional operand.  Logic for higher dimension
# operands is not currently part of the package.  
# The candidate structure for 3-D data is:
#    3 P R C a[0,0,0] a[0,0,1] a[0,0,2] ... a[0,0,C-1] ... a[P-1,R-1,C-1]
#    P=planes R=rows, C=cols
#    An intuitive view is that the data is multiple 2-D planes of rows 
#    and columns listed in order from the 0 plane to the P-1 plane.  
#    The index into the underlying Tcl list for a[i,j,k] is
#            set index [expr 4 + $i*$R*$c + $j*$c + $k]
#
#
#
#  La operands, explained by example:
#
# <s>     a one element scalar (llength == 1)
# eg., 3.1415
# vectors should use the first dimension of 2
# 2 1 0 3.1415    ;# a vector of length 1
# 2 0 1 3.1415    ;# error - vector should use first dimension only
# 2 3 0 1 2 3     ;# vector {1 2 3}
# 2 3 0 1 2 3 4   ;# error - vector has additional element
# 2 0 0           ;# error - cannot distinguish 1-D or 2-D 
# 2 0 0 <s>       ;# error - simpler to have scalars always llength 1
# a 2D array of N rows, M columns
# 2 N M a[0,0] a[0,1] a[0,2] ... a[0,M] a[1,0] a[1,1] ... a[1,M] ... a[N,M]  
#

# getting started
#
# % source "la.tcl"           ;# advanced users will use "package require La"
# % namespace import La::*    ;# makes names like "La::show" useable as "show"
# set x {2 2 3 1 2 3 4 5 6}
# % show $x
# 1 2 3
# 4 5 6
# % show [transpose $x]
# 1 4
# 2 5
# 3 6
# % set tx [transpose $x]
# 2 3 2 1 4 2 5 3 6
# % show [join_cols $tx $tx]
# 1 4 1 4
# 2 5 2 5
# 3 6 3 6
# % show [moffset $x 0.2]
# 1.2 2.2 3.2
# 4.2 5.2 6.2
# % set v {2 3 0 1 2 3}
# % show [mmult $x $v]
# 14.0
# 32.0



############################### Package Source Code #####################

# Advanced Users:
# The package namespace name is defined here by default.
# You can specify a different name before using the package.
global La_namespace
if { ![info exists La_namespace] } { set La_namespace La }

namespace eval $La_namespace {

  namespace export \
dim dim_br \
dotprod dotprod_br \
join_cols join_cols_br \
join_rows join_rows_br\
lassign lassign_br \
madd mdiv mprod msub mat_binary_op mat_binary_op_br \
mathprec \
mcols mcols_br\
mdiag \
mdingdong \
mevsvd mevsvd_br \
mhilbert \
mident \
mlssvd \
mmult mmult_br \
mnorms mnorms_br mnormalize mnormalize_br \
mrange mrange_br \
mround mround_ij \
mrows mrows_br\
msolve msolve_br \
msum msum_br \
msvd msvd_br\
mat_unary_op madjust moffset mscale mat_unary_op_br \
promote promote_br demote demote_br\
show show_br \
transpose transpose_br \
vdiag vdiag_br \
vtrim vtrim_br

  global La_namespace
  if { $La_namespace != "La" } {
      # help advanced users by moving possible compiled code to
      # their optionally chosen namespace name.
      if { [info commands La::assign_br] == "La::assign_br" } {
          rename La:assign_br $La_namespace::assign_br
          }
      }

############################# dim ####################################
############################# dim ####################################
############################# dim ####################################
#  return the dimension of the argument
#   and to some degree verify a proper format
# {} - null
# 0 - scalar
# 1 - vector
# 2 - array
# N - higher
proc dim x {
    return [dim_br x]
    }

proc dim_br name_x {
   upvar $name_x x
   set llen [llength $x]
   if { $llen == 0 } {return {}}
   if { $llen == 1 && [string is double -strict $x] } { return 0 }
   set n [lindex $x 0]
   if { $n == 2 } {
       set d1 [lindex $x 1]
       set d2 [lindex $x 2]
       if {[catch {
           if { $d1 != int($d1) || $d2 != int($d2) } {error {} }
           }] } { error "improper La operand format" }
       if { $d2 == 0 } {
           if { $llen == [expr {$d1 + 3}] } {return 1 } ;# 2 n 0 <n>
           error "improper length of vector"
           }
       if { $llen == [expr {$d1 * $d2 + 3}] } { return 2 }
       error "improper length of matrix"
       }
    if { $n >= 3 && [string is integer -strict $n] } { 
        # TBD verify llen as prod of dims + dims + 1 
        return $n 
        }
    error "improper La operand format"
    }



############################ dotprod #################################
############################ dotprod #################################
############################ dotprod #################################
#
# Dot Product = sum over i, Ai * Bi
#
# Can work columns or rows in matrices
# because indexing increment is passed.
# Default arguments work with vectors, Nx1, or 1XN matrices
#
proc dotprod {a b {N {}} {a0index 3} {b0index 3} {ainc 1} {binc 1}} {
    return [dotprod_br a b $N $a0index $b0index $ainc $binc]
    }

#
# always returns a scalar result
#
proc dotprod_br {a_in b_in {N {}} {a0index 3} {b0index 3} {ainc 1} {binc 1}} {
    upvar $a_in a
    upvar $b_in b
    set z 0.0
    if { $N == {}} { ;# default length assumes vector, 1xN, or Nx1 first argument
        set len [llength $a]
        set N [expr {$len-3}]
        # and then see if that matches one of b's dimensions
        if { [lindex $b 1] != $N && [lindex $b 2] != $N } {
            error "Assumed length does not seem proper."
            }
        }
    for {set i 0} {$i < $N} {incr i} {
        set aval [lindex $a $a0index]
        set bval [lindex $b $b0index]
        set z [expr { $z + $aval * $bval}]
        incr a0index $ainc
        incr b0index $binc
        }
    return $z
    }


########################### join_cols ################################
########################### join_cols ################################
########################### join_cols ################################
#
# combine vector or matrix args as added columns
#
proc join_cols {a b} {
    set dimA [dim_br a]
    set dimB [dim_br b]
    if { $dimA > 2 || $dimB > 2 } { 
        error "cannot handle > 2D data"
        }
    if { $dimA < 2 } {promote_br a}
    if { $dimB < 2 } {promote_br b}
    join_cols_br a b
    return $a
    }    

    
proc join_cols_br {a_in b_in {c_out {}}} {
    upvar $a_in a
    upvar $b_in b
    if { $c_out == {}} {set c_out $a_in }
    upvar $c_out c
    # have to match row counts
    set rows [lindex $a 1]
    set acols [lindex $a 2]
    set bcols [lindex $b 2]
    if { [lindex $b 1] != $rows } {
        error "cannot append columns with inequal rows A\[$rows,\] + B\[[lindex $b 1],\]"
        }
    set total_cols [expr {$acols + $bcols}]
    set result [list 2 $rows $total_cols]
    for {set row 0} {$row < $rows} {incr row} {
        for {set col 0} {$col < $acols} {incr col} {
            set i_a [expr {3 + $row * $acols + $col}]
            lappend result [lindex $a $i_a]
            }
        for {set col 0} {$col < $bcols} {incr col} {
            set i_b [expr {3 + $row * $bcols + $col}]
            lappend result [lindex $b $i_b]
            }
        }
    set c $result
    return $c_out
    }
    
########################### join_rows ################################
########################### join_rows ################################
########################### join_rows ################################
#
# combine vector or matrix args as added rows
#
proc join_rows {a b} {
    set dimA [dim_br a]
    set dimB [dim_br b]
    if { $dimA > 2 || $dimB > 2 } { 
        error "cannot handle > 2D data"
        }
    if { $dimA < 2 } {promote_br a}
    if { $dimB < 2 } {promote_br b}
    join_rows_br a b
    return $a
    }    

    
proc join_rows_br {a_in b_in {c_out {}}} {
    upvar $a_in a
    upvar $b_in b
    if { $c_out == {}} {set c_out $a_in }
    upvar $c_out c
    # have to match col counts
    set arows [lindex $a 1]
    set brows [lindex $b 1]
    set acols [lindex $a 2]
    set bcols [lindex $b 2]
    if { $acols != $bcols } {
        error "cannot append rows with inequal columns A\[,$acols\] + B\[,$bcols\]"
        }
    set total_rows [expr {$arows + $brows}]
    set c [concat 2 $total_rows $acols [lrange $a 3 end] [lrange $b 3 end]]
    return $c_out
    }
    
########################### lassign ##################################
########################### lassign ##################################
########################### lassign ##################################

#
# replace a single value in a Tcl list 
# This call assumes the caller knows the desired position in the list.
#
proc lassign {x index value} {
    return [lreplace $x $index $index $value]
    }

# the mysterious K procedure improves efficiency because of how list
# manipulation is performed (thank you Kevin Kenny)
proc K {x y} {set x}

# prototype list element replace by value - this can be done in C code
# without copying the list so it can be efficient, especially since Tcl
# will maintain an internal representation of a list already
# split into elements.  
# lreplace forces you to copy the whole list.
if { [info commands ${La_namespace}::lassign_br] == {}} {
proc lassign_br {listname index value} {
    upvar $listname thelist
    set thelist [lreplace [K $thelist [set thelist {}]] $index $index $value]
    return $listname
    }
  }


################ madd,mdiv, mprod, msub, mat_binary_op ##############
################ madd,mdiv, mprod, msub, mat_binary_op ##############
################ madd,mdiv, mprod, msub, mat_binary_op ##############
#
# vector/matrix binary operations on elements like addition
#
proc mat_binary_op {a b op} {
    set dimA [dim_br a]
    set dimB [dim_br b]
    if { $dimA > 2 || $dimB > 2 } { 
        error "cannot handle > 2D data"
        }
    if { $dimA < 1 || $dimA < $dimB } {promote_br a}
    if { $dimB < 1 || $dimB < $dimA } {promote_br b}
    mat_binary_op_br a b $op c
    return $c
    }

# matrix addition and similar
proc madd {a b} { return [mat_binary_op $a $b +] }
proc mdiv {a b} { return [mat_binary_op $a $b {*1.0/}] }
proc mprod {a b} { return [mat_binary_op $a $b *] }
proc msub {a b} { return [mat_binary_op $a $b -] }


# perform a binary operation using corresponding elements of
# two matrices or vectors - op is commonly +, -, *, / 
# but you can also add formula constructions as in 
# " *2.0 + 3.4*"
# Be careful with division since Tcl tries to perform integer 
# division.  You may want to use "*1.0 / " to insure floating point math
#
proc mat_binary_op_br {a_in b_in op {c_out {}}} {
    upvar $a_in a
    upvar $b_in b
    if { $c_out == {}} { set c_out $a_in }
    upvar $c_out c
    set NR [lindex $a 1]
    set MC [lindex $a 2]
    set NR2 [lindex $b 1]
    set MC2 [lindex $b 2]
    if { $NR != $NR2 || $MC != $MC2 } { 
        error "arguments are not conformable A\[$NR,$MC\] vs B\[$NR2,$MC2\]"
        }
    set result [list 2 $NR $MC]
    set imax [llength $a]
    for {set i 3} {$i < $imax} {incr i} {
        set ai [lindex $a $i]
        set bi [lindex $b $i]
        lappend result [expr $ai $op $bi]
        }
    set c $result
    return $c_out
    }


########################## mathprec ################################
########################## mathprec ################################
########################## mathprec ################################
#
# test math precision

# what is the smallest number such that 
#   1+epsilon > 1
# does the machine truncate or round-off
# what is the radix and number of digits
#
# Intel PC (IEEE 8 byte floating point) results:
# radix=2.0 digits=53 epsilon=2.22044604925e-016 method=truncation
#
proc mathprec {{puts puts}} {
    for {set test 1.0} {$test + 1.0 != $test} {set test [expr {$test*2.0}]} { }
    for {set diff 1.0} {$test + $diff == $test} {set diff [expr {$diff + 1.0}]} { }
    set radix [expr {($test + $diff) - $test}]
    if { $diff < $radix } { set method rounding } else {set method truncation }
    set digits 0
    for {set test 1.0} {$test + 1 != $test} {set test [expr {$test*$radix}]} {
       incr digits
       }
    set epsilon [expr {pow($radix,(1.0-$digits))}]
    if {$puts != {}} {
        $puts "radix=$radix digits=$digits epsilon=$epsilon method=$method"
        }
    return $epsilon
    }
    
############################ mcols #################################
############################ mcols #################################
############################ mcols #################################
#
# you can use mcols and mrows to determine the number of columns or
# rows and keep your code isolated from the details of the 
# data representation
#
proc mcols m { return [lindex $m 2] }
proc mcols_br name_m { upvar $name_m m ; return [lindex $m 2] }
proc mrows m { return [lindex $m 1] } 
proc mrows_br name_m { upvar $name_m m ; return [lindex $m 1] }

############################ mdiag #################################
############################ mdiag #################################
############################ mdiag #################################
#
# create a diagonal matrix from a vector, an Nx1 matrix, or a 1XN matrix
# 
proc mdiag v {
   set dim [dim_br v] 
   if { $dim == 2 } { demote_br v }
   if { [lindex $v 2] != 0 } { error "expected vector argument" }
   set N [lindex $v 1]
   set result [list 2 $N $N]
   for {set row 0} {$row < $N} {incr row} {
       for {set col 0} {$col < $N} {incr col} {
           set iv [expr {3 + $row}]
           if { $row == $col } { lappend result [lindex $v $iv] }\
           else { lappend result 0 }
           }
       }
   return $result
   }

############################ mdingdong ##############################
############################ mdingdong ##############################
############################ mdingdong ##############################
#
# create the Ding Dong test matrix, a Cauchy matrix that is 
# represented inexactly in the machine, but very stable for
# inversion by elimination methods
#
proc mdingdong N {
   if { ![string is integer -strict $N] } { error "improper size \"$N\"" }
   set result [list 2 $N $N]
   for {set row 0} {$row < $N} {incr row} {
       for {set col 0} {$col < $N} {incr col} {
           lappend result [expr {0.5/($N - $col - $row - 0.5)}]
           }
       }
   return $result
   }

############################ mevsvd #################################
############################ mevsvd #################################
############################ mevsvd #################################
#
# eigenvectors/eigenvalues of a real symmetric matrix by
# singular value decomposition
#
proc mevsvd {A {eps 2.3e-16}} {
    puts "here"
    mevsvd_br A evals
puts "eigenvalues=[show $evals]"
    return $A
    }


# the eigenvectors are returned in place as the columns of A
# the eigenvalues are returned as a vector 

proc mevsvd_br {A_in_out evals_out {eps 2.3e-16}} {
    upvar $A_in_out A
    upvar $evals_out evals
    set n [lindex $A 1]
    if { [lindex $A 2] != $n } { error "expecting square matrix" }
    for {set i 0} {$i < $n} {incr i} {
        set ii [expr {3 + $i*$n + $i}]
        set v [lindex $A $ii]
        for {set j 0} {$j < $n} {incr j} {
            if { $i != $j } {
                set ij [expr {3 + $i*$n + $j}]
                set Aij [lindex $A $ij]
                set v [expr {$v - abs($Aij)}]
                }
             }
        if { ![info exists h] } { set h $v }\
        elseif { $v < $h } { set h $v }
        }
    # h is lower bound on Gershgorin region
    if { $h <= $eps } {
        set h [expr {$h - sqrt($eps)}]
        # try to make smallest eigenvalue positive and not too small
        for {set i 0} {$i < $n} {incr i} {
            set ii [expr {3 + $i*$n + $i}]
            set Aii [lindex $A $ii]
            lassign_br A $ii [expr {$Aii - $h}]
            }
        }\
    else {
        set h 0.0
        }

    set count 0
  # top of the iteration
  for {set isweep 0} {$isweep < 30 && $count < $n*($n-1)/2} {incr isweep} {
    set count 0   ;# count of rotations in a sweep

    for {set j 0} {$j < [expr {$n-1}]} {incr j} {
        for {set k [expr {$j+1}]} {$k < $n} {incr k} {
            set p [set q [set r 0.0]]
            for {set i 0} {$i < $n} {incr i} {
                set ij [expr {3+$i*$n+$j}]
                set ik [expr {3+$i*$n+$k}]
                set Aij [lindex $A $ij]
                set Aik [lindex $A $ik]
                set p [expr {$p + $Aij*$Aik}]
                set q [expr {$q + $Aij*$Aij}]
                set r [expr {$r + $Aik*$Aik}]
                }
             if { 1.0 >= 1.0 + abs($p/sqrt($q*$r)) } {
                 if { $q >= $r } {
                     incr count
                     # no rotation needed
                     continue
                     }
                 }
             set q [expr {$q-$r}]
             set v [expr {sqrt(4.0*$p*$p + $q*$q)}]
             if { $v == 0.0 } continue
             if { $q >= 0.0 } {
                 set c [expr {sqrt(($v+$q)/(2.0*$v))}]
                 set s [expr {$p/($v*$c)}]
                 }\
             else {
                 set s [expr {sqrt(($v-$q)/(2.0*$v))}]
                 if { $p < 0.0 } { set s [expr {0.0-$s}] }
                 set c [expr {$p/($v*$s)}]
                 }
             # rotation
             for {set i 0} {$i < $n} {incr i} {
                set ij [expr {3+$i*$n+$j}]
                set ik [expr {3+$i*$n+$k}]
                set Aij [lindex $A $ij]
                set Aik [lindex $A $ik]  
                lassign_br A $ij [expr {$Aij*$c + $Aik*$s}]  
                lassign_br A $ik [expr {$Aik*$c - $Aij*$s}]  
                }
            } ;# k
        } ;# j
    #puts "pass=$isweep skipped rotations=$count"
    } ;# isweep
    
    # now columns are orthogonal, rescale 
    # and flip signs if all negative or zero 
    set evals [list 2 $n 0]
    for {set j 0} {$j < $n} {incr j} {
        set s 0.0
        set notpositive 0
        for {set i 0} {$i < $n} {incr i} {
            set ij [expr {3+$i*$n+$j}]
            set Aij [lindex $A $ij]
            if { $Aij <= 0.0 } { incr notpositive }
            set s [expr {$s + $Aij*$Aij}]
            }
        set s [expr {sqrt($s)}]
        if { $notpositive == $n } { set sf [expr {0.0-$s}] }\
        else { set sf $s }
        for {set i 0} {$i < $n} {incr i} {
            set ij [expr {3+$i*$n+$j}]
            set Aij [lindex $A $ij]
            lassign_br A $ij [expr {$Aij/$sf}]
            }
        lappend evals [expr {$s+$h}]
        }
     return $A_in_out
     }            

############################ mhilbert ###############################
############################ mhilbert ###############################
############################ mhilbert ###############################
#
# create the Hilbert test matrix which is notorious for being 
# ill conditioned for eigenvector/eigenvalue solutions
#
proc mhilbert N {
   if { ![string is integer -strict $N] } { error "improper size \"$N\"" }
   set result [list 2 $N $N]
   for {set row 0} {$row < $N} {incr row} {
       for {set col 0} {$col < $N} {incr col} {
           lappend result [expr {1.0/($col + $row +1.0)}]
           }
       }
   return $result
   }

############################ mident #################################
############################ mident #################################
############################ mident #################################
#
# create an identity matrix of order N
# 
proc mident N {
   if { ![string is integer -strict $N] } { error "improper size \"$N\"" }
   set result [list 2 $N $N]
   for {set row 0} {$row < $N} {incr row} {
       for {set col 0} {$col < $N} {incr col} {
           if { $row == $col } { lappend result 1 }\
           else { lappend result 0 }
           }
       }
   return $result
   }

############################### mlssvd #########################
############################### mlssvd #########################
############################### mlssvd #########################
#
# solve the over-determined linear equations in the least squares
# sense using SVD
#
#  Ax ~ y
#   y[m] is dependent variable such as a measured outcome
#   x[n] vector of independent variables 
#   A[m,n] - each row is a set dependent variable values,
#       the first column is usually 1 to allow for a constant
#       in the regression 
#
# q is the minimum singular value, lessor values are treated as 0
#
proc mlssvd {A y {q 0.0} {puts puts} {epsilon 2.3e-16}} {
    # A[m,n]
    set m [lindex $A 1]
    set n [lindex $A 2]
    promote_br y  ;# now expect y[m,1]
    if { [lindex $y 1] != $m } {
        error "cannot conform A\[$m,$n\]*X\[$n] = y\[[lindex $y 1],[lindex $y 2]]"
        }
    msvd_br A S V  
    # now A has been transformed to U[m,n]
    # S[n], V[n,n]
    if { $puts != {}} {
        $puts singular\ values=[show $S %.6g]
        }
    set tol [expr {$epsilon * $epsilon * $n * $n}]
    # form Utrans*y into g
    set g [list 2 $n 0]
    for {set j 0} {$j < $n} {incr j} {
        set s 0.0
        for {set i 0} {$i < $m} {incr i} {
            set ij [expr {3 + $i*$n + $j}]
            set Aij [lindex $A $ij]
            set yi [lindex $y [expr {3 + $i}]]
            set s [expr {$s + $Aij*$yi}]     
            }
        lappend g $s ;# g[j] = $s
        }
    # form VS+g = VS+Utrans*g
    set x [list 2 $n 0]
    for {set j 0} {$j < $n} {incr j} {
        set s 0.0
        for {set i 0} {$i < $n} {incr i} {
            set iindex [expr {$i+3}]
            set zi [lindex $S $iindex]
            if { $zi > $q } {
                set ji [expr {3 + $j*$n+$i}]
                set Vji [lindex $V $ji]
                set gi [lindex $g $iindex]
                set s [expr {$s + $Vji*$gi/$zi}]
                }
            }
        lappend x $s
        }               
    return $x                
    }

############################ mmult ##################################
############################ mmult ##################################
############################ mmult ##################################
#
# matrix multiplication
#
#  A[p,q] x B[q,r] = C[p,r]
#   Cij = dot product A[i,] row x B[,j] col
#   Cij = sum   Aik * Bkj  ; k=0..q-1
#
#  rules: Vector arguments are always promoted to Nx1 arrays.
#         So chances are if you are using one as a left operand
#         you probably intend to use the transpose of it (1xN),
#         which is easily done using the transpose procedure.
#    
proc mmult {a b} {
    set dimA [dim_br a]
    set dimB [dim_br b]
    if { $dimA > 2 || $dimB > 2 } { 
        error "cannot handle > 2D data"
        }
    if { $dimA < 2 } {promote_br a}
    if { $dimB < 2 } {promote_br b}
    mmult_br a b c
    return $c
    }

# caller is responsible to provide 2D arguments
proc mmult_br {a_in b_in c_out} {
    upvar $a_in a
    upvar $b_in b
    upvar $c_out c
    set p [lindex $a 1]
    set q [lindex $a 2]
    set q2 [lindex $b 1]
    set r [lindex $b 2]
    if { $q2 != $q } { 
        error "matrices are not conformable A\[$p,$q\] x B\[$q2,$r\]"
        }
    # Cij = dot product A[i,] row x B[,j] col
    set c [list 2 $p $r]
    set a_base 3
    set ainc 1   
    set binc $r
    for {set row 0} {$row < $p} {incr row} {
        set b_base 3
        for {set col 0} {$col < $r} {incr col} {
            lappend c [dotprod_br a b $q $a_base $b_base $ainc $binc]
            incr b_base 1
            }
        incr a_base $q
        }

    return $c_out
    }

######################## mnorms, mnormalize ##########################
######################## mnorms, mnormalize ##########################
######################## mnorms, mnormalize ##########################
#
# compute the means and sigmas of each column
#
#
proc mnorms a {
    mnorms_br a means sigmas
    return [transpose [join_cols $means $sigmas]]
    }

proc mnorms_br {a_in means_out sigmas_out} {
    upvar $a_in a
    if { [dim_br a] != 2 } { error "expecting matrix" }
    set NR [lindex $a 1]
    if { $NR < 2 } { error "insufficient rows to calculate sigmas" }
    set NC [lindex $a 2]
    # vector results are returned
    set means [set sigmas [list 2 $NC 0]]
    for {set i 0} {$i < $NC} {incr i} {
        mrange_br a colvect $i $i
        # mean, stddev are Hume C-code Tcl commands
        set coldata [lrange $colvect 3 end]
        lappend means [mean $coldata]
        lappend sigmas [stddev $coldata]
        }
    upvar $means_out m
    set m $means
    upvar $sigmas_out s
    set s $sigmas
    return [list $means_out $sigmas_out]
    }

#
# normalize each column by subtracting the corresponding mean and then
# dividing by the corresponding sigma
#
    
proc mnormalize {a means sigmas} {
   if { [dim_br a] == 1 } { promote_br a }
   mnormalize_br a means sigmas
   return $a
   }

#
# the code will work with means and sigmas as vectors, Nx1, or 1xN matrices
#
proc mnormalize_br {a_in means_in sigmas_in {c_out {}}} {
    upvar $a_in a
    upvar $means_in means
    upvar $sigmas_in sigmas
    if { [dim_br a] != 2 } { error "expecting matrix" }
    set NR [lindex $a 1]
    set NC [lindex $a 2]
    set nm [llength $means] 
    set ns [llength $sigmas]
    if { $nm != $ns || $nm != [expr 3+$NC] } {
        error "non-conformable arguments a[$NR,$NC] means[expr {$nm-3}] sigmas[expr {$ns-3}]"
        }
    set result [list 2 $NR $NC]
    for {set row 0} {$row < $NR} {incr row} {
        for {set i 0} {$i < $NC} {incr i} {
            set vindex [expr {3 + $i}]
            set mean [lindex $means $vindex]
            set sigma [lindex $sigmas $vindex]
            set aindex [expr {3 + $row*$NC + $i}]
            set val [lindex $a $aindex]
            lappend result [expr { (0.0 + $val - $mean)/$sigma }]
            }
        }
    if { $c_out == {}} { set c_out $a_in }
    upvar $c_out c
    set c $result
    return $c_out
   }

# the Hume dmh_wish has these in C code, here they are for other shells
if { [info commands mean] != "mean" } {
proc ::mean numlist {
    set len [llength $numlist]
    if { $len == 0 } { return {}}
    if { $len == 1 } { return [expr {double($numlist)}] }
    set s 0.0
    for {set i 0} {$i < $len} {incr i} {
        set x [lindex $numlist $i]
        set s [expr {$s + $x}]
        }
    return [expr {$s/$len}]
    }
}
if { [info commands stddev] != "stddev" } {
proc ::stddev numlist {
    set len [llength $numlist]
    if { $len < 2 } { return {} }
    set sx [set sxx 0.0]
    for {set i 0} {$i < $len} {incr i} {
        set x [lindex $numlist $i]
        set sx [expr {$sx + $x}]
        set sxx [expr {$sxx + $x * $x}]
        }
    # the real world is full of surprises like stdev {2.41 2.41 2.41}
    # causing a domain error, so you get strange code
    set diff [expr {$sxx - $sx*$sx/$len}]
    if { $diff < 0.0 || ($diff + 0.0 != $diff) } {return 0.0 }
    return [expr {sqrt($diff/($len-1.0))}]
    }
}

#################### mrange ##########################################
#################### mrange ##########################################
#################### mrange ##########################################
#
#
#  return a subset of selected columns, selected rows
#  indexing always starts with 0.  "end" can be used as an index.
#  You are specifying reverse ordering for the selection when start>last.
#
proc mrange {m colstart collast {rowstart 0} {rowlast end}} {
    set dim [dim_br m]
    if { $dim < 2 } { promote_br m }
    mrange_br m m $colstart $collast $rowstart $rowlast
    return $m
    }

#
proc mrange_br {a_in c_out colstart collast {rowstart 0} {rowlast end}} {
   upvar $a_in a
   # matrix is assumed
   set NR [lindex $a 1]
   set NC [lindex $a 2]
   foreach var {colstart collast} {
       if { [set $var] == {end} } { set $var [expr {$NC - 1}] ; continue }
       set $var [expr int([set $var])]
       if { [set $var] < 0 } { error "column index should not be -ve"}\
       elseif { [set $var] >= $NC } { error "column index [set $var] > [expr {$NC - 1}]" }
       }
   if { $colstart <= $collast } { set colstep 1 } else {set colstep -1}
   foreach var {rowstart rowlast} { 
       if { [set $var] == {end} } { set $var [expr {$NR - 1}] ; continue }
       set $var [expr int([set $var])]
       if { [set $var] < 0 } { error "row index should not be -ve"}\
       elseif { [set $var] >= $NR } { error "row index [set $var] > [expr {$NR - 1}]" }
       }
   if { $rowstart <= $rowlast } { set rowstep 1 } else {set rowstep -1}
   set newNR [expr {$rowstep*($rowlast - $rowstart) + 1}]
   set newNC [expr {$colstep*($collast - $colstart) + 1}]
   set result [list 2 $newNR $newNC]
   for {set row $rowstart} {1} {incr row $rowstep} {
       for {set col $colstart} {1} {incr col $colstep} {
           set index [expr {3 + $row*$NC + $col}]
           lappend result [lindex $a $index]
           if { $col == $collast } break
           }
       if { $row == $rowlast } break
       }
   upvar $c_out c
   set c $result
   return $c_out
   }

############################ msolve #################################
############################ msolve #################################
############################ msolve #################################

#
# solve the matrix problem Ax = p for x, where p may be multiple columns
#
# when p is the identity matrix, the solution x, is the inverse of A
#


proc msolve {A p} {
   promote_br p    ;# a vector becomes Nx1
 
   set N [lindex $A 2]
   # paste the rhs columns on the end
   join_cols_br A p

   msolve_br A 

   # now peel off the solution which replaced the rhs columns
   mrange_br A x $N end
   return $x
   }

#
# Gauss elimination with partial pivoting
#
proc msolve_br {A_in {tolerance 2.3e-16}} {
    upvar $A_in A
    set n [lindex $A 1]
    set NC [lindex $A 2]
    set p [expr {$NC - $n}]
    set Det 1.0
    for {set j 0} {$j <= [expr {$n - 2}]} {incr j} {
        set j_plus1 [expr {$j + 1}]
        set jj [expr {3 + $j * $NC + $j}]
        set Ajj [lindex $A $jj]
        set s [expr {abs($Ajj)}]
        set k $j
        for {set h $j_plus1} {$h < $n} {incr h} {
           set Ahj [lindex $A [expr {3+ $h * $NC + $j}]]
           set Ahj [expr {abs($Ahj)}]
           if { $Ahj > $s } {
               set s $Ahj
               set k $h
               }
           } ;# end pivot search
        if { $k != $j } { ;# row interchange
            for {set i $j} {$i < $NC} {incr i} {
                set ki [expr {3+$k*$NC + $i}]
                set ji [expr {3+$j*$NC + $i}]
                set Aki [lindex $A $ki]
                set Aji [lindex $A $ji]
                lassign_br A $ki $Aji
                lassign_br A $ji $Aki
                }
            set Det [expr {0.0-$Det}]
            }
        set Ajj [lindex $A $jj]
        set Det [expr {$Det * $Ajj}]
        if { abs($Ajj) <= $tolerance } { 
            error "matrix is computationally singular at a tolerance of $tolerance"
            }
        for {set k $j_plus1} {$k < $n} {incr k} {
            set kj [expr {3+$k*$NC + $j}]
            set Akj [lindex $A $kj]
            set Akj [expr {double($Akj)/$Ajj}]
            lassign_br A $kj $Akj ;# to form multiplier mkj
            for {set i $j_plus1} {$i < $NC} {incr i} {
                set ki [expr {3+$k*$NC + $i}]
                set Aki [lindex $A $ki]
                set ji [expr {3+$j*$NC + $i}]
                set Aji [lindex $A $ji]
                lassign_br A $ki [expr {$Aki - double($Akj)*$Aji}]
                }
            } ;# k
        } ;# j
    set nn [expr {3+($n-1)*$NC+($n-1)}]  ;# n-1,n-1
    set Ann [lindex $A $nn]
    set Det [expr {$Det * $Ann}]
    if { abs($Ann) < $tolerance } {
        error "matrix is computationally singular at a tolerance of $tolerance"
        }
    # factoring completed
    #    mij stored in i>j
    #    Rij stored in i<=j < n
    #    fij stored in j= n ... n+p-1   
    # Back substitution for Rx=f
    for {set i $n} {$i < $NC} {incr i} {
        set ni [expr {3+($n-1)*$NC + $i}]  ;# really n-1,i
        set Ani [lindex $A $ni]
        set Ani [expr {double($Ani)/$Ann}]
        lassign_br A $ni $Ani
        for {set j [expr {$n-2}]} {$j >= 0} {incr j -1} {
            set ji [expr {3+$j*$NC+$i}]
            set s [lindex $A $ji]
            for {set k [expr {$j+1}]} {$k < $n} {incr k} {
                set jk [expr {3+$j*$NC+$k}]
                set Ajk [lindex $A $jk]
                set ki [expr {3+$k*$NC+$i}]
                set Aki [lindex $A $ki]
                set s [expr {$s - $Ajk*$Aki}]
                }
            set jj [expr {3+$j*$NC+$j}]
            set Ajj [lindex $A $jj]
            lassign_br A $ji [expr {double($s)/$Ajj}]
            }
        }
    return $A_in
    }  
             
 
########################## msum #####################################
########################## msum #####################################
########################## msum #####################################
#
# compute the sums of each column, return a vector or scalar
# call twice to get sum of columns and rows (set total [msum [msum $a]])
#
# a is a vector or matrix
# result is a vector or scalar
proc msum a {
    if { [dim_br a] == 1 } { promote_br a }
    msum_br a sums
    # convert vector to scalar if only 1 column
    return [demote $sums] 
    }

# a_in is matrix
# sums_out writes a vector
proc msum_br {a_in sums_out} {
    upvar $a_in a
    if { [dim_br a] != 2 } { error "expecting matrix" }
    set NR [lindex $a 1]
    set NC [lindex $a 2]
    # vector results are returned
    set sums [list 2 $NC 0]
    for {set j 0} {$j < $NC} {incr j} {
        set sum 0.0
        set index [expr {3 + $j}]
        for {set i 0} {$i < $NR} {incr i} {
            set Aij [lindex $a $index]
            set sum [expr {$sum + $Aij}]
            incr index $NC
            }
        lappend sums $sum
        }
    upvar $sums_out m
    set m $sums
    return $sums_out
    }


################################ msvd ##########################
################################ msvd ##########################
################################ msvd ##########################
#
# Singular Value Decomposition
# decompose matrix A into (U)(S)(Vtrans) where
#     A[m,n] is the original matrix
#     U[m,n] has orthogonal columns (Ut)(U) = (1(k)
#       and multiplies to an identity matrix     ...
#       supplemented with zeroes if needed        0(n-k))
#     V[n,n] is orthogonal   (V)(Vtran) = [mident $n]
#         the eigenvectors representing the principal components
#     S is diagonal with the positive singular values of A
#
# when you have V,S,U
#    A[m,n]V[n,n] = B[m,n]  is a transformation of A to orthogonal columns, B
#    B[m,n] = U[m,n]S[n,n]  
#    square S and divide by (m-1) to get PCA eigenvalues
#
proc msvd A {
    msvd_br A S V
puts U:\n[show $A %12.4f]  
puts \nS:\n[show $S %12.4f]
puts \nV:\n[show $V %12.4f]
    return $V
    }
#
proc msvd_br {A_in_U_out S_out V_out {epsilon 2.3e-16}} {
    upvar $A_in_U_out A
    set m [lindex $A 1]
    set n [lindex $A 2]
    set tolerance [expr {$epsilon * $epsilon* $m * $n}]
    upvar $V_out V
    set V [mident $n]
    upvar $S_out z

  # top of the iteration
  set count 1
  for {set isweep 0} {$isweep < 30 && $count > 0} {incr isweep} {
    set count [expr {$n*($n-1)/2}] ;# count of rotations in a sweep
    for {set j 0} {$j < [expr {$n-1}]} {incr j} {
        for {set k [expr {$j+1}]} {$k < $n} {incr k} {
            set p [set q [set r 0.0]]
            for {set i 0} {$i < $m} {incr i} {
                set ij [expr {3+$i*$n+$j}]
                set ik [expr {3+$i*$n+$k}]
                set Aij [lindex $A $ij]
                set Aik [lindex $A $ik]
                set p [expr {$p + $Aij*$Aik}]
                set q [expr {$q + $Aij*$Aij}]
                set r [expr {$r + $Aik*$Aik}]
                }
            if { $q < $r } {
                set c 0.0
                set s 1.0
                }\
            elseif { $q * $r == 0.0 } { ;# underflow of small elements
                incr count -1
                continue
                }\
            elseif { ($p*$p)/($q*$r) < $tolerance } { ;# cols j,k are orthogonal
                incr count -1
                continue
                }\
            else {
                set q [expr {$q-$r}]
                set v [expr {sqrt(4.0*$p*$p + $q*$q)}]
                set c [expr {sqrt(($v+$q)/(2.0*$v))}]
                set s [expr {$p/($v*$c)}]
                # s == sine of rotation angle, c == cosine
                }
            # rotation of A
            for {set i 0} {$i < $m} {incr i} {
                set ij [expr {3+$i*$n+$j}]
                set ik [expr {3+$i*$n+$k}]
                set Aij [lindex $A $ij]
                set Aik [lindex $A $ik]  
                lassign_br A $ij [expr {$Aij*$c + $Aik*$s}]  
                lassign_br A $ik [expr {$Aik*$c - $Aij*$s}]  
                }
            # rotation of V
            for {set i 0} {$i < $n} {incr i} {
                set ij [expr {3+$i*$n+$j}]
                set ik [expr {3+$i*$n+$k}]
                set Vij [lindex $V $ij]
                set Vik [lindex $V $ik]
                lassign_br V $ij [expr {$Vij*$c + $Vik*$s}]  
                lassign_br V $ik [expr {$Vik*$c - $Vij*$s}]  
                }
            } ;# k
        } ;# j
    #puts "pass=$isweep skipped rotations=$count"
    } ;# isweep

    set z [list 2 $n 0]
    for {set j 0} {$j < $n} {incr j} {
        set q 0.0
        for {set i 0} {$i < $m} {incr i} {
            set ij [expr {3+$i*$n+$j}]
            set Aij [lindex $A $ij]
            set q [expr {$q+$Aij*$Aij}]
            }
        set q [expr {sqrt($q)}]
        lappend z $q
        if { $q >= $tolerance } {
            for {set i 0} {$i < $m} {incr i} {
                set ij [expr {3+$i*$n+$j}]
                set Aij [lindex $A $ij]
                lassign_br A $ij [expr {$Aij/$q}]
                }
            }
        } ;# j
    return [list $A_in_U_out $S_out $V_out]
    }


######### mat_unary_op, madjust, moffset, mscale ####################
######### mat_unary_op, madjust, moffset, mscale ####################
######### mat_unary_op, madjust, moffset, mscale ####################
#
# matrix unary operations like scaling
#
proc mat_unary_op {a op} {
    set dimA [dim_br a]
    if { $dimA > 2} { 
        error "cannot handle > 2D data"
        }
    if { $dimA < 1 } { promote_br a }
    mat_unary_op_br a $op c
    return $c
    }

proc moffset {a delta} {return [mat_unary_op $a "expr $delta +"] }
proc mscale {a delta} { return [mat_unary_op $a "expr $delta *"] }
proc madjust {a scale offset} { return [mat_unary_op $a "expr $offset + $scale *"] }
proc mround_ij {eps aij} {  ;# round off a number if close to an integer
    if { abs($aij - int($aij)) <= $eps } { return [expr {double(int($aij))}] }
    return $aij
    }
proc mround {a {eps 1.0e-8}} { return [mat_unary_op $a [list mround_ij $eps]] }

# perform a unary operation on the elements of a vector or matrix
# the value is lappended to your op and the result is executed
# "expr 0.5 *" divides by 2
# you can use a procedure for more complex formulae
#
proc mat_unary_op_br {a_in op {c_out {}}} {
    upvar $a_in a
    if { $c_out == {}} { set c_out $a_in }
    upvar $c_out c
    set NR [lindex $a 1]
    set MC [lindex $a 2]
    set result [list 2 $NR $MC]
    set imax [llength $a]
    for {set i 3} {$i < $imax} {incr i} {
        lappend result [eval [concat $op [lindex $a $i]]]
        }
    set c $result ;# only changes a_in,c if no error
    return $c_out
    }
########################### promote #################################    
########################### promote ################################# 
########################### promote ################################# 
#
# promote a scalar or vector to an array
# vec[N] promoted to Nx1 array
#
proc promote x {
    promote_br x
    return $x
    }

proc promote_br {name_in {name_out {}}} {
    upvar $name_in x
    if { $name_out == {}} {set name_out $name_in }
    upvar $name_out y
    set dim [dim_br x]
    if { $dim == 0 } {
        set y [list 2 1 1 $x]
        return $name_out
        }
    if { $dim == 1 } { ;# 2 n 0 v1 v2 ... 
        set n [lindex $x 1]
        if { $name_out != $name_in } { set y [lreplace $x 2 2 1] }\
        else { lassign_br y 2 1 }
        return $name_out
        }
    if { $name_out != $name_in } {set y $x}
    return $name_out
    }


# demote an Nx1 or 1XN matrix to a vect[n]
# demote a vect[1] to a scalar
# call twice to demote 1x1 matrix to scalar
proc demote x {
    demote_br x
    return $x
    }

proc demote_br {a_in {c_out {}}} {
    upvar $a_in a
    if { $c_out == {}} { set c_out $a_in }
    upvar $c_out c
    if { $c_out != $a_in } { set c $a }
    # demote 1XN or Nx1 to vector
    set dim [dim_br a]
    if { $dim == 0 || $dim > 2 } {return $c_out }
    set nr [lindex $a 1]
    if { $dim == 2 } {
        set mc [lindex $a 2]
        if { $mc == 1 } {
            set c [lreplace $a 2 2 0]
            }\
        elseif { $nr == 1 } {
            set c [lreplace $a 1 2 $mc 0]
            }
        return $c_out
        }
    # demote vector1 to scalar
    if { $dim == 1 } {
       if { $nr == 1 } { 
           set c [lindex $a 3]
           }
       }
    return $c_out
    }

############################ show ###################################
############################ show ###################################
############################ show ###################################

# 
# Return a formatted string representation for an lalf operand
# the format argument is optionally used per-element
# with the Tcl format command.  Examples:
#  %.4f     gives 4 decimal digit fixed point
#  %12.4e   provides 4 decimal digit exponential format with a fixed
#           width  
#
# The col_join argument can be used to separate values with commas
# and tabs, etc.  Similarly the row_join allows you to define the
# row separation string.
#
proc show {x {format {}} {col_join { }} {row_join \n}} {
    show_br x x $format $col_join $row_join
    return $x
    }

proc show_br {name_in {name_out {}} {format {}} {col_join { }} {row_join \n}} {
    upvar $name_in x 
    if { $name_out == {}} { set name_out $name_in }
    upvar $name_out result
    set dim [dim_br x]
    if { $dim == 0 } { 
        if { $format != {} } {set result [format $format $x] }\
        else { set result $x }
        return $name_out}
    if { $dim == 1 } { ;# 2 n 0 v1 v2 ... 
        if { $format == {}} {
            set result [join [lrange $x 3 end] $col_join]
            }\
        else {
            set temp [list]
            foreach item [lrange $x 3 end] {
                lappend temp [format $format $item]
                }
            set result [join $temp $col_join]
            }
        return $name_out
        }
    if { $dim == 2} {
        if { $format != {}} {
            mat_unary_op_br x [list format $format] result
            } \
        else {
            set result $x
            }
        set rows [list]
        set NR [lindex $result 1]
        set MC [lindex $result 2]
        for {set row 0} {$row < $NR} {incr row} {
            set start [expr {3 + $row*$MC}]
            set end [expr {$start + $MC - 1}]
            lappend rows [join [lrange $result $start $end] $col_join]
            }
        set result [join $rows $row_join]
        return $name_out
        }

    error "> 2-D TBD"
    }
     

########################### transpose ###############################
########################### transpose ###############################
########################### transpose ###############################
#
# exchange [i,j] with [j,i] 
#
# a vector is promoted to a 1xN array by transpose (1 row x N col)
#
#% set a {2 2 3 01 02 03 11 12 13}
#% show $a
# 01 02 03
# 11 12 13
#% transpose $a
# 2 3 2 01 11 02 12 03 13
#% show [transpose $a]
# 01 11
# 02 12
# 03 13
#
proc transpose x { ;# call by value
    transpose_br x
    return $x
    }

# transpose
# call by reference - default output is to overwrite input data
#
proc transpose_br {name_in {name_out {}}} {
    upvar $name_in x
    if { $name_out == {}} {set name_out $name_in }
    upvar $name_out y
    set dim [dim_br x]
    if { $dim == 0 } { return }
    if { $dim == 1 } { ;# {2 N 0 ...} vector to 1xN
        set N [lindex $x 1]
        set y [lreplace $x 1 2 1 $N]
        return $name_out
        }
    if { $dim == 2 } {
        
        set NR [lindex $x 1]
        set MC [lindex $x 2]
        set result [list 2 $MC $NR]
        for {set j 0} {$j < $MC} {incr j} {  ;# j == col in source
            for { set i 0 } { $i < $NR } { incr i } { ;# i == row in source
                lappend result [lindex $x [expr {3 +  $i * $MC + $j}]]
                }
            }
        set y $result 
        return $name_out
        }
    error "transpose cannot handle >2D data"
    }

########################### vdiag #############################
########################### vdiag #############################
########################### vdiag #############################
#
# Create a vector from the diagonal elements of a matrix
#
proc vdiag m {
    promote_br m
    vdiag_br m v
    return $v
    }

proc vdiag_br {name_m_in name_v_out} {
    upvar $name_m_in m
    upvar $name_v_out v
    set N [lindex $m 1]
    set C [lindex $m 2]
    if { $C < $N } { set N $C }
    set v [list 2 $N 0]
    for {set i 0} {$i < $N} {incr i} {
        set ij [expr {$i*($C + 1) + 3}]
        lappend v [lindex $m $ij]
        }
   
   return $name_v_out
   }
 

########################### vtrim #############################
########################### vtrim #############################
########################### vtrim #############################

#
# for a vector, just return the actual vector elements
#   trim away the dimension and size data in the front
#
proc vtrim x {
    vtrim_br x
    return $x
    }

proc vtrim_br {name_in {name_out {}}} {
    upvar $name_in x
    if { $name_out == {}} {set name_out $name_in }
    upvar $name_out y
    set y [lrange $x 3 end]
    return $name_out
    }











} ;# end of namespace
