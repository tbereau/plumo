#
# Perform Free-Energy Perturbation between lambda_i and lambda_j
# on all peptide beads present in the simulation.
# Author: Tristan Bereau
# Spring 2013
#

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::fep {
  variable fep
  variable file_fep
  namespace export printav_fep
  namespace export setup_fep
  namespace export analyze_fep
  namespace export resetav_fep
}

proc ::cgtools::analysis::fep::resetav_fep { } {
  # Do nothing 
}

proc ::cgtools::analysis::fep::printav_fep { } {
  # Do nothing
}

proc ::cgtools::analysis::fep::setup_fep { lambda_i lambda_j args } {
  global ::cgtools::analysis::outputdir
  global ::cgtools::analysis::iotype
  variable file_fep
  variable summedExp
  variable counter
  variable lmbd_i
  variable lmbd_j
  variable ::peptideb::lambda_coupling

  set lmbd_i $lambda_i
  set lmbd_j $lambda_j
  set summedExp 0.0
  set counter   0

  if { [file exists "$outputdir/fep.dat"] } {
    set newfile 0
  } else {
    set newfile 1
  }
  set file_fep [open "$outputdir/fep.dat" $iotype]
  if { $newfile || $iotype == "w"} {
    puts $f_dispb "\# time E(lambda_i) E(lambda_i+1) deltaE <deltaF>"
  }
  close $file_fep

  # Set force field -- simulate at lambda_i
  namespace eval :: {
    set ::peptideb::nb_interactions ""
    set ::peptideb::lambda_coupling $lmbd_i
    source [file join [file dirname [info script]] ../topologies/peptide.tcl]
    set_nb_interactions $::peptideb::nb_interactions
  }
  integrate 0

}

proc ::cgtools::analysis::fep::analyze_fep { } {
  variable ::cgtools::analysis::topology
  global ::cgtools::analysis::outputdir
  global ::cgtools::analysis::iotype
  variable file_fep
  variable summedExp
  variable counter
  variable lmbd_i
  variable lmbd_j
  variable ::cgtools::systemtemp
  variable ::peptideb::lambda_coupling

  set epot_i [expr [analyze energy total] - [analyze energy kinetic]]
  # Now change of force field with other lambda coupling
  namespace eval :: {
    set ::peptideb::nb_interactions ""
    set ::peptideb::lambda_coupling $lmbd_j
    source [file join [file dirname [info script]] ../topologies/peptide.tcl]
    set_nb_interactions $::peptideb::nb_interactions
  }
  integrate 0
  set epot_j [expr [analyze energy total] - [analyze energy kinetic]]
  # Now revert back to original force field with lambda_i
  namespace eval :: {
    set ::peptideb::nb_interactions ""
    set ::peptideb::lambda_coupling $lmbd_i
    source [file join [file dirname [info script]] ../topologies/peptide.tcl]
    set_nb_interactions $::peptideb::nb_interactions
  }
  integrate 0

  set deltaE [expr $epot_j - $epot_i]
  set summedExp [expr $summedExp + exp($deltaE/$systemtemp)]
  incr counter

  set file_fep [open "$outputdir/fep.dat" a]
  puts $file_fep "%16.4f \t %9.4f \t %9.4f \t %9.4f \t %9.4f" \
    [setmd time] $epot_i $epot_j $deltaE \
    [expr -$systemtemp * log($summedExp/$counter)] ]
  close $file_fep

  return $dis_partbilayer
}

