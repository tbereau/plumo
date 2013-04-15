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

proc ::cgtools::analysis::fep::setup_fep { lambda_i lambda_j {eqtime 0} {delta ""} } {
  global ::cgtools::analysis::outputdir
  global ::cgtools::analysis::iotype
  global ::cgtools::analysis::suffix
  variable file_fep
  variable summedExp
  variable counter
  variable fepequitime
  variable lmbd_i
  variable lmbd_j
  variable ::cgtools::forcefield::peptideb::lambda_coupling
  variable ::cgtools::forcefield::peptideb::softcore_flag
  variable ::cgtools::forcefield::peptideb::softcore_delta

  set lmbd_i $lambda_i
  set lmbd_j $lambda_j
  set summedExp 0.0
  set counter   0
  set fepequitime $eqtime

  set file_fep "$outputdir/time_vs_fep$suffix"
  if { [file exists $file_fep] } {
    set newfile 0
  } else {
    set newfile 1
  }
  set pipe [open $file_fep $iotype]
  if { $newfile || $iotype == "w"} {
    puts $pipe "\# time E(lambda_i) E(lambda_i+1) deltaE <deltaF>"
  }
  close $pipe

  set softcore_flag 1

  if { $delta != "" } {
    set softcore_delta $delta
  }
  # Set force field -- simulate at lambda_i
  ::cgtools::forcefield::update_peptide_ff $lmbd_i
}

proc ::cgtools::analysis::fep::analyze_fep { } {
  variable ::cgtools::analysis::topology
  global ::cgtools::analysis::outputdir
  global ::cgtools::analysis::iotype
  global ::cgtools::analysis::suffix
  variable file_fep
  variable summedExp
  variable counter
  variable fepequitime
  variable lmbd_i
  variable lmbd_j
  variable ::cgtools::systemtemp
  variable ::cgtools::forcefield::peptideb::lambda_coupling

  set currenttime [setmd time]
  if { $currenttime >= $fepequitime } {
    set epot_i [expr [analyze energy total] - [analyze energy kinetic]]
    # Now change of force field with other lambda coupling
    ::cgtools::forcefield::update_peptide_ff $lmbd_j

    set epot_j [expr [analyze energy total] - [analyze energy kinetic]]
    # Now revert back to original force field with lambda_i
    ::cgtools::forcefield::update_peptide_ff $lmbd_i
    
    set deltaE [expr $epot_j - $epot_i]
    set summedExp [expr $summedExp + exp(-$deltaE/$systemtemp)]
    incr counter

    set pipe [open $file_fep a]
    puts $pipe [format "%16.4f \t %9.4f \t %9.4f \t %9.4f \t %9.4f" \
      $currenttime $epot_i $epot_j $deltaE \
      [expr -$systemtemp * log($summedExp/$counter)] ]
    close $pipe
  }
}

