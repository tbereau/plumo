#
# Perform Thermodynamic integration at lambda_i
# Author: Tristan Bereau
# November 2015
#

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::ti {
  variable ti
  variable file_ti
  namespace export printav_ti
  namespace export setup_ti
  namespace export analyze_ti
  namespace export resetav_ti
}

proc ::cgtools::analysis::ti::resetav_ti { } {
  # Do nothing
}

proc ::cgtools::analysis::ti::printav_ti { } {
  # Do nothing
}

proc ::cgtools::analysis::ti::setup_ti {  dlambda {eqtime 0} {delta ""} } {
  global ::cgtools::analysis::outputdir
  global ::cgtools::analysis::iotype
  global ::cgtools::analysis::suffix
  variable file_ti
  variable tiequitime
  variable lmbd
  variable dlmbd $dlambda
  variable ::cgtools::lambda_coupling
  variable ::cgtools::softcore_flag
  variable ::cgtools::softcore_delta

  set lmbd $lambda_coupling
  set tiequitime $eqtime

  # Check that dlambda is strictly positive
  if {[expr $dlmbd <= 0.0]} {
    ::mmsg::err [namespace current] "dlambda must be strictly positive"
  }
  # Check that lambda +/- dlambda is within [0,1]
  if { [expr $lmbd - 0.5*$dlmbd < 0.0] || [expr $lmbd + 0.5*$dlmbd > 1.0]} {
    ::mmsg::err [namespace current] "lambda +/- dlambda not within \[0;1\]"
  }

  set file_ti "$outputdir/time_vs_ti_$lmbd"
  if { [file exists $file_ti] } {
    set newfile 0
  } else {
    set newfile 1
  }
  set pipe [open $file_ti $iotype]
  if { $newfile || $iotype == "w"} {
    puts $pipe "\# time dE/dlambda"
  }
  close $pipe

  set softcore_flag 1

  if { $delta != "" } {
    set softcore_delta $delta
  }
  # Set force field -- simulate at lambda_i
  ::cgtools::forcefield::update_peptide_ff $lmbd
}

proc ::cgtools::analysis::ti::analyze_ti { } {
  variable ::cgtools::analysis::topology
  global ::cgtools::analysis::outputdir
  global ::cgtools::analysis::iotype
  global ::cgtools::analysis::suffix
  variable file_ti
  variable tiequitime
  variable lmbd
  variable dlmbd
  variable ::cgtools::systemtemp
  variable ::cgtools::lambda_coupling

  set currenttime [setmd time]
  if { $currenttime >= $tiequitime } {
    # Evaluate force field with lambda+dlambda
    ::cgtools::forcefield::update_peptide_ff [expr $lmbd + 0.5*$dlmbd]
    set epotp 0.0
    # Sum up only peptide-lipid contributions
    # lipid beads
    for {set lb 0} {$lb < 8} {incr lb} {
      # peptide beads
      for {set pb 8} {$pb < 43} {incr pb} {
        if { [catch {set epair [analyze energy nonbonded $lb $pb]}] } {
          set epair 0.0
        }
        set epotp [expr "$epotp + $epair"]
      }
    }
    # Evaluate force field with lambda-dlambda
    ::cgtools::forcefield::update_peptide_ff [expr $lmbd - 0.5*$dlmbd]
    set epotm 0.0
    # Sum up only peptide-lipid contributions
    # lipid beads
    for {set lb 0} {$lb < 8} {incr lb} {
      # peptide beads
      for {set pb 8} {$pb < 43} {incr pb} {
        if { [catch {set epair [analyze energy nonbonded $lb $pb]}] } {
          set epair 0.0
        }
        set epotm [expr "$epotm + $epair"]
      }
    }

    # Now revert back to original force field with lambda_i
    ::cgtools::forcefield::update_peptide_ff $lmbd

    set pipe [open $file_ti a]
    puts $pipe [format "%16.4f %9.4f" $currenttime [expr ($epotp-$epotm)/$dlmbd]]
    close $pipe
  }
}
