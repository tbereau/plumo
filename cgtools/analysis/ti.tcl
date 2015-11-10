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

proc ::cgtools::analysis::ti::setup_ti { {eqtime 0} {delta ""} } {
  global ::cgtools::analysis::outputdir
  global ::cgtools::analysis::iotype
  global ::cgtools::analysis::suffix
  variable file_ti
  variable tiequitime
  variable lmbd
  variable ::cgtools::lambda_coupling
  variable ::cgtools::softcore_flag
  variable ::cgtools::softcore_delta

  set lmbd $lambda_coupling
  set tiequitime $eqtime

  set file_ti "$outputdir/time_vs_ti_$lmbd"
  if { [file exists $file_ti] } {
    set newfile 0
  } else {
    set newfile 1
  }
  set pipe [open $file_ti $iotype]
  if { $newfile || $iotype == "w"} {
    puts $pipe "\# time <dE/dlambda>_lambda"
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
  variable ::cgtools::systemtemp
  variable ::cgtools::lambda_coupling

  set currenttime [setmd time]
  if { $currenttime >= $tiequitime } {
    # Now change of force field with original Hamiltonian (lambda = 1)
    ::cgtools::forcefield::update_peptide_ff 1.0

    set epot 0.0
    # Sum up only peptide-lipid contributions
    # lipid beads
    for {set lb 0} {$lb < 8} {incr lb} {
      # peptide beads
      for {set pb 8} {$pb < 43} {incr pb} {
        if { [catch {set epair [analyze energy nonbonded $lb $pb]}] } {
          set epair 0.0
        }
        set epot [expr "$epot + $epair"]
      }
    }
    # Now revert back to original force field with lambda_i
    ::cgtools::forcefield::update_peptide_ff $lmbd

    set pipe [open $file_ti a]
    puts $pipe [format "%16.4f \t %9.4f" $currenttime $epot]
    close $pipe
  }
}
