#
# Evaluate Multi-state Bennett acceptance ratio at lambda_i
# Author: Tristan Bereau
# November 2015
#

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::mbar {
  variable mbar
  variable file_mbar
  namespace export printav_mbar
  namespace export setup_mbar
  namespace export analyze_mbar
  namespace export resetav_mbar
}

proc ::cgtools::analysis::mbar::resetav_mbar { } {
  # Do nothing
}

proc ::cgtools::analysis::mbar::printav_mbar { } {
  # Do nothing
}

proc ::cgtools::analysis::mbar::setup_mbar {  lambdas {eqtime 0} {delta ""} } {
  global ::cgtools::analysis::outputdir
  global ::cgtools::analysis::iotype
  global ::cgtools::analysis::suffix
  variable file_mbar
  variable mbarequitime
  variable lmbds $lambdas
  variable ::cgtools::lambda_coupling
  variable ::cgtools::softcore_flag
  variable ::cgtools::softcore_delta

  set mbarequitime $eqtime

  # Check that each dlambda is strictly positive
  foreach li $lmbds {
    if {[expr $li < 0.0] || [expr $li > 1.0]} {
      ::mmsg::err [namespace current] "Each lambda must belong to [0,1]"
    }
  }

  set file_mbar "$outputdir/time_vs_mbar"
  if { [file exists $file_mbar] } {
    set newfile 0
  } else {
    set newfile 1
  }
  set pipe [open $file_mbar $iotype]
  if { $newfile || $iotype == "w"} {
    puts -nonewline $pipe "\# time "
    foreach li $lmbds {
      puts -nonewline $pipe " E_$li"
    }
  }
  puts $pipe ""
  close $pipe

  set softcore_flag 1

  if { $delta != "" } {
    set softcore_delta $delta
  }
}

proc ::cgtools::analysis::mbar::analyze_mbar { } {
  variable ::cgtools::analysis::topology
  global ::cgtools::analysis::outputdir
  global ::cgtools::analysis::iotype
  global ::cgtools::analysis::suffix
  variable file_mbar
  variable mbarequitime
  variable lmbds
  variable ::cgtools::systemtemp
  variable ::cgtools::lambda_coupling

  set currenttime [setmd time]
  if { $currenttime >= $mbarequitime } {
    set pipe [open $file_mbar a]
    puts -nonewline $pipe [format "%16.4f " $currenttime]
    foreach li $lmbds {
      # Evaluate force field with li
      ::cgtools::forcefield::update_peptide_ff $li
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

      puts -nonewline $pipe [format "%9.4f " $epot]
    }
    puts $pipe ""
    close $pipe
    # Now revert back to original force field with lambda_i
    ::cgtools::forcefield::update_peptide_ff $lambda_coupling
  }
}
