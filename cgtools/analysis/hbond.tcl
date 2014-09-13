#
# Analyze peptide hydrogen bonds, as well as helix and beta-sheet content.
# Author: Tristan Bereau
# September 2014
#

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::hbond {
    variable dist_partbilayer
    variable dist_partbilayer_i
    variable f_dispb
    namespace export printav_hbond
    namespace export setup_hbond
    namespace export analyze_hbond
    namespace export resetav_hbond
}

proc ::cgtools::analysis::hbond::resetav_hbond { } {
    variable dist_partbilayer
    variable dist_partbilayer_i
    variable ::cgtools::virtual_sites

    set dist_partbilayer ""
    set dist_partbilayer_i 0
    foreach molVir $virtual_sites {
        lappend dist_partbilayer 0.0
    }
}

proc ::cgtools::analysis::hbond::printav_hbond { } {
    global ::cgtools::analysis::time
    global ::cgtools::analysis::outputdir
    variable ::cgtools::virtual_sites
    variable dist_partbilayer
    variable dist_partbilayer_i

    set f_dispb [open "$outputdir/hbond.dat" a]
    puts -nonewline $f_dispb [format "%15.4f " $time]

    set eleCounter 0
    foreach molVir $virtual_sites {
        puts -nonewline $f_dispb [format "%7.4f " [expr [lindex $dist_partbilayer $eleCounter]/$dist_partbilayer_i]]
        incr eleCounter
    }
    puts $f_dispb ""
    close $f_dispb
}

proc ::cgtools::analysis::hbond::setup_hbond { args } {
    global ::cgtools::analysis::outputdir
    global ::cgtools::analysis::iotype
    variable ::cgtools::virtual_sites
    variable f_dispb

    if { [file exists "$outputdir/hbond.dat"] } {
        set newfile 0
    } else {
        set newfile 1
    }
    set f_dispb [open "$outputdir/hbond.dat" $iotype]
    if { $newfile || $iotype == "w"} {
        puts -nonewline $f_dispb "\# Time "
        foreach i $virtual_sites {
            puts -nonewline $f_dispb "distance(particule-bilayer) "
        }
        puts $f_dispb ""
    }
    close $f_dispb

    resetav_hbond
    
}

proc ::cgtools::analysis::hbond::analyze_hbond { } {
    variable ::cgtools::analysis::topology
    global ::cgtools::analysis::outputdir
    global ::cgtools::analysis::iotype
    variable dist_partbilayer
    variable dist_partbilayer_i
    variable ::cgtools::virtual_sites
    variable ::cgtools::partID_membrane_midplane
    global ::cgtools::analysis::time
    variable f_dispb

    if { [llength $virtual_sites] == 0 } { return 0 }

    set memcomz [lindex [part $partID_membrane_midplane print pos] 2]
    # Update fake particle position
    part $partID_membrane_midplane pos 1 1 $memcomz
    set eleCounter 0

    foreach molVir $virtual_sites {
        set posVirSit [part $molVir print pos]
        set posz [lindex $posVirSit 2]
        lset dist_partbilayer $eleCounter [expr [lindex $dist_partbilayer $eleCounter]+abs($posz - $memcomz)]
        # Check that the virtual site is close to the center of mass of the
        # molecule. Otherwise exit (parallelization issue).
        set molid [part $molVir print mol]
        set posCoM "0.0 0.0 0.0"
        set comIdx 0
        for {set j 0} {$j < [setmd n_part]} {incr j} {
            if {[part $j print mol] == $molid} {
                set posCoM [::cgtools::utils::add2vec $posCoM [part $j print pos]]
                incr comIdx
            }
        }
        set posCoM [::cgtools::utils::scalevec $posCoM [expr 1./$comIdx]]
        set disCoMVir [::cgtools::utils::distance $posCoM $posVirSit]
        if { $disCoMVir > 5.0 } {
            ::mmsg::err [namespace current] "Virtual site of molecule $molid away from its center of mass (dis=$disCoMVir)."
        }
        incr eleCounter
    }
    incr dist_partbilayer_i

    return
}

