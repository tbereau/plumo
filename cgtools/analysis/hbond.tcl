#
# Analyze peptide hydrogen bonds, as well as helix and beta-sheet content.
# Author: Tristan Bereau
# September 2014
#

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::hbond {
    variable helicity
    variable helicity_i
    variable hbondE
    variable hbondE_i
    variable f_hbond
    variable f_hbond_file
    namespace export printav_hbond
    namespace export setup_hbond
    namespace export analyze_hbond
    namespace export resetav_hbond
}

proc ::cgtools::analysis::hbond::resetav_hbond { } {
    variable helicity
    variable helicity_i
    variable hbondE
    variable hbondE_i

    set helicity 0.0
    set helicity_i 0
    set hbondE 0.0
    set hbondE_i 0

}

proc ::cgtools::analysis::hbond::printav_hbond { } {
    global ::cgtools::analysis::time
    global ::cgtools::analysis::outputdir
    variable helicity
    variable helicity_i
    variable hbondE
    variable hbondE_i
    variable f_hbond_file

    set f_hbond [open $f_hbond_file a]
    
    puts -nonewline $f_hbond [format "%15.4f " $time]
    puts -nonewline $f_hbond [format "%7.4f " [expr $helicity/$helicity_i]]
    puts $f_hbond [format "%7.4f " [expr $hbondE/$hbondE_i]]

    close $f_hbond
}

proc ::cgtools::analysis::hbond::setup_hbond { args } {
    global ::cgtools::analysis::iotype
    variable ::cgtools::virtual_sites
    variable ::cgtools::ragtime_path
    global ::cgtools::analysis::outputdir
    variable f_hbond
    variable f_hbond_file

    set f_hbond_file "$outputdir/time_vs_hbond.dat"

    # Check that the ragtime package can be found.
    if { ![file exists $ragtime_path ] } {
        ::mmsg::err $this "Ragtime package wasn't found.\n
        Should be $ragtime_path"
    }

    if { [file exists $f_hbond_file] } {
        set newfile 0
    } else {
        set newfile 1
    }
    set f_hbond [open $f_hbond_file $iotype]
    if { $newfile || $iotype == "w"} {
        puts $f_hbond "\# Time   Helicity_ratio   H-bond_energy"
    }
    close $f_hbond

    resetav_hbond
    
}

proc ::cgtools::analysis::hbond::analyze_hbond { } {
    variable ::cgtools::analysis::topology
    global ::cgtools::analysis::iotype
    variable helicity
    variable helicity_i
    variable hbondE
    variable hbondE_i
    global ::cgtools::analysis::time
    variable ::cgtools::espresso::pdb_output
    variable ::cgtools::ragtime_path

    # Helicity
    # run ragtime. it's possible we get an error, e.g. there's no hbonds. catch it.
    if {![catch {set output [exec $ragtime_path $pdb_output]} errmsg]} {
        # Helicity percentage
        set aggr_list [lsearch -all $output "percentage"]
        foreach position $aggr_list {
            if {[lindex $output [expr $position-1]]=="Helicity"} {
                set helicity [expr $helicity + \
                    [lindex $output [expr $position+2]]]
            }
        }  
    } else {
        ::mmsg::send [namespace current] "ragtime encountered an error:"
        ::mmsg::send [namespace current] $errmsg
        ::mmsg::send [namespace current] "The parameter will be set *Arbitrarily* to 0."
        set param1 0.
    }                      
    incr helicity_i

    # H-bond energy
    set hbondE [expr $hbondE + [analyze energy nonbonded 8 11]]
    incr hbondE_i

    return
}

