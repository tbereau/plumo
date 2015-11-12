# Author: Tristan Bereau
# April 2013

package require ::mmsg 1.0.0
package provide ::cgtools::forcefield 1.0.0

namespace eval ::cgtools::forcefield {
    proc source_all_ff { } {
        # lipid
        source_lipid_ff
        # peptide
        source_peptide_ff
    }


    proc source_peptide_ff { } {
        variable ::cgtools::forcefield::partlist_per_res3letter
        variable ::cgtools::forcefield::resparttypelist
        variable ::cgtools::forcefield::respartcharmmbeadlist
        variable ::cgtools::moltypelists
        variable ::cgtools::bonded_parms
        variable ::cgtools::nb_interactions

        # Peptide parameters
        source $::cgtools::cgtoolsdir/forcefield/peptide_parameters.tcl
        source $::cgtools::cgtoolsdir/forcefield/peptide_sc_parameters.tcl
        source $::cgtools::cgtoolsdir/forcefield/peptide.tcl
    }

    proc update_peptide_ff { {lambda_i "1.0"} } {
        variable ::cgtools::nb_interactions
        variable ::cgtools::hremd
        variable ::cgtools::softcore_flag
        variable ::cgtools::lambda_coupling
        # Source the peptide forcefield.
        # Optional argument: lambda_i for FEP calculation.
        set nb_interactions ""
        set lambda_coupling $lambda_i
        source $::cgtools::cgtoolsdir/forcefield/peptide.tcl
        ::cgtools::utils::set_nb_interactions $nb_interactions
        source_hbond_ff
        # Optionally, read additional commands from user
        variable ::cgtools::userscript
        if { [info exists userscript] } {
          ::mmsg::send [namespace current] "Reading user script $userscript"
          source $userscript
        }
        integrate 0
    }

    proc source_hbond_ff { } {
        # Is called by ::cgtools::generation::placemol
        variable ::cgtools::hremd
        variable ::cgtools::forcefield::peptideb::lambda_coupling
        source $::cgtools::cgtoolsdir/forcefield/hbond_inter.tcl
    }

    proc source_lipid_ff { } {
        variable bonded_parms
        variable ::cgtools::moltypelists
        variable ::cgtools::bonded_parms
        variable ::cgtools::nb_interactions
        # Lipid force field
        source $::cgtools::cgtoolsdir/forcefield/popc.tcl
        source $::cgtools::cgtoolsdir/forcefield/dppc.tcl
        source $::cgtools::cgtoolsdir/forcefield/dopc.tcl
    }

}
