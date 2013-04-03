# Author: Tristan Bereau
# April 2013

package require ::mmsg 1.0.0
package provide ::cgtools::forcefield 1.0.0

namespace eval ::cgtools::forcefield {
    proc source_peptide_ff { } {
        # Peptide parameters
        source [file join [file dirname [info script]] peptide_parameters.tcl]
        source [file join [file dirname [info script]] peptide_sc_parameters.tcl]
        source [file join [file dirname [info script]] peptide.tcl]   
    }

    proc update_peptide_ff { lambda_i "1.0" } {
        # Source the peptide forcefield.
        # Optional argument: lambda_i for FEP calculation.
        namespace eval :: {
            set peptideb::nb_interactions ""
            set peptideb::lambda_coupling $lambda_i
            source [file join [file dirname [info script]] peptide.tcl]
            set_nb_interactions $peptideb::nb_interactions
        }
    }

    proc source_hbond_ff { } {
        source [file join [file dirname [info script]] hbond_inter.tcl]
    }

    proc source_popc_ff { } {
        # Lipid force field
        source [file join [file dirname [info script]] popc.tcl]
    }

}

