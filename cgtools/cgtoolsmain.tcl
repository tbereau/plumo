# Zun-Jing Wang 1st version:  2008 Sep. 23 - Oct. 9, 
#		modification: 2009 Sep.02
# Tristan Bereau: 2013. Introduce FEP.
#------------------------------------------------------------#
# Script for running a CG peptide-membrane simulation
#
# Implements settings contained in a parameter file which should be
# given as an argument.
#
# -----------------------------------------------------------#
# Get the name of the current namespace
set this [namespace current]

# --- ensure required packages are available  ---------------#
set result [package require ::mmsg]
::mmsg::send $this ""
::mmsg::send $this "               ***********************  "
::mmsg::send $this "              *                       * "
::mmsg::send $this "              *          PLUM         * "
::mmsg::send $this "              *     ZUN-JING WANG     * "
::mmsg::send $this "              *     TRISTAN BEREAU    * "
::mmsg::send $this "              *      (2008-2013)      * "
::mmsg::send $this "              *                       * "
::mmsg::send $this "               ***********************  "
::mmsg::send $this ""
::mmsg::send $this "loaded version [package require cgtools] of cgtools"
::mmsg::send $this "loaded version [package require cmdline] of cmdline"
::mmsg::send $this "loaded version [package require ::cgtools::analysis] of analysis"
::mmsg::send $this "loaded version [package require ::cgtools::generation] of generation"
::mmsg::send $this "loaded version [package require ::cgtools::utils] of utils"
::mmsg::send $this "loaded version [package require ::cgtools::espresso] of espresso"
::mmsg::send $this "loaded version [package require ::cgtools::forcefield] of forcefield"
::mmsg::send $this "loaded version [package require ::mmsg] of mmsg."
::mmsg::send $this ""
# ---------------------------------------------------------- #

# ---- Process Command Line Args --------------------------- #
set options { 
    {n.arg      1    set the number of processors }
}

set usage "Usage: 
\tcgtoolsmain.tcl <CONFIG_NAME>
     or 
\tEspresso cgtoolsmain.tcl <CONFIG_NAME> 
     * Possible options :
\t-replica \[-connect HOST\] \t\t starts a replica exchange simulation
\t-hremd \[-connect HOST\] \[-port PORT\] \t\t Hamiltonian Replica Exchange MD
\t-hybrid \t\t starts a MC-MD hybrid simulation
\t-annealing \t\t starts a annealing simulation
\t-new \t\t start a new simulation rather than starting from the last checkpoint (depending on CONFIG_FILE information)\n"
array set params [::cmdline::getoptions argv $options $usage]

if { $argc<1} {
    ::mmsg::err $this "Parameter file missing. $usage"
}

namespace eval ::cgtools {
    variable moltypelists 

    set this [namespace current]
    ::mmsg::setnamespaces $this

    # Set default replica exchange contoller
    set replica 0
    set replica_connect 0
    set tcp_port 12000

    # Set default HREMD exchange controller
    variable hremd 0
    set hremd_connect 0
    # Implicit membrane environment will trigger the softcore flag
    variable implicit_membrane 0
    # Restraints on bilayer
    variable membrane_restraint       0
    variable membrane_restraint_k     0.5
    variable membrane_restraint_dist -1.0
    # resume simulation
    variable pdb_resume ""
    # Multi-timestepping
    variable multitimestep 0
    # fake particle that monitors the membrane midplane
    variable partID_membrane_midplane -1
    # Restraints for umbrella sampling. In the form: {molID {{k z0} {k z0}}}
    # for restraints on the first two particles of molecule ID $molID.
    variable umbrella_restraints ""
    # Store all virtual sites
    variable virtual_sites ""
    # Peptide's temperature
    variable temp_peptide -1.

    # Set default system parameters.
    set thermo Langevin
    set warmup_temp 0.7
    set warmup_freq 1
    set warmsteps 20
    set warmtimes 20
    set free_warmsteps 20
    set free_warmtimes 20
    set startmdtime 0
    set npt "off"
    set mgrid 8
    set stray_cut_off 1000.0
    set use_vmd "offline"
    set tablenames ""
    set overlapnames ""
    set trappedmols ""
    set vmdcommands ""
    set userfixedparts ""
    set moltypelists ""

    # Get directory of cgtools
    set cgtoolsdir ""
    foreach dir $auto_path {
        if { [string first cgtools $dir] > -1 } {
            set cgtoolsdir $dir
            if { [string first "~" $cgtoolsdir] == 0 } {
                set cgtoolsdir "$::env(HOME)"
                append cgtoolsdir [string range $dir 1 end]
            }
        }
    }
    if { $cgtoolsdir == "" } {
        ::mmsg::err "Can't find cgtools directory in \$auto_path variable."
    }

    # Directories now stored in cgtools
    set tabledir "$cgtoolsdir/forcefield/forcetables"
    set overlapdir "$cgtoolsdir/forcefield/overlapcoffs"

    # Read the parameter file. 
    set paramsfile [lindex $argv 0]
    set nprocessors $params(n)
    ::mmsg::send $this "using paramsfile: $paramsfile"
    source $paramsfile

    # Make sure we don't have multitimestep with NPT
    if { $multitimestep > 0 } {
        require_feature MULTI_TIMESTEP
        ::mmsg::send $this "Turning on multi-timestepping algorithm."
        ::mmsg::send $this [format "Main time step: %7.4f" $main_time_step]
        ::mmsg::send $this [format "Smaller time step: %7.4f" [expr $main_time_step/$multitimestep]]
    }

    if { $npt == "on" } {
        require_feature NPT
    }

    # MPI distribution
    # MPI distribution of bilayer
    if { [lindex geometry 1] != "random"} {
        if {[setmd n_nodes]==2} {
            setmd node_grid 2 1 1
        } elseif {[setmd n_nodes]==4} {
            setmd node_grid 2 2 1
        } elseif {[setmd n_nodes]==6} {
            setmd node_grid 3 2 1
        } elseif {[setmd n_nodes]==8} {
            setmd node_grid 4 2 1
        } elseif {[setmd n_nodes]==16} {
            setmd node_grid 4 4 1
        } elseif {[setmd n_nodes]==32} {
            setmd node_grid 8 4 1
        } elseif {[setmd n_nodes]==64} {
            setmd node_grid 8 8 1
        }
    }

    # MPI distribution of random 
    if { [lindex geometry 1] == "random"} {
        if {[setmd n_nodes]==8} {
            setmd node_grid 2 2 2
        } elseif {[setmd n_nodes]==64} {
            setmd node_grid 4 4 4
        }
    }

    # Cmdline arguments
    if { $argc >= 2 } {
        for { set k 1 } { $k < $argc } { incr k } {
            switch -- [lindex $argv $k] "-replica" {
                set replica 1
                ::mmsg::send $this "Replica exchange algorithm has been turned on."
                if {[lindex $argv [expr $k+1]]== "-connect"} {
                    set replica_connect [lindex $argv [expr $k+2]]
                    incr k 2
                }
                if {[lindex $argv [expr $k+1]] == "-port" } {
                    set tcp_port [lindex $argv [expr $k+2]]
                    incr k 2
                } 
            } "-hremd" {
                set hremd 1
                ::mmsg::send $this "Hamiltonian replica exchange MD turned on"
                if {[lindex $argv [expr $k+1]] == "-connect" } {
                    set hremd_connect [lindex $argv [expr $k+2]]
                    incr k 2
                }
                if {[lindex $argv [expr $k+1]] == "-port" } {
                    set tcp_port [lindex $argv [expr $k+2]]
                    incr k 2
                } 
            } "-new" {
                set newcomp 1
                ::mmsg::send $this "Start a new simulation, i.e., will NOT be resumed at the last checkpoint."
                incr k
            } "-hybrid" {
                set hybrid 1
                ::mmsg::send $this "Start a hybrid MC-MD simulation."
                incr k
            } "-annealing" {
                set annealing 1
                ::mmsg::send $this "Start an annealing MD simulation."
                incr k
            }
        }
    } 
    ::mmsg::send $this "Finished parsing options."

    # ---------------------------------------------------------- #
    # Allow children namespaces that we can explicitly allow messages from these
    ::mmsg::send $this "Assigning namespaces"
    if { [ catch { ::mmsg::setnamespaces \
        "{:: [namespace children ::cgtools] [namespace children ::parallel_tempering]}" } errmsg ] } {
        puts "Warning: $errmsg"
    }
    set message_allowlist { :: $this ::cgtools ::cgtools::utils ::cgtools::generation ::cgtools::analysis ::cgtools::espresso}
    ::mmsg::send $this "Assigning children namespaces"
    set children [namespace children ::cgtools::analysis]
    foreach child $children {
        lappend message_allowlist $child
    }
    set children [namespace children ::cgtools::generation]
    foreach child $children {
        lappend message_allowlist $child
    }
    set children [namespace children ::cgtools::utils]
    foreach child $children {
        lappend message_allowlist $child
    }
    # Set the namespaces from which messages will be printed
    if { [ catch { ::mmsg::setnamespaces $message_allowlist } errmsg ] } {
        puts "Warning: $errmsg"
    }

    # Enable debug messages
    #::mmsg::enable debug

    # Read forcefield
    ::mmsg::send $this "Reading force field"
    ::cgtools::forcefield::source_all_ff

    # Erase all the directory we're not resuming
    if { [info exists newcomp] } {
        foreach file [glob -nocomplain -directory $ident $ident.vmd*] {
            catch {exec rm -rf $ident}
        }
    }

    # ----------------- Start the script ----------------------- #
    if { $replica == 1 } {
        # Replica exchange is turned on

        # If <ident> doesn't exist then create it
        catch { exec mkdir -p $ident }

        ::mmsg::send $this "Starting parallel tempering at temperatures: $replica_temps"

        if {$replica_connect == 0} {
            # attempt to run PT 3 times if it catches an error
            set n_tries 0
            set max_n_tries 3
            while { $n_tries <= $max_n_tries } {
                set file_hn [open "hostfile" w]
                puts $file_hn "[exec hostname] $tcp_port"
                close $file_hn
                ::mmsg::send $this "host [exec hostname], port $tcp_port"
                if { [ catch { parallel_tempering::main -values $replica_temps -rounds $replica_rounds \
                    -init cgtools::espresso::replica_init -swap cgtools::espresso::replica_swap \
                    -perform cgtools::espresso::replica_perform -port $tcp_port -info comm } errmsg ] } {
                    ::mmsg::send $this $errmsg
                } else {
                    break
                }
                incr tcp_port
                incr n_tries
            }
        } else {
            # attempt to run PT 3 times if it catches an error
            set n_tries 0
            set max_n_tries 3
            while { $n_tries <= $max_n_tries } {
                ::mmsg::send $this "Connection attempt [expr $n_tries+1]/[expr $max_n_tries]..."
                # Attempt to connect 10 seconds after the hostname file has been written
                set long_enough 0
                while { $long_enough == 0 } {
                    if { [file exists hostfile] } {
                        after 10000
                        set file_hn [open "hostfile" r]
                        set data [read $file_hn]
                        set data [split $data " "]
                        close $file_hn
                        set replica_connect [lindex $data 0]
                        set tcp_port [lindex $data 1]
                        set long_enough 1
                    }
                }
                ::mmsg::send $this "host [exec hostname], port $tcp_port"
                if { [catch { parallel_tempering::main -connect $replica_connect \
                    -init cgtools::espresso::replica_init \
                    -swap cgtools::espresso::replica_swap \
                    -perform cgtools::espresso::replica_perform -port $tcp_port -info comm } errmsg ] } {
                    mmsg::send $this $errmsg
                } else {
                    break
                }
                incr n_tries
            }
        }
    } elseif { $hremd != 0 } {
        # HREMD is turned on
        # If <ident> doesn't exist then create it
        catch { exec mkdir -p $ident }
        ::mmsg::send $this "Starting HREMD simulation with coupling lambda:\n  $lambda_values"
        if {$hremd_connect == 0} {
            # attempt to run HREMD 3 times if it catches an error
            set n_tries 0
            set max_n_tries 3
            while { $n_tries <= $max_n_tries } {
                set file_hn [open "hostfile" w]
                puts $file_hn "[exec hostname] $tcp_port"
                close $file_hn
                ::mmsg::send $this "host [exec hostname], port $tcp_port"
                if { [ catch { parallel_tempering::main -values $lambda_values -rounds $replica_rounds \
                    -init cgtools::espresso::hremd_init -swap cgtools::espresso::hremd_swap \
                    -perform cgtools::espresso::hremd_perform -port $tcp_port -info comm } errmsg ] } {
                    ::mmsg::send $this $errmsg
                } else {
                    break
                }
                incr tcp_port
                incr n_tries
            }
        } else {
            # attempt to run PT 3 times if it catches an error
            set n_tries 0
            set max_n_tries 3
            while { $n_tries <= $max_n_tries } {
                ::mmsg::send $this "Connection attempt [expr $n_tries+1]/[expr $max_n_tries]..."
                # Attempt to connect 10 seconds after the hostname file has been written
                set long_enough 0
                while { $long_enough == 0 } {
                    if { [file exists hostfile] } {
                        after 10000
                        set file_hn [open "hostfile" r]
                        set data [read $file_hn]
                        set data [split $data " "]
                        close $file_hn
                        set replica_connect [lindex $data 0]
                        set tcp_port [lindex $data 1]
                        set long_enough 1
                    }
                }
                ::mmsg::send $this "host [exec hostname], port $tcp_port"
                if { [catch { parallel_tempering::main -connect $hremd_connect \
                    -init cgtools::espresso::hremd_init \
                    -swap cgtools::espresso::hremd_swap \
                    -perform cgtools::espresso::hremd_perform -port $tcp_port -info comm } errmsg ] } {
                    mmsg::send $this $errmsg
                } else {
                    break
                }
                incr n_tries
            }
        }
    } else {
        # Replica exchange / HREMD is off

        if { [info exists hybrid] } {
            # Start hybrid simulation 
            ::mmsg::send $this "Starting an Hybrid computation at temperature : $systemtemp"
            espresso::hybrid_init
        } else {
            if { [info exists annealing] }  {
                # Start Annealing
                ::mmsg::send $this "Starting an Annealing computation at temperature : $systemtemp"
                espresso::annealing_init
            } else {
                # Normal MD simulation 
                ::mmsg::send $this "Starting an Espresso computation at temperature : $systemtemp"
                espresso::espresso_init 	
            }
        }
    }

    # terminate program
    mmsg::send $this "\n\nfinished"
    exit 1
}

