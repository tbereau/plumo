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
::mmsg::send $this "loaded version [package require ::cgtools::utils] of espresso"
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
\t-hybrid \t\t starts a MC-MD hybrid simulation
\t-annealing \t\t starts a annealing simulation
\t-new \t\t start a new simulation rather than starting from the last checkpoint (depending on CONFIG_FILE information)\n"
array set params [::cmdline::getoptions argv $options $usage]

if { $argc<1} {
    ::mmsg::err $this "Parameter file missing. $usage"
}

namespace eval ::cgtools {
    set this [namespace current]
    ::mmsg::setnamespaces $this

    # Set default replica exchange contoller
    set replica 0
    set replica_connect 0

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

    # Directories now stored in cgtools
    set tabledir "./forcefield/forcetables"
    set overlapdir "./forcefield/overlapcoffs"

    # Read the parameter file. 
    set paramsfile [lindex $argv 0]
    set nprocessors $params(n)
    ::mmsg::send $this "using paramsfile: $paramsfile"
    source $paramsfile

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

    # ---------------------------------------------------------- #
    # Allow children namespaces that we can explicitly allow messages from these
    catch {::mmsg::setnamespaces ":: [namespace children ::cgtools] [namespace children ::parallel_tempering]"}
    set message_allowlist { :: ::cgtools::utils ::cgtools::generation ::cgtools::analysis }
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
    catch { ::mmsg::setnamespaces $message_allowlist }


    # Enable debug messages
    #::mmsg::enable debug

    # ----------------- Start the script ----------------------- #
    if { $replica == 1 } {
        # Replica exchange is turned on

        # Erase all the directory we're not resuming
        if { [info exists newcomp] } {
            foreach file [glob -nocomplain -directory $ident $ident.vmd*] {
                catch {exec rm -rf $ident}
            }
        }

        # If <ident> doesn't exist then create it
        catch { exec mkdir $ident }

        ::mmsg::send $this "Starting parallel tempering at temperatures : $replica_temps"

        if {$replica_connect == 0} {
            parallel_tempering::main -values $replica_temps -rounds $replica_rounds \
                -init cgtools::espresso::replica_init -swap cgtools::espresso::replica_swap \
                -perform cgtools::espresso::replica_perform -info comm
            #		-perform cgtools::espresso::replica_perform -info all 
        } else {
            parallel_tempering::main -connect $replica_connect -init cgtools::espresso::replica_init \
                -swap cgtools::espresso::replica_swap \
                -perform cgtools::espresso::replica_perform -info comm
            #		-perform cgtools::espresso::replica_perform -info all 
        }
    } else {
        # by default, Replica exchange is off
        
        # Erase all the directory we're not resuming
        if { [info exists newcomp] } {
            foreach file [glob -nocomplain -directory $ident $ident.vmd*] {
                catch {exec rm -rf $ident}	
            }
        }

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

