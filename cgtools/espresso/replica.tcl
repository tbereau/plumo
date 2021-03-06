# ::cgtools::espresso
#
#------------------------------------------------------------#
# Routines for the replica exchange algorithm (selfassembly)
# Author: Zun-Jing Wang 
# Oct 13 2009
#------------------------------------------------------------#
#
namespace eval cgtools {
    namespace eval espresso {

        #################################################################################
        # proc replica_init
        # Initializes the system before replica exchange.
        # Arguments: 
        # - id:   identity of the system, only used when multiple temperatures are dilivered to one processer
        # - temp: temperature

        proc replica_init {id temp} {

            #############################################################
            # local variable:   
            #           f_topo
            #                   starttime
            #                   startj
            #           startk
            #                   cutoff 
            #                   i
            #                   timingstart 
            #                   timingcurr
            #                   elapsedtime
            #############################################################
            global jjjjjj_[set temp]
            global topology_[set temp]

            variable topology

            variable timingstart
            variable this
            variable jjjjjj
            variable kkkkkk
            variable initial_temp $temp
            variable prev_temp $temp
            variable checkpointexists
            variable folder

            variable ::cgtools::temp_peptide
            variable ::cgtools::generation::peptide_parts            

            set this [namespace current]
            mmsg::send $this "Starting replica exchange instance at temperature $temp"

            # Resume old simulation
            set mdinit $cgtools::startmdtime
            set resuming 0
            # Set the output folders 
            set folder "$cgtools::outputdir/temp$temp"
            if { [catch {exec mkdir $folder}] } {
                # Directory exists. Try to resume the simulation (see below).
                set resuming 1
            } else {
                puts "No existing directory. Fresh start."
            }
            set file_f [join [list $folder/observables.dat ] ""]
            set file_h [join [list $folder /histogram.dat ] ""]

            # Attempt to read a checkpoint file
            set checkpointexists [::cgtools::utils::readcheckpoint $folder ]

            ########## Start a new computation, if no checkpoint ##########
            if { !$checkpointexists } {
                # No checkpoint exists so we need to setup everything from scratch
                set startj 0
                set startk 0

                if { $resuming == 1 } {
                    set resuming 0
                    variable ::cgtools::pdb_resume
                    # Resuming old simulation. Look for latest PDB file
                    set tclfiles [glob -nocomplain -directory $folder *.pdb]
                    set lastFile ""
                    if { [llength $tclfiles] > 0} {
                        set lastDate "000000000"
                        foreach f $tclfiles {
                            set fdate [file mtime $f]
                            if { [expr {$fdate > $lastDate}] && [string first "warm" $f] == "-1" } {
                                set lastFile $f
                                set lastDate $fdate
                            }
                        }
                    }
                    if { $lastFile != "" } {
                        # Read rough time step from file id
                        set vmdId [string first $cgtools::ident.vmd $lastFile]
                        if { $vmdId != -1 } {
                            set vmdLgth [string length "$cgtools::ident.vmd"]
                            set initTime [scan [string range $lastFile [expr $vmdId+$vmdLgth] end-4] %d]
                            set startk [expr $initTime + 1]
                            set startj [expr $initTime + 1]
                            set mdinit [expr $mdinit + \
                                ($initTime+1) * $cgtools::replica_timestep * $cgtools::main_time_step * \
                                    $cgtools::write_frequency]
                        }
                        set ::cgtools::pdb_resume $lastFile
                        set resuming 1
                        puts "Resuming simulation from PDB: $pdb_resume"
                        puts "init time step $mdinit" 
                    } else {
                        puts "No PDB to resume from"
                    }
                }

                if { $resuming == 0 } {
                    # Creating observable data file
                    set f [open $file_f w]
                    puts $f "\# Replica Observables at Temperature $temp"
                    puts $f "\# MD-Time \t Potential Energy \t replica id"
                    close $f

                    # Creating histogram data file
                    set f [open $file_h w]
                    puts $f "\# Histogram of Energies at Temperature $temp"
                    puts $f "\# Potential Energy \t Number of hits"
                    close $f

                    # Gnuplot script to look at all histograms at once
                    set f [open "$cgtools::outputdir/histogram.gnu" w]
                    puts $f "\# Gnuplot script to read histograms of all temperatures"
                    puts $f ""
                    puts $f "p " nonewline
                    foreach temperature $cgtools::replica_temps {
                        puts $f "\'temp$temperature/histogram.dat\' w histeps" nonewline
                        if {$temperature != [lindex $cgtools::replica_temps end]} {
                            puts $f ", " nonewline
                        }
                    }
                    close $f

                    # Gnuplot script to look at all observables at once
                    set f [open "$cgtools::outputdir/observables.energy.gnu" w]
                    puts $f "\# Gnuplot script to read observables of all temperatures"
                    puts $f ""
                    puts $f "p " nonewline
                    foreach temperature $cgtools::replica_temps {
                        puts $f "\'temp$temperature/observables.dat\' w line" nonewline
                        if {$temperature != [lindex $cgtools::replica_temps end]} {
                            puts $f ", " nonewline
                        }
                    }
                    close $f

                    set f [open "$cgtools::outputdir/observables.replica.gnu" w]
                    puts $f "\# Gnuplot script to read observables of all temperatures"
                    puts $f ""
                    puts $f "p " nonewline
                    foreach temperature $cgtools::replica_temps {
                        puts $f "\'temp$temperature/observables.dat\' u 1:3 w step" nonewline
                        if {$temperature != [lindex $cgtools::replica_temps end]} {
                            puts $f ", " nonewline
                        }
                    }
                    close $f

                }


                # Setup the output directory by creating it and copying forcetables and overlapped potcoffs to it
                #::cgtools::utils::setup_outputdir  $cgtools::outputdir -paramsfile $cgtools::paramsfile \
                    -tabdir $cgtools::tabledir -tabnames $cgtools::tablenames -coffdir $cgtools::overlapdir \
                    -coffnames $cgtools::overlapnames -readpdbdir $cgtools::readpdbdir \
                    -readpdbname $cgtools::readpdbname
                ::cgtools::utils::setup_outputdir  $cgtools::outputdir -paramsfile $cgtools::paramsfile \
                    -tabdir $cgtools::tabledir -tabnames $cgtools::tablenames -coffdir $cgtools::overlapdir \
                    -coffnames $cgtools::overlapnames -readpdbdir $cgtools::readpdbdir \
                    -readpdbname $cgtools::readpdbname

                # Set the box dimensions
                setmd box_l [lindex $cgtools::setbox_l 0] [lindex $cgtools::setbox_l 1] [lindex $cgtools::setbox_l 2]
                
                # Specify the bonded interactions
                #puts "bonded_parms $cgtools::bonded_parms"
                ::cgtools::utils::set_bonded_interactions $cgtools::bonded_parms

                # Specify any other non-bonded interactions
                if { [ catch { ::cgtools::utils::set_nb_interactions $cgtools::nb_interactions } ] } {
                    mmsg::send $this "no non-bonded interactions used"
                }

                set cutoff [setmd max_cut] 
                puts "max_cut is $cutoff"

                # Initialize variable moltypelists in the namespace ::cgtools::utils
                ::cgtools::utils::initmoltypeskey $cgtools::moltypelists 
                #puts "$cgtools::moltypelists"
                #puts "$cgtools::ident"
                #puts "$cgtools::system_specs"
                #puts "$cgtools::setbox_l"

                # Initialize topology
                set topology [::cgtools::generation::generate_system $cgtools::system_specs $cgtools::setbox_l]
                # Set the generated topology into the internals of espresso.
                
                #puts "HelloHelloHelloHelloHelloHelloHelloHelloHello"

                ::cgtools::utils::set_topology $topology
                #puts "$topology"
                
                # See if there is any fixed molecules 
                set cgtools::trappedmols [::cgtools::generation::get_trappedmols]
                if { $resuming == 0 } {
                    # Fix molecules if necessary
                    if { $cgtools::trappedmols != -1 } {
                        ::cgtools::utils::trap_mols $cgtools::trappedmols
                    }
                }

                # set exclusions for the bonded particles 1: 2-bodyinteraction  2: 3-body interaction
                part auto_exclusions 1
                
                #Initialise Random Number Generator
                ::cgtools::utils::init_random $cgtools::nprocessors



                # ----------- Integration Parameters before warmup -----------#
                setmd periodic 1 1 1
                setmd time_step $cgtools::warm_time_step
                setmd skin      $cgtools::verlet_skin
                if { $temp_peptide > 0. } {
                    # Set all peptide particles at temperature $temp
                    for {set i 0} {$i < [llength $peptide_parts]} {incr i} {
                        part [lindex $peptide_parts $i] temp $temp
                    }
                    # The rest of the system is at $systemtemp
                    thermostat langevin $cgtools::systemtemp $cgtools::langevin_gamma
                } else {
                    thermostat langevin $temp $cgtools::langevin_gamma    
                }
                
                
                # Set the topology and molecule information
                #----------------------------------------------------------#
                #write topology file
                set f_topo [open "$folder/$cgtools::ident.top" w]
                blockfile_write_topology $f_topo write topology   
                close $f_topo

            # Check if there are any extra vmdcommands and if not initialize a default
            ::cgtools::utils::initialize_vmd $cgtools::use_vmd $folder $cgtools::ident \
                $topology -extracommands $cgtools::vmdcommands

                #Perform the warm up integration
                #----------------------------------------------------------#
                # Warm up containing fixed particles 
                mmsg::send $this "warming up at [setmd temp]"
                ::cgtools::utils::warmup  $cgtools::warmsteps $cgtools::warmtimes $topology -cfgs 10 \
                    -outputdir $folder

                # If the membrane have any fixed particles, unfix them after warmup
                set cgtools::userfixedparts [::cgtools::generation::get_userfixedparts ]
                for {set i 0} { $i <  [setmd n_part] } {incr i} {
                    if { [lsearch $cgtools::userfixedparts $i ] == -1 } {
                        if { ([part $i print type] > 7 && $::cgtools::implicit_membrane == 1) || \
                                $::cgtools::implicit_membrane == 0 } {
                            part [expr $i] fix 0 0 0
                        }
                    }
                }

                # ----------- Integration Parameters before warmup without any fixed particles -----------#
                setmd time_step $cgtools::main_time_step
                thermostat langevin  $temp $cgtools::langevin_gamma

                # Warm up without any fixed particle 
                ::mmsg::send $this "warming up again at  [setmd temp]"
                ::cgtools::utils::warmup $cgtools::free_warmsteps $cgtools::free_warmtimes $topology \
                    -startcap 1000 -outputdir $folder
                
                # Setup analysis
                ::cgtools::analysis::setup_analysis $cgtools::analysis_flags -outputdir $folder \
                    -g $cgtools::mgrid -str $cgtools::stray_cut_off
                

                # Reset the time to a starttime (usually zero) after warmup
                setmd time $mdinit 

            }
            
            # Resume a computation, if exists checkpoint
            if { $checkpointexists } {
                # A checkpoint exists so all we need to do is reset the moltypelists, topology and setup analysis again
                
                set topology [set topology_[set temp]]
                #puts "$topology"
                ::cgtools::utils::initmoltypeskey $cgtools::moltypelists 
                ::cgtools::utils::read_topology "$folder/$cgtools::topofile" 

                # set exclustions for the bonded particles 1: 2-bodyinteraction  2: 3-body interaction
                part auto_exclusions 1
                
                # Setup analysis
                ::cgtools::analysis::setup_analysis $cgtools::analysis_flags -outputdir $folder \
                    -g $cgtools::mgrid -str $cgtools::stray_cut_off
                
                set topology [analyze set] 
                ::cgtools::utils::initialize_vmd $cgtools::use_vmd $folder $cgtools::ident $topology

                
                set jjjjjj [set jjjjjj_[set temp]]

                set kkkkkk [expr $cgtools::analysis_write_frequency - 1]

                # Make sure that we start exactly from where the checkpoint was written
                set startj [set jjjjjj]
                set startk [expr [set kkkkkk] + 1]
            }

            set timingstart [clock clicks -milliseconds]
            set jjjjjj $startj
            set kkkkkk $startk
            ### ONly necessary if each instance handles more than one configuration, 
            ### e.g. 300 temperatures in 10 parallel processers
            #global config
            #set config($id) "{[part]} [setmd time]"
        }

        #################################################################################
        # proc replica_perform
        # performs MD integration
        # Arguments : - id:   identity of system
        #         - temp: temperature
        proc replica_perform {id temp} {

            global jjjjjj_[set temp]
            global topology_[set temp]

            variable topology

            variable timingstart
            variable this
            variable jjjjjj
            variable kkkkkk

            variable initial_temp
            variable prev_temp

            variable ::cgtools::temp_peptide
            variable ::cgtools::generation::peptide_parts

            ### ONly necessary if each instance handles more than one configuration, 
            ### e.g. 300 temperatures in 10 parallel processers
            #global config
            #foreach p [lindex $config($id) 0] { eval part $p }
            #setmd time [lindex $config($id) 1]


            set folder "$cgtools::outputdir/temp$temp"
            set file_f [join [list $folder /observables.dat ] ""]
            set file_h [join [list $folder /histogram.dat   ] ""]
            
            #Main Integration                                          #
            #----------------------------------------------------------#

            # ----------- Integration Parameters after warmup -----------#
            setmd time_step $cgtools::main_time_step
            if { $temp_peptide > 0. } {
                # Set all peptide particles at temperature $temp
                for {set i 0} {$i < [llength $peptide_parts]} {incr i} {
                    part [lindex $peptide_parts $i] temp $temp
                }
                thermostat langevin $cgtools::systemtemp $cgtools::langevin_gamma
            } else {            
                thermostat langevin $temp $cgtools::langevin_gamma
            }

            # Rescale velocities

            global ::cgtools::analysis::n_particles
            set facvec [expr sqrt($temp/$prev_temp)]
            for { set part_no 0 } {$part_no < $n_particles } { incr part_no } {
                if { ($temp_peptide > 0. && [lsearch $peptide_parts $part_no] != -1) || \
                    $temp_peptide < 0.} {
                    set partvel [::cgtools::utils::scalevec [part $part_no print v] $facvec]
                    part $part_no v [lindex $partvel 0] [lindex $partvel 1] [lindex $partvel 2]
                }
            }
            set prev_temp $temp




            variable ::cgtools::implicit_membrane
            if { $::cgtools::implicit_membrane == 1 } {
                variable ::cgtools::forcefield::peptideb::softcore_flag
                set softcore_flag 1
                puts "Running implicit membrane"
                ::cgtools::forcefield::update_peptide_ff 0.0
                # Freezing all lipids -- loop over all particles and look for types {0-7}
                global ::cgtools::analysis::n_particles
                for { set part_no 0 } {$part_no < $n_particles} {incr part_no} {
                    set ptype [part $part_no print type]
                    if { $ptype < 8 } {
                        part $part_no fix 1 1 1
                    }
                }
            }


            if { $cgtools::thermo == "DPD" } {
                thermostat off
                set dpd_r_cut [setmd max_cut]
                if { $temp_peptide > 0. } {
                    thermostat set dpd $cgtools::systemtemp $cgtools::dpd_gamma $dpd_r_cut    
                } else {
                    thermostat set dpd $temp $cgtools::dpd_gamma $dpd_r_cut
                }
                mmsg::send $this "DPD thermostat has been set"
                mmsg::send $this "Thermostat is: [thermostat]"
            }

            variable ::cgtools::membrane_restraint
            if { $cgtools::npt == "on" && $::cgtools::implicit_membrane == 0 } {
                if { ($::cgtools::membrane_restraint == 1 && $temp == [lindex $cgtools::replica_temps 0]) \
                    || $::cgtools::membrane_restraint == 0} {
                    # membrane_restraint -> NPT only in lowest replica
                    integrate set npt_isotropic $cgtools::p_ext $cgtools::piston_mass 1 1 0
                    mmsg::send $this "npt integrator has been set"
                    flush stdout
                    #-cubic_box
                    if { $temp_peptide > 0. } {
                        thermostat set npt_isotropic $cgtools::systemtemp  $cgtools::gamma_0  $cgtools::gamma_v
                    } else {
                        thermostat set npt_isotropic $temp  $cgtools::gamma_0  $cgtools::gamma_v                    
                    }
                }
            }

            # mmsg::send $this "run [set kkkkkk] at time=[setmd time]"

            # Do the real work of integrating equations of motion
            # mmsg::send $this "starting integration: run $cgtools::replica_timestep steps"
            integrate $cgtools::replica_timestep

            ## get the id of temperatures in parallel computation
            set label 0
            foreach tnow $cgtools::replica_temps {
                if {$tnow == $initial_temp} { break }
                incr label
            }

            set current_energy [expr [analyze energy total]-[analyze energy kinetic]]

            ::cgtools::utils::append_obs $file_f $current_energy $label
            ::cgtools::utils::write_histogram $file_h $current_energy   

            # Setup analysis
            ::cgtools::analysis::setup_analysis $cgtools::analysis_flags -outputdir $folder \
                -g $cgtools::mgrid -str $cgtools::stray_cut_off


            # Call all of the analyze routines that we specified when setting up our analysis
            ::cgtools::analysis::do_analysis

            # If kkkkkk is a multiple of analysis_write_frequency then write the analysis results to file
            if { [expr [set kkkkkk] + 1] % $cgtools::analysis_write_frequency ==0 } {
                ::cgtools::analysis::print_averages
                #::cgtools::utils::update_force $rdfcglist $rdfaalist $tabledir $tablenames
            }

            # If kkkkkk is a multiple of write_frequency then write out a full particle configuration
            if { [expr [set kkkkkk] + 1] % $cgtools::write_frequency ==0 } {
                if { $cgtools::use_vmd == "offline" } {
                    ::cgtools::utils::writepdb_charmm \
                        "$folder/$cgtools::ident.vmd[format %04d [set jjjjjj]].pdb" \
                        $topology -periodbox 1
                }

                incr jjjjjj
            }
            #end of if { [expr [set kkkkkk] + 1] % $cgtools::write_frequency ==0 }
            
            ::cgtools::utils::update_midplane_pos $topology

            incr kkkkkk

            # Set the elapsed CPU time in computation, do not count that used for warm up
            set timingcurr [clock clicks -milliseconds]
            set elapsedtime [expr  ($timingcurr - $timingstart)/1000.]
            set widthrounds [string length $::cgtools::replica_rounds]
            ::mmsg::send $this [format "Frame: %[set widthrounds]d/%[set widthrounds]d; elapsed time: %9.0f seconds" \
                $kkkkkk $cgtools::replica_rounds $elapsedtime]

            ### ONly necessary if each instance handles more than one configuration, 
            ### e.g. 300 temperatures in 10 parallel processers
            #set config($id) "{[part]} [setmd time]"

        }

        #################################################################################
        # proc replica_swap
        # Calculates the swap probabilities between two replicas
        # Arguments : - id: identity of the system to evaluate
        #         - temperature of system 1 : temp1
        #             - temperature of system 2 : temp2
        # Returns two probabilities.
        proc replica_swap {id temp1 temp2} {
            ### ONly necessary if each instance handles more than one configuration, 
            ### e.g. 300 temperatures in 10 parallel processers
            #global config
            #foreach p [lindex $config($id) 0] { eval part $p }

            #set epot [expr [analyze energy total] - [analyze energy kinetic]]
            
            # Replica-Exchange Solute Tempering (Liu et al. PNAS _102_ (2005))
            # Sum up all energy contributions protein-protein and protein-lipid
            # lipid has IDs: 0-7, protein: 8-43
            set epot 0.0
            for {set beadprot 8} {$beadprot < 42} {incr beadprot} {
                for {set beadall 0} {$beadall < 8} {incr beadall} {
                    set epot [expr {$epot + 0.5*[analyze energy nonbonded $beadall $beadprot]}]
                }
                for {set beadall 8} {$beadall < 42} {incr beadall} {      
                    set epot [expr {$epot + [analyze energy nonbonded $beadprot $beadall]}]
                }
            }
            return "[expr $epot/$temp1] [expr $epot/$temp2]"
        }
    }
}
