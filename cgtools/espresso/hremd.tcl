# ::cgtools::espresso
#
#------------------------------------------------------------#
# Routines for the Hamiltonian Replica Exchange MD algorithm
# Author: Tristan Bereau
# Apr 14 2013
#------------------------------------------------------------#
#
namespace eval cgtools {
    namespace eval espresso {

        #################################################################################
        # proc hremd_init
        # Initializes the system before HREMD.
        # Arguments: 
        # - id:     identity of the system, only used when multiple replicas are dilivered
        #           to one processer
        # - lambda: Hamiltonian coupling

        proc hremd_init {id lambda} {

            #############################################################
            # local variable:   
            #       f_topo
            #                   starttime
            #                   startj
            #       startk
            #                   cutoff 
            #                   i
            #                   timingstart 
            #                   timingcurr
            #                   elapsedtime
            #############################################################
            global jjjjjj_[set lambda]
            global topology_[set lambda]

            variable topology

            variable timingstart
            variable this
            variable jjjjjj
            variable kkkkkk
            variable initial_lambda $lambda
            variable checkpointexists
            variable folder
            variable ::cgtools::forcefield::peptideb::softcore_flag
            variable ::cgtools::hremd
            # No softcore. We're scaling the Hbond strength only.  
            set softcore_flag 0
            if { $cgtools::hremd == 2 } {
	      # Softcore only for hremd==2
	      set softcore_flag 1
            }
            
            set this [namespace current]
            mmsg::send $this "Starting HREMD instance at lambda $lambda"

            # Set the output folders 
            set folder "$cgtools::outputdir/lambda$lambda"
            catch {exec mkdir $folder}
            set file_f [join [list $folder/observables.dat ] ""]
            set file_h [join [list $folder /histogram.dat ] ""]

            # Attempt to read a checkpoint file
            set checkpointexists [::cgtools::utils::readcheckpoint $folder ]

            ########## Start a new computation, if no checkpoint ##########
            if { !$checkpointexists } {
                # Creating observable data file
                set f [open $file_f w]
                puts $f "\# HREMD Observables at Hamiltonian coupling $lambda"
                puts $f "\# MD-Time \t Potential Energy \t replica id"
                close $f

                # Creating histogram data file
                set f [open $file_h w]
                puts $f "\# Histogram of Energies at Hamiltonian coupling $lambda"
                puts $f "\# Potential Energy \t Number of hits"
                close $f

                # Gnuplot script to look at all histograms at once
                set f [open "$cgtools::outputdir/histogram.gnu" w]
                puts $f "\# Gnuplot script to read histograms at all Hamiltonian couplings"
                puts $f ""
                puts $f "p " nonewline
                foreach lambda_i $cgtools::lambda_values {
                    puts $f "\'lambda$lambda_i/histogram.dat\' w histeps" nonewline
                    if {$lambda_i != [lindex $cgtools::lambda_values end]} {
                        puts $f ", " nonewline
                    }
                }
                close $f

                # Gnuplot script to look at all observables at once
                set f [open "$cgtools::outputdir/observables.energy.gnu" w]
                puts $f "\# Gnuplot script to read observables at all Hamiltonian couplings"
                puts $f ""
                puts $f "p " nonewline
                foreach lambda_i $cgtools::lambda_values {
                    puts $f "\'lambda$lambda_i/observables.dat\' w line" nonewline
                    if {$lambda_i != [lindex $cgtools::lambda_values end]} {
                        puts $f ", " nonewline
                    }
                }
                close $f

                set f [open "$cgtools::outputdir/observables.hremd.gnu" w]
                puts $f "\# Gnuplot script to read observables at all Hamiltonian couplings"
                puts $f ""
                puts $f "p " nonewline
                foreach lambda_i $cgtools::lambda_values {
                    puts $f "\'lambda$lambda_i/observables.dat\' u 1:3 w step" nonewline
                    if {$lambda_i != [lindex $cgtools::lambda_values end]} {
                        puts $f ", " nonewline
                    }
                }
                close $f

                # No checkpoint exists so we need to setup everything from scratch
                set startj 0
                set startk 0

                # Setup system
                ::cgtools::utils::setup_outputdir  $cgtools::outputdir -paramsfile $cgtools::paramsfile \
                    -tabdir $cgtools::tabledir -tabnames $cgtools::tablenames -coffdir $cgtools::overlapdir \
                    -coffnames $cgtools::overlapnames -readpdbdir $cgtools::readpdbdir \
                    -readpdbname $cgtools::readpdbname
                # Construct a directory for checkpoint backups inside $cgtools::outputdir/lambda$lambda 
                catch { exec rmdir $cgtools::outputdir/checkpoint_bak }    
                catch { exec mkdir $folder/checkpoint_bak }    

                # Set the box dimensions
                setmd box_l [lindex $cgtools::setbox_l 0] [lindex $cgtools::setbox_l 1] [lindex $cgtools::setbox_l 2]
                
                # Specify the bonded interactions
                #puts "bonded_parms $cgtools::bonded_parms"
                ::cgtools::utils::set_bonded_interactions $cgtools::bonded_parms

                # Specify any other non-bonded interactions
                if { [ catch { ::cgtools::utils::set_nb_interactions $cgtools::nb_interactions } ] } {
                    mmsg::warn $this "No non-bonded interactions used."
                }

                set cutoff [setmd max_cut] 
                ::mmsg::send $this "Cutoff used: $cutoff"

                # Initialize variable moltypelists in the namespace ::cgtools::utils
                ::cgtools::utils::initmoltypeskey $cgtools::moltypelists 

                # Initialize topology
                set topology [::cgtools::generation::generate_system $cgtools::system_specs $cgtools::setbox_l]
                # Set the generated topology into the internals of espresso.

                ::cgtools::utils::set_topology $topology
                
                # See if there is any fixed molecules 
                set cgtools::trappedmols [::cgtools::generation::get_trappedmols]
                # Fix molecules if necessary
                if { $cgtools::trappedmols != -1 } {
                    ::cgtools::utils::trap_mols $cgtools::trappedmols
                }

                # set exclusions for the bonded particles 1: 2-bodyinteraction  2: 3-body interaction
                part auto_exclusions 1
                
                #Initialise Random Number Generator
                ::cgtools::utils::init_random $cgtools::nprocessors

                # ----------- Integration Parameters before warmup -----------#
                setmd periodic 1 1 1
                setmd time_step $cgtools::warm_time_step
                setmd skin      $cgtools::verlet_skin
                thermostat langevin $cgtools::warmup_temp $cgtools::langevin_gamma
                
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
                        part [expr $i] fix 0 0 0
                    }
                }

                # ----------- Integration Parameters before warmup without any fixed particles -----------#
                setmd time_step $cgtools::main_time_step
                thermostat langevin  $cgtools::warmup_temp $cgtools::langevin_gamma

                # Warm up without any fixed particle 
                ::mmsg::send $this "warming up again at  [setmd temp]"
                ::cgtools::utils::warmup $cgtools::free_warmsteps $cgtools::free_warmtimes $topology \
                    -startcap 1000 -outputdir $folder
                
                # Setup analysis
                ::cgtools::analysis::setup_analysis $cgtools::analysis_flags -outputdir $folder \
                    -g $cgtools::mgrid -str $cgtools::stray_cut_off
                
                thermostat langevin $cgtools::systemtemp $cgtools::langevin_gamma

                # Reset the time to a starttime (usually zero) after warmup
                setmd time $cgtools::startmdtime 
            }
            
            # Resume a computation, if exists checkpoint
            if { $checkpointexists } {
                # A checkpoint exists so all we need to do is reset the moltypelists, topology and setup analysis again
                
                set topology [set topology_[set lambda]]
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

                
                set jjjjjj [set jjjjjj_[set lambda]]

                set kkkkkk [expr $cgtools::analysis_write_frequency - 1]

                # Make sure that we start exactly from where the checkpoint was written
                set startj [set jjjjjj]
                set startk [expr [set kkkkkk] + 1]
            }

            set timingstart [clock clicks -milliseconds]
            set jjjjjj $startj
            set kkkkkk $startk
            ### ONly necessary if each instance handles more than one configuration, 
            ### e.g. 300 replicas in 10 parallel processers
            #global config
            #set config($id) "{[part]} [setmd time]"
        }

        #################################################################################
        # proc hremd_perform
        # performs MD integration
        # Arguments : - id:   identity of system
        #         - lambda: Hamiltonian coupling
        proc hremd_perform {id lambda} {

            global jjjjjj_[set lambda]
            global topology_[set lambda]

            variable topology

            variable timingstart
            variable this
            variable jjjjjj
            variable kkkkkk

            variable initial_lambda

            ### ONly necessary if each instance handles more than one configuration, 
            ### e.g. 300 replicas in 10 parallel processers
            #global config
            #foreach p [lindex $config($id) 0] { eval part $p }
            #setmd time [lindex $config($id) 1]


            set folder "$cgtools::outputdir/lambda$lambda"
            set file_f [join [list $folder /observables.dat ] ""]
            set file_h [join [list $folder /histogram.dat   ] ""]
            
            #Main Integration                                          #
            #----------------------------------------------------------#

            # ----------- Integration Parameters after warmup -----------#
            setmd time_step $cgtools::main_time_step
            thermostat langevin $cgtools::systemtemp $cgtools::langevin_gamma

            if { $cgtools::thermo == "DPD" } {
                thermostat off
                set dpd_r_cut [setmd max_cut]
                thermostat set dpd $cgtools::systemtemp $cgtools::dpd_gamma $dpd_r_cut
                mmsg::send $this "DPD thermostat has been set"
                mmsg::send $this "Thermostat is: [thermostat]"
            }
            if { $cgtools::npt == "on" } {
                integrate set npt_isotropic $cgtools::p_ext $cgtools::piston_mass 1 1 0
                mmsg::send $this "npt integrator has been set"
                flush stdout
                #-cubic_box
                thermostat set npt_isotropic $cgtools::systemtemp  $cgtools::gamma_0  $cgtools::gamma_v
            }

            # Set force field with Hamiltonian coupling $lambda
            ::cgtools::forcefield::update_peptide_ff $lambda

            mmsg::send $this "run [set kkkkkk] at time=[setmd time]"

            # Call all of the analyze routines that we specified when setting up our analysis
            ::cgtools::analysis::do_analysis

            # If kkkkkk is a multiple of analysis_write_frequency then write the analysis results to file
            if { [expr [set kkkkkk] + 1] % $cgtools::analysis_write_frequency ==0 } {
                ::cgtools::analysis::print_averages
                #::cgtools::utils::update_force $rdfcglist $rdfaalist $tabledir $tablenames
            }

            # If kkkkkk is a multiple of write_frequency then write out a full particle configuration
            if { [expr [set kkkkkk] + 1] % $cgtools::write_frequency ==0 } {
                # polyBlockWrite "$folder/$cgtools::ident.[format %04d [set jjjjjj]].out" \
                #     {time box_l npt_p_diff } \
                #     {id pos type mass v f molecule} 
                # mmsg::send $this "wrote file $folder/$cgtools::ident.[format %04d [set jjjjjj]].out " 
                # flush stdout

                if { $cgtools::use_vmd == "offline" } {
                    ::cgtools::utils::writepdb_charmm \
                        "$folder/$cgtools::ident.vmd[format %04d [set jjjjjj]].pdb" \
			$topology -periodbox 1
                }

                incr jjjjjj

                # # Write a checkpoint to allow restarting.  Overwrites previous checkpoint
                # set jjjjjj_$lambda [set jjjjjj]
                # set topology_$lambda [set topology]
                # mmsg::send $this "setting checkpoint_$lambda [set kkkkkk] [setmd time] [set jjjjjj]"   
                # catch { exec rm -f $folder/checkpoint.latest.chk} 
                # checkpoint_set "$folder/checkpoint.latest.out"
                # # Try to copy a checkpoint to the backup checkpoint folder
                # # Usefull if the program crashes while writing a checkpoint
                # if { [ catch { exec cp -f $folder/checkpoint.latest.out \
                #                    $folder/checkpoint_bak/checkpoint.latest.out } ] } {
                #     mmsg::warn $this "warning: couldn't copy backup checkpoint"
                # }
            }
            #end of if { [expr [set kkkkkk] + 1] % $cgtools::write_frequency ==0 }

            # Do the real work of integrating equations of motion
            mmsg::send $this "starting integration: run $cgtools::hremd_timestep steps"
            integrate $cgtools::hremd_timestep

            ## get the replica ID in parallel computation
            set label 0
            foreach lnow $cgtools::lambda_values {
                if {$lnow == $initial_lambda} { break }
                incr label
            }

            set current_energy [expr [analyze energy total]-[analyze energy kinetic]]

            ::cgtools::utils::append_obs $file_f $current_energy $label
            ::cgtools::utils::write_histogram $file_h $current_energy        
            
            incr kkkkkk

            # Set the elapsed CPU time in computation, do not count that used for warm up
            set timingcurr [clock clicks -milliseconds]
            set elapsedtime [expr  $timingcurr - $timingstart]
            ::mmsg::send $this "elapsed time: $elapsedtime"

            ### ONly necessary if each instance handles more than one configuration, 
            ### e.g. 300 temperatures in 10 parallel processers
            #set config($id) "{[part]} [setmd time]"

        }

        #################################################################################
        # proc hremd_swap
        # Calculates the swap probabilities between two replicas
        # Arguments : - id: identity of the system to evaluate
        #             - Hamiltonian coupling of system 1 : lambda1
        #             - Hamiltonian coupling of system 2 : lambda2
        # Returns two probabilities.
        proc hremd_swap {id lambda1 lambda2} {
            ### ONly necessary if each instance handles more than one configuration, 
            ### e.g. 300 replicas in 10 parallel processers
            #global config
            #foreach p [lindex $config($id) 0] { eval part $p }

            # calculate potential energies for both force fields
            ::cgtools::forcefield::update_peptide_ff $lambda1
            set epot1 [expr [analyze energy total] - [analyze energy kinetic]]
            ::cgtools::forcefield::update_peptide_ff $lambda2
            set epot2 [expr [analyze energy total] - [analyze energy kinetic]]
            
            return "[expr $epot1/$cgtools::systemtemp] [expr $epot2/$cgtools::systemtemp]"
            
        }
    }
}
