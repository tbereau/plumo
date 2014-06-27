# 72 POPC bilayer with WALP23 peptide. HREMD simulation: decoupling of
# peptide interactions with membrane: 
#  lambda = 1: original Hamiltonian.
#  lambda = 0: peptide decoupled from membrane.
#
::mmsg::send [namespace current] "loading parameter file .. " nonewline
flush stdout

# Specify the name of the job <ident> the location of forcetables and
# the output and scripts directories, dir and filename of the initial positions
#set currentrundir "runcode"
set ident "walp23-folded"
set outputdir "./$ident"
set topofile "$ident.top"
set readpdbdir "./readfiles"

# lipid and peptide input files
set numlipids 72 
set lipidfile "bilayer_"
append lipidfile $numlipids
append lipidfile ".init.crd"
set numpeptides 1
set resilist1 {NAP GLY TRP TRP LEU ALA LEU ALA \
    LEU ALA LEU ALA LEU ALA LEU ALA LEU ALA LEU ALA LEU TRP TRP ALA CAP}
set peptide1file "walp23.inserted.folded.crd"

set readpdbname [list "$lipidfile" "$peptide1file"]
# Specify how we want to use vmd for visualization: allowable values
# are "interactive" "offline" "none".  interactive rarely works
set use_vmd "offline"

#########################################################
########## Specify the geometry and composition ##########
#########################################################
# --- Specify the lipid+protein geometry and composition ----#
### Set the geometry 
set geometry "geometry fromfile $ident/$lipidfile"
### Set the n_molslist 
#n_molslist structure:   { n_molslist {  { 0 72 } } }
set component 0
set n_molslist "n_molslist [list [list [list $component $numlipids]]]"
### Now bundle the above info into a list
lappend lipidspec $geometry
lappend lipidspec $n_molslist
unset geometry n_molslist

# Set the box size 
# Notice, should be the same as the values in the bilayer readfile
set lengthx 51;#54.53
lappend setbox_l $lengthx
lappend setbox_l $lengthx
lappend setbox_l 166.9462 
unset lengthx

# --- Specify geometry and compositions of the peptide #1 ----#
### Set the geometry 
set geometry "geometry fromfile $ident/$peptide1file"
### Set the n_molslist 
set component 1
set n_molslist "n_molslist [list [list [list $component $numpeptides]]]"
set sequence "sequence [list $resilist1]"
### Now bundle the above info into a list
lappend peptidespec $geometry
lappend peptidespec $n_molslist
lappend peptidespec $sequence
unset geometry n_molslist


# Now group the lipidspec with peptidespec into a list of such
# systems (we can have multiple systems if we like each with different
# composition of molecule types
lappend system_specs $lipidspec
lappend system_specs $peptidespec
unset lipidspec peptidespec


# Warmup parameters
#----------------------------------------------------------#
set warm_time_step 0.1
set warmup_temp 1.
set warmup_freq 5
set warmsteps 5;#50
set warmtimes 20
set free_warmsteps 1;#10
set free_warmtimes 20

# ------ Integration parameters -----------------#
set main_time_step 0.01
set verlet_skin 2.0 
# -------Constant Temperature-----------------#
set langevin_gamma 0.2
set systemtemp 1.
# -------Constant Pressure-----------------#
set npt "on"
set p_ext 0.000
set piston_mass 0.0005
set gamma_0 0.2
set gamma_v 0.00004


# -------DPD-----------------#
#set thermo "DPD"
#set dpd_gamma 1.0 

# Simulation time (in units of tau) while peptide is frozen
set fix_time 10.0

# -------- Set the espresso integration steps and times 
# The number of steps to integrate with each call to integrate
set int_steps   1000
# The number of times to call integrate
set int_n_times 1000

# -------- Frequency of backup and analysis
# backup frequency 
set write_frequency 10
# analysis frequency 
set analysis_write_frequency 1


##########################################################
# Analysis Parameters
#----------------------------------------------------------# 
set mgrid 8
set stray_cut_off 30.

# These are are parameters that will be passed to the setup_analysis
lappend analysis_flags boxl
lappend analysis_flags temperature

::mmsg::send [namespace current] "Parameter file parsing is done"
