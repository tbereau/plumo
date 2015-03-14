# Configuration file #

The configuration file describes the simulation protocol to PLUMO. Have a look at examples in the directory
```
/PATH/TO/plumo/example.scripts
```
As an example, we will follow `/PATH/TO/plumo/example.scripts/walp.23/walp23-folding.tcl`.

The following describes a number of important options for a simulation:
### General information ###
```
set outputdir "walp23-folded"
```
will assign the name of the output directory.
```
set readpdbdir "./readfiles"
```
is the location of the input PDBs.

### Input files ###
```
set numlipids 72 
set lipidfile "bilayer_"
append lipidfile $numlipids
append lipidfile ".init.crd"

set numpeptides 1
set resilist1 {NAP GLY TRP TRP LEU ALA LEU ALA \
    LEU ALA LEU ALA LEU ALA LEU ALA LEU ALA LEU ALA LEU TRP TRP ALA CAP}
set peptide1file "walp23.inserted.folded.crd"

set readpdbname [list "$lipidfile" "$peptide1file"]
```
will load a membrane structure from `CRD` file `bilayer_72.init.crd` and read 1 peptide in the `CRD` file `walp23.inserted.folded.crd`. Note that `PDB` files can also be read.  The sequence of the peptide in 3-letter code needs to be provided. Note that we've added N- and C-termini caps, `NAP` and `CAP`, respectively.

### Geometry and composition ###
```
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
```
generates a lipid membrane with $numlipids (i.e., 72 here) lipids from file `$ident/$lipidfile`. Lipids have component 0 in PLUM, while peptides have component 1.

### Box size ###
```
# Set the box size 
# Notice, should be the same as the values in the bilayer readfile
set lengthx 51.0
lappend setbox_l $lengthx
lappend setbox_l $lengthx
lappend setbox_l 170.0 
unset lengthx
```
Sets the box size to 51\*51\*170.

## Thermostat/barostat ##
```
# Thermostat
set langevin_gamma 0.2; # Langevin friction
set systemtemp 1.; # temperature in units of kT_body
# Barostat
set npt "on"; # "on/off" flag for constant-pressure
set p_ext 0.000; # reference pressure
set piston_mass 0.0005; # piston mass
set gamma_0 0.2; # Langevin friction
set gamma_v 0.00004; # box friction
```

## Warmup parameters ##
```
set warm_time_step 0.1; # warmup time step
set warmup_temp 1.; # warmup temperature
set warmup_freq 5; # warmup frequency of file storage
set warmsteps 5; # Number of steps to integrate at each warmup phase
set warmtimes 20; # Number of warmup phases
set free_warmsteps 1; # Second stage of warmup.
set free_warmtimes 20
```
time step for warmup can be set to `0.1 tau`, since the peptide will be frozen anyway.

## Other parameters ##
```
set main_time_step 0.01; # main time step. 
set verlet_skin 2.0 
# Simulation time (in units of tau) while peptide is frozen
set fix_time 10.0
# The number of steps to integrate with each call to integrate
set int_steps   1000
# The number of times to call integrate
set int_n_times 1000
# backup frequency (store output PDB)
set write_frequency 10
# analysis frequency (compute observables)
set analysis_write_frequency 1
```
The main time step can't be larger than 0.01 for a stable integration of the peptide degrees of freedom.