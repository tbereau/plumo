# Hello Peptide-Membrane Simulation #

Here's a quick guide to run your first peptide-membrane simulation.

Head over to the the example scripts
```
cd /PATH/TO/plumo/example.scripts
```

Copy the `walp.23` directory outside the PLUM directory, to avoid creating new files in the `plumo` subdirectory. Go to the new directory.

Inside `walp.23`, the `walp23-folded.tcl` describes the configuration file to run the folded, transmembrane WALP23 peptide in a small patch of POPC lipids at physiological conditions.

Run the simulation by calling
```
Espresso /PATH/TO/plumo/cgtools/cgtoolsmain.tcl walp23-folded.tcl
```

If you're architecture supports message passing interface (MPI) for parallel calculations, use the following to speed up the calculation on, e.g., 4 cores
```
mpirun -np 4 Espresso /PATH/TO/plumo/cgtools/cgtoolsmain.tcl -n 4 walp23-folded.tcl
```
Note that you'll need a much larger membrane to justify more than 4 cores per simulation.

The simulation should generate a subdirectory `walp23-folded` with a number of files:
  * `time_vs_*` are time evolutions of the observables you've asked PLUM to compute (specified in the config file)
  * `vmd_animation.script` is a `VMD` script to automatically load the PSF and PDB files of the run (for more information on this script, see https://github.com/tbereau/peptideB/wiki/Installation of the [peptideB](https://github.com/tbereau/peptideB) package)
  * `warm.*` correspond to the initial warmup phase