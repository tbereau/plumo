# Generate initial structures #

## Membrane ##

The directory
```
/PATH/TO/plumo/example.scripts/membrane/readfiles
```
contains a number of pre-equilibrated membrane patches.  Information about their (approximate) size should be included in the header. We recommend simulations in the NPT (constant temperature and pressure) ensemble with 0 lateral tension.


## Peptides ##

We describe two methods to import peptide structures:
  * import an existing structure from an existing PDB file (e.g., atomistic simulation or experimental structure)
  * generate a conformation from the peptide coarse-grained model itself

### Import existing atomistic/experimental structure ###

Any PDB file can a priori be loaded in PLUMO. It will, however, require some adpatation. It needs to be stripped down such that only the following atoms are present: `N, CA, CB, C, O`, where all atoms but ` CB` belong to the backbone, and `CB` is the (first or "beta") carbon of the side chain. They need to be listed in that order of appearance. PLUMO does not read the information about the atom IDs, and thus the PDB lines can be swapped/deleted to achieve the appropriate format.  Here's a snippet of a PDB file that matches the desired format
```
ATOM   1942 N    GLY     2      50.382   2.973 105.699  1.00  0.00       P1
ATOM   1943 CA   GLY     2      50.906   2.805 104.263  1.00  0.00       P1
ATOM   1944 CB   GLY     2      51.714   3.204 104.012  1.00  0.00       P1
ATOM   1945 C    GLY     2      50.052   1.866 103.299  1.00  0.00       P1
ATOM   1946 O    GLY     2      48.846   1.608 103.231  1.00  0.00       P1
ATOM   1947 N    TRP     3      50.899   1.289 102.430  1.00  0.00       P1
ATOM   1948 CA   TRP     3      50.570   0.213 101.504  1.00  0.00       P1
ATOM   1949 CB   TRP     3      48.885   0.005 101.660  1.00  0.00       P1
ATOM   1950 C    TRP     3      50.888   0.316 100.051  1.00  0.00       P1
ATOM   1951 O    TRP     3      51.078   0.001  99.265  1.00  0.00       P1
```
where GLY and TRP correspond to residue numbers 2 and 3 of the peptide.

### Generate conformation ###

One can alternatively efficiently generate conformations of peptide(s) in **water** by using the [`peptideB`](https://github.com/tbereau/peptideB) package, which implements the peptide-only force field `[`[Bereau and Deserno (2009)](http://scitation.aip.org/content/aip/journal/jcp/130/23/10.1063/1.3152842)`]`. Output PDB files from this package can be imported into PLUMO once amide hydrogens have been deleted.  The following script will automatically delete all hydrogen entries
```
awk -f /PATH/TO/plumo/misc.scripts/clean_h.awk file.pdb
```

### Translate peptide in the _z_ direction ###

Using one of the two abovementioned methods, one may need to finetune height of the peptide(s) with respect to the membrane (e.g., outside, almost adsorbing, or inserted).  This can be readily achieved using the awk script
```
awk -f /PATH/TO/plumo/misc.scripts/translate_z.awk -v z=10 file.pdb
```
where we've translated the system by _delta z = 10 Angstroem_.