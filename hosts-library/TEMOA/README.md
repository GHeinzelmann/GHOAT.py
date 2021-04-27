## Adding tetra-endo-methyl octa-acid (TEMOA) systems to the GHOAT workflow

Here we show how to configure GHOAT.py for guest-host calculations on the TEMOA host, using the initial coordinates and parameters from the proposed benchmark sets at https://github.com/MobleyLab/benchmarksets [1].

### Partial charges and atom types

The .mol2 files for the TEMOA host and its proposed guests will provide the partial charges and atom types, and can be found inside the following folder from the benchmarks distribution above:

*benchmarksets-master/input_files/gdcc-set1/mol2*

The host and guests .mol2 files needed for the binding free energy calculations have to be copied to the *GHOAT/parameters* folder inside the GHOAT package.  Their names should match the host and guest names in the GHOAT input file (see the tutorial).

### Initial structures

The initial coordinates files (.rst7) with full system parameter files (.prmtop) can be found inside the folders:

*benchmarksets-master/input_files/gdcc-set1/prmtop-rst7*

To generate the complex structure file for each guest bound to the TEMOA host in .pdb format, the *cpptraj* program from Ambertools can be used. Type, inside the prmtop-rst7 folders above:

cpptraj -p temoa-3.prmtop -y temoa-3.rst7 -x host-temoa-guest-3.pdb > cpptraj.log

, for calculations on the TEMOA host in complex with guest number 3. Copy the generated structure (host-temoa-guest-3.pdb) to the *GHOAT/structures* folder inside the GHOAT package. 

### Running the calculations

With the files above and the TEMOA-specific GHOAT input parameters included in the host-TEMOA.in file, the calculations can now be performed as shown in the example from the GHOAT.py tutorial. Additional information can be found in the GHOAT tutorial and user guide.


### References

[1] D.L. Mobley, G. Heinzelmann, N.M. Henriksen, and M.K. Gilson. "Predicting binding free energies: Frontiers and benchmarks (a perpetual review)" - https://escholarship.org/uc/item/9p37m6bq.
