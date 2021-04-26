## Adding cucurbit[7]uril (CB7) systems to the GHOAT workflow

Here we show how to configure GHOAT.py for guest-host calculations on the CB7 host, using the coordinate and restart files from the proposed benchmark sets at https://github.com/MobleyLab/benchmarksets.

### Partial charges and atom types

The .mol2 files for the CB7 host and its proposed guests will provide the partial charges and atom types, and can be found inside the following folders:

*benchmarksets-master/input_files/cb7-set1/mol2*

*benchmarksets-master/input_files/cb7-set2/mol2*

The host and guest mol2 files needed for the binding free energy calculations should be added to the *GHOAT/parameters* folder inside the GHOAT package.

### Initial structures

The needed structure (or restart) files can be found inside the folders:

*benchmarksets-master/input_files/cb7-set1/prmtop-rst7*

*benchmarksets-master/input_files/cb7-set2/prmtop-rst7*

To generate the complex structure file for each guest bound to the CB7 host, the *cpptraj* program from Ambertools can be used. After making sure *cpptraj* is in your path (running amber.sh), type inside the prmtop-rst7 folders above:

cpptraj -p cb7-1.prmtop -y cb7-1.rst7 -x host-cb7-guest-1.pdb > cpptraj.log

, for calculations on the CB7 host in complex with guest number one. Transfer the generated structure (host-cb7-guest-1.pdb) to the *GHOAT/structures* folder inside the GHOAT package. 

### Running the calculations

With the files above and the CB7-specific input parameters included in the host-CB7.in file, the calculations can now be performed as shown in the example from the GHOAT.py tutorial.