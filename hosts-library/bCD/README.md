## Adding beta-cyclodextrin (bCD) systems to the GHOAT workflow

Here we show how to configure GHOAT.py for guest-host calculations on the bCD host, using the coordinate and restart files from the proposed benchmark sets at https://github.com/MobleyLab/benchmarksets [1].

### Partial charges and atom types

The .mol2 files for the bCD host and its proposed guests will provide the partial charges and atom types, and can be found inside the following folder from the benchmark distribution:

*benchmarksets-master/input_files/cd-set2/mol2*

Here we use modified versions of these host and guest mol2 files, so that they are in the correct format for use with GHOAT. We describe the modifications below.

#### Host mol2 file

The modified host file (host-bcd.mol2) has only one host residue (MGO), and is already included in the  *hosts-library/bCD/parameters* folder from the GHOAT distribution. The parameters are identical to the ones provided by the benchmarks, both obtained from Ref. [2], but they are adjusted so GHOAT can read them and correctly build the systems.


#### Guest mol2 files

The guest mol2 files also need to be modified in order to match the system parameter (.prmtop) files from [1]. Here the atom types need to be converted to the General Amber Force Field (GAFF) format. This can be done for guest one using the *antechamber* program from Ambertools: 

antechamber -i guest-1.mol2 -fi mol2 -o guest-1-gaff.mol2 -fo mol2 -s 2 -at gaff -dr no

, and similarly for all other guests. Once in the correct format, the guests and host parameter files are then copied to the *GHOAT/parameters* folder, so they can be used in the binding free energy calculations.

### Additional files

A few additional files are needed to build the bCD systems, since this host has multiple identical residues and uses parameters from Ref. [2]. They are already included in the *hosts-library/bCD/parameters* and *hosts-library/bCD/amber_files* folders, and should be copied to the *GHOAT/parameters* and *GHOAT/amber_files* folders, respectively. More details on these files can be found in the user guide.  

### Initial structures

The needed structure (or restart) files can be found inside the folders:

*benchmarksets-master/input_files/cd-set2/prmtop-rst7*

To generate the complex structure file for each guest bound to the bCD host, the *cpptraj* program from Ambertools can be used. Type, inside the prmtop-rst7 folders above:

cpptraj -p bcd-1-s.prmtop -y bcd-1-s.rst7 -x host-bcd-guest-1s.pdb > cpptraj.log

, for calculations on the bCD host in complex with guest number one in the 's' orientation. Transfer the generated structure (host-bcd-guest-1s.pdb) to the *GHOAT/structures* folder inside the GHOAT package.

### Running the calculations

With the files above and the aCD-specific input parameters included in the host-aCD.in file, the calculations can now be performed as shown in the example from the GHOAT.py tutorial. Additional information can be found in the user guide.

### References

[1] D.L. Mobley, G. Heinzelmann, N.M. Henriksen, and M.K. Gilson. "Predicting binding free energies: Frontiers and benchmarks (a perpetual review)" - https://escholarship.org/uc/item/9p37m6bq.

[2] C. Cézard, X. Trivelli, F. Aubry, F. Djedaïni-Pilard and F.-Y. Dupradeau (2011).  “Molecular dynamics studies of native and substituted cyclodextrins in different media: 1. Charge derivation and force field performances”. Physical Chemistry Chemical Physics, **13**, 15103-15121.