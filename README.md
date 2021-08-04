# O-glycosylation Starting Structure Formation
This repository contains the PyRosetta script for building the starting structure of an 
elongation reaction of mucin-type (O-)glycosylation reaction. 
Currently it has been tested for core-1 to core-2 chain formation
reaction. 

# Use:
Code: ```O_glycosylation_starting_structure_formation.py```
Flags: 
```-pepseq```:      Input the one letter amino acid code of the acceptor peptide sequence.

```-sugar```:       Input the name of the sugar you want to use, e.g. core1.

```-gly_loc```:     Input the index (starting from 1) of the amino acid in the acceptor peptide sequence that will be glycosylated by the input sugar.

```-ref```:         If you want the Benchmarking to be done, input the location of the reference pdb file.

```-rI``:           Input if you want to repack the interface region of the acceptor peptide and the enzyme.

```-enzyme```:      Input the location of the enzyme + donor pdb file.

```-cnstr```:       Input the location of the constraints file.

```-decoy```:       Input the number of independent decoys the user wants.

```-out_pdb```:     Enter True, if you want the pdb of the output file.

```-out_dir```:     Input the direcotory where you want your output_pdb file to be saved.

```-out_name```:    Input the name of the output_pdb file.