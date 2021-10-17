# Peptide Glycosylation Starting Structure Generation
This repository contains the PyRosetta script for building the starting structure of an 
elongation reaction of mucin-type (O-)glycosylation reaction. 
Currently it has been tested for core-1 to core-2 chain formation
reactions. 

The PyRosetta script can be used directly from command line using ```python3 PeptideGlycosylationStartingStructureGenerator.py <necessary flags>``` or by importing the ```.py``` file in another PyRosetta instance and then modifying the parameters from the class itself.

# Code: 
```PeptideGlycosylationStartingStructureGenerator.py```

# Flags
```-pepseq```:      Input the one letter amino acid code of the acceptor peptide sequence.

```-sugar```:       Input the name of the sugar you want to use, e.g. core1.

```-gly_loc```:     Input the index (starting from 1) of the amino acid in the acceptor peptide sequence that will be glycosylated by the input sugar.

```-ref```:         If you want the Benchmarking to be done, input the location of the reference pdb file.

```-rI```:           Input if you want to repack the interface region of the acceptor peptide and the enzyme.

```-enzyme```:      Input the location of the enzyme + donor pdb file.

```-cnstr```:       Input the location of the constraints file.

```-decoy```:       Input the number of independent decoys the user wants.

```-out_pdb```:     Enter True, if you want the pdb of the output file.

```-out_dir```:     Input the direcotory where you want your output_pdb file to be saved.

```-out_name```:    Input the name of the output_pdb file.

# Execution Example
## From Command Line 
```cd EXAMPLE```
```python3 ../PeptideGlycosylationStartingStructureGenerator.py -pepseq ETTSHST -sugar core1 -gly_loc 3 -ref 3OTK-closed-monomer-alpha-GlcNAc_2GAM-GalBGalNAc.pdb -rI True -enzyme 3OTK-closed-monomer-alpha-GlcNAc-S217C_0005_598_manually_removed.pdb -cnstr 3OTK_constraints_file_ETTSHST.cst -decoy 2 -out_pdb True```

This command builds a ETTSHST peptide, glycosylates the third residue of that peptide using core-2 sugar and makes the starting structure of a closed monomer of core 2 beta1,6-n-acetylglucosaminyltransferase (RCSB id: 3otk).

Output: ```PeptideGlycosylationStartingStructureGenerator_ETTSHST.txt```

## Using PyRosetta Script
```cd EXAMPLE```
```python3 run_example.py```

Does the same job as above, just with STP peptide.

Output: ```PeptideGlycosylationStartingStructureGenerator_STP.txt```
