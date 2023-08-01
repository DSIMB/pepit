# PEPIT

\(C\) Frédéric Guyon, Gautier Moroy

**Non-sequential alignment of binding sites for fast peptide screening**


## Abstract

The PepIT program proposes peptides that may interact with a given protein. PepIT is based on a non-sequential alignment algorithm to identify peptide binding sites that share geometrical and physicochemical properties with regions on the target protein. PepIT compares the entire surface of the target protein with the peptide binding sites in the Propedia dataset, which contains more than 19,000 high-resolution protein-peptide structures. Once a peptide binding site similar to a portion of the protein surface is found, the peptide inserted into the binding site is repositioned on the corresponding portion of the protein surface. 

## Dependencies

The R program is required as Pepit is a R package.
A few R scripts are included that make use of the pepit R package.
No R programming is required.

## Installation

The pepit package is installed on github deposit https://github.com/DSIMB/pepit.

Invoke R and type:

```R
install.packages("devtools")
```

```R
library(devtools)
```

```R
install_github("DSIMB/pepit")
```

## Usage

Two scripts make use of the package pepit.

BuiltBSBank.R   builds a bank of protein-peptide complexes with the format expected by pepit

pepit.R performs a search of a binding site of the bank similar to the given protein target.


### Bank of complexes preparation

As input, the script takes a bank description file, that is a 3 columns file with a header, that can be read by
R as a data frame. 

The second input is a directory name where are written the ouput complexe pdb files.

Example:
```bash
mkdir Bank
```

For example, with the following file nanobank.dat
```bash
PDB Receptor.Chain Peptide.Chain
148l E S
1a07 A C
1a07 B D
1a08 A C
```

```bash
Rscript ./BuiltBSBank.R nanobank.dat Bank
```

the script produces 4 pdb files  `148lE:S.pdb`  `1a07A:C.pdb`  `1a07B:D.pdb`  `1a08A:C.pdb`
utiliser le script pepit.R

### Peptide search
use the script `pepit.R`

```bash
Rscript ./pepit.R target target_chain BSBank prefix
```

- target is a pdb file or a pdb id. In this case, `pepit.R` download the pdb file from the Protein Data Bank RCSB PDB.
- target_chain gives the receptor chain of the pdb target file. It can be a single character or a list of chains as H,L or A,B,C or "*" meaning all the protein chains (the double quote are mandatory)
- BSBank is binding site bank created by `BuiltBSBank.R`
- prefix is a string used as a prefix for the names of all the output files.

By default, `pepit.R` produces 1 output file : <prefix>.score

It is possible to change the default values of the program parameters at the beginning of the script file.

For example, if the parameter POSE is set to TRUE, an alignment file that gives the mapping between the protein target and the selected binding sites is output.
This file name is  <prefix>.al. PDB files with the selected peptides posed onto the target are also produced.

### Example

```bash
Rscript ./pepit.R 148l E Bank prefix test
```

output : test.score

```bash
index bs target precision bslen alen rmsd coverage meandist score
1 sampleBank/148lE:S.pdb 148l 1 12 12 0 1 0 144
```
Only one hit is found that is the binding site of 148l (148lE:S.pdb biding site of the peptide S onto the chain E of 148l).

The match is perfect: binding site of 12 atoms

alignment length: 12 matches

rmsd between binding site and matched atom on the protein target: 0

coverage: 1 (100% of binding site atoms matched to a protein surface atom)

mean distortion: 0

pepit score (weight of the matching of the correspondance): 144 (12**2)
