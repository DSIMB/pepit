# output
The script creates a .dat file and a .pdb file for each line in the <Bank> directory.
Each .dat file contains the data of the receptor atoms in contact with the peptide (by default at less than 5A).
The colums of the dataframe are:
- eleno : Atom serial number
- elety : Atom name
- resid : Residue name
- chain : Chain identifier
- resno : Residue sequence number
- insert : Code for insertions of residues
- x, y, z : Atom coordinates
- pc : phycico-chemical class of residues
  - 1 Nonpolar Hydrophilic (Ala, Phe, Ile,Met, Leu, Pro, Val)
  - 2 Nonpolar Hydrophobic (Gly, Trp)
  - 3 Polar-uncharged  Hydrophilic (Cys)
  - 4 Polar-uncharged  Hydrophobic (Asn, Gln, Ser, Thr, Tyr)
  - 5 Negatively-charged  Hydrophobic (Asp, Glu)
  - 6 Positively-charged Hydrophobic (Lys, His, Arg)

Atoms are annotated A,C,O,N,a,b,c,o,n and placed in column named elety.

The file is named idX:Y.dat, with X the receptor chain or list of chains and Y the peptide chain.

The ligands/peptides are contained in separate .pdb files named idY.pdb
(id is the pdb identifier of the receptor protein)

# bash
echo "PDB Receptor.Chain Peptide.Chain Peptide.Size" > bank.dat
echo "148l E S" >> bank.dat
echo "1a07 A C" >> bank.dat
echo "1a07 B D" >> bank.dat
echo "1a08 A C" >> bank.dat

mkdir Bank
Rscript ./BuiltBSBank.R  bank.dat Bank

make files: Bank/148lE:S.dat, Bank/1a07A:C.dat, Bank/1a07B:D.dat, Bank/1a08A:C.dat
and
Bank/148lS.pdb, Bank/1a07C.pdb, Bank/1a07D.pdb, Bank/1a08C.pdb

