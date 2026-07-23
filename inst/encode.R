#!/usr/bin/Rscript
# 
# encode a pdb file to a data file used as input for pepit
# inputs:
# a pdb file or a pdb id
# list of chains to be encoded
# 
# output:
# a data file containing a table with columns
# eleno elety resid chain resno insert x y z hseu hsed

args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  cat("usage: encode.R id.pdb (or id.cif)\n")
  q()
}
pdbfile = args[1] # file or id
tchain = args[2] # a chain or a chain list A,B,C or * for all chains
#
# user defined parameters
#
library(pepit)

pdb = bio3d::read.cif(pdbfile)

chainlist = unique(pdb$atom$chain)
tchain = unlist(strsplit(tchain,split="")) # could be ",", but works fine as is
if (tchain[1]=="*") tchain = chainlist
tchain = intersect(tchain, chainlist)

pdb = bio3d::trim.pdb(pdb, string = "protein")
pdb = bio3d::trim.pdb(pdb, chain=tchain)

D = encode(pdb)

outfile = paste(tools::file_path_sans_ext(pdbfile),tchain,".dat",sep="")
write.table(D, quote=FALSE, row.name=FALSE, file=outfile)
