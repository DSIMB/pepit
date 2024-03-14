#!/usr/bin/Rscript

# 
# Create a bank of binding sites
# inputs:
# a data file containing a table: pdbid target_chains peptide_chain
# a bank directory for the bs pdb file
# 
# output:
# a set of encode .dat files containing a binding site (chains not renamed) and
# .pdb files containing peptides ligand (chain renamed P). A dat file and a pdb
# file with the same id form a complex

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  cat("usage: BuiltBSBank.R infile.dat Bankdirectory\n")
  q()
}
infile=args[1]
BSBANK=args[2]
 
library(pepit)

#
# user defined parameters
#
set.pepit("CONTACT", 8.0)
set.pepit("ADD", "calpha")
set.pepit("HSECUTOFF", 13)

D=read.table(infile, header=TRUE)

for (i in 1:nrow(D)) {
  id=as.character(D[i,1])
  pdb=bio3d::read.pdb(id)
  tchain=as.character(D[i,2])# a target chain can be a list of chains (e.g. H,L)
  lchain=as.character(D[i,3])
  lchain=substring(lchain,1,1) #only one ligand chain
  outfile=paste(BSBANK,"/",id,tchain,":",lchain,".dat",sep="")

  chainlist=unique(pdb$atom$chain)
  tchain=unlist(strsplit(tchain,split=""))
  tchain=intersect(tchain, chainlist)
  tchain=setdiff(tchain, lchain)

  pdb=bio3d::trim.pdb(pdb, string="protein")
  pdb=bio3d::trim.pdb(pdb, string="noh")
  
  inds = get_binding_sites_2(pdb, target_chains=tchain, ligand_chains=lchain, add=get.pepit("ADD"))
  
  for (ch in tchain) {
    target.pdb=bio3d::trim.pdb(pdb, chain=ch)
    target.data = encode(target.pdb)
    eleno = pdb$atom$eleno[inds$atom]
    target.data = target.data[target.data$eleno%in%eleno,]
    col = ifelse(file.exists(outfile), FALSE, TRUE)
    write.table(target.data, quote=FALSE, col.names = col, row.names=FALSE, file=outfile) #
  }
  id = substring(id,1,4)
  ligand.pdb = bio3d::trim.pdb(pdb, chain=lchain)
  outfile=paste(BSBANK,"/",id,lchain,".pdb",sep="")
  bio3d::write.pdb(ligand.pdb, file=outfile)
}

