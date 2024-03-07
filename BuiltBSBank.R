#!/usr/bin/Rscript

# 
# Create a bank of binding sites
# inputs:
# a data file containing a table: pdbid target_chains peptide_chain
# a bank directory for the bs pdb file
# 
# output:
# a set of pdb file containing a binding site (chains not renamed) and a peptide ligand (chain renamed P)

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
set.pepit("ACC", 10)
set.pepit("radius", 1.4)

D=read.table(infile, header=TRUE)

for (i in 1:nrow(D)) {
  id=as.character(D[i,1])
  pdb=bio3d::read.pdb(id)
  tchain=as.character(D[i,2])# a target chain is can be a list of chains (e.g. H,L)
  lchain=as.character(D[i,3])
  lchain=substring(lchain,1,1) #only one ligand chain
  outfile=paste(BSBANK,"/",id,tchain,":",lchain,".pdb",sep="")
  
  chainlist=unique(pdb$atom$chain)
  tchain=unlist(strsplit(tchain,split=""))
  tchain=intersect(tchain, chainlist)
  tchain=setdiff(tchain, lchain)

  pdb=bio3d::trim.pdb(pdb, string="protein")
  pdb=bio3d::trim.pdb(pdb, string="noh")
  pdb=type_atoms(pdb, ligchain=lchain) # ligchain is not typed
  holo=access_atoms(pdb, chains=tchain, add=get.pepit("ADD"), minacc=get.pepit("ACC"), select=FALSE)
  bsnlig=get_binding_sites(holo, target_chains=tchain, ligand_chains=lchain, add=get.pepit("ADD"), add.ligand=TRUE, access=get.pepit("ACC"))
  #bsnlig$atom$chain[bsnlig$atom$chain==lchain]="P"
  if (!is.null(bsnlig)) bio3d::write.pdb(bsnlig, file=outfile)
}