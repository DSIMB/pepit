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
set.pepit("PROTEIN", TRUE)
set.pepit("CONTACT", 6.0)
set.pepit("ADD", "calpha")
set.pepit("HSECUTOFF", 13)

D=read.table(infile, header=TRUE)

resid.mode = FALSE
if (ncol(D) == 4) resid.mode = TRUE

for (i in 1:nrow(D)) {
  id=as.character(D[i,1])
  pdb=bio3d::read.pdb(id)
  tchain=as.character(D[i,2])# a target chain can be a list of chains (e.g. H,L)
  lchain=as.character(D[i,3])
  lchain=substring(lchain,1,1) #only one ligand chain
  outfile=paste(BSBANK,"/",id,tchain,":",lchain,".dat",sep="")

  #pdb=bio3d::trim.pdb(pdb, string="protein")
  pdb=bio3d::trim.pdb(pdb, string="noh")
  chains=unique(pdb$atom$chain)
  tchain=intersect(tchain, chains)
  if (length(tchain)==0) {
    next
  }
  
  target.pdb = bio3d::trim.pdb(pdb, chain=tchain, string="protein")
  res = NULL
  if (resid.mode) {
    res = as.integer(unlist(strsplit(D[i,4], split=",")))
    outfile=paste(BSBANK,"/",id,tchain,":",lchain, res[1],".dat",sep="")
  }
  if (get.pepit("PROTEIN")) {
    ligand.pdb = bio3d::trim.pdb(pdb, chain=lchain, resno=res, string="protein")
  } else {
    ligand.pdb = bio3d::trim.pdb(pdb, chain=lchain, resno=res)
  }
  
  inds = get_binding_sites(target.pdb, ligand.pdb, add=get.pepit("ADD"))
  
  for (ch in tchain) {
    pdb=bio3d::trim.pdb(target.pdb, chain=ch)
    target.data = encode(pdb)
    eleno = pdb$atom$eleno[inds$atom]
    target.data = target.data[target.data$eleno%in%eleno,]
    col = ifelse(file.exists(outfile), FALSE, TRUE)
    write.table(target.data, quote=FALSE, col.names = col, row.names=FALSE, file=outfile, append=TRUE) #
  }
  id = substring(id,1,4)
  outfile=paste(tools::file_path_sans_ext(outfile),".pdb",sep="")
  bio3d::write.pdb(ligand.pdb, file=outfile)
}

