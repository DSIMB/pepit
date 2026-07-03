#!/usr/bin/Rscript

# 
# Create a bank of complexes 
# inputs:
# a data file containing a table: id target_chains peptide_chain [pdb_file] with header
# header = id tchain pchain [pdb]
# a bank directory for the bs pdb file
# 
# output:
# a set of pdb files with target and peptide chains written to Bankdirectory

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  cat("usage: getBank.R infile.dat Bankdirectory\n")
  q()
}
infile = args[1]
BSBANK = args[2]

library(pepit)

D = read.table(infile, colClasses="character", header=TRUE, tryLogical=FALSE)

for (i in 1:nrow(D)) {
  id = as.character(D$id[i])
  pdbfile = D$pdb[i]
  if (is.null(D$pdb[i])) pdbfile = ""
  pdb = NULL
  if (!file.exists(pdbfile)) 
     cat("pdb file", pdbfile, "not found\n")
  if (file.exists(pdbfile)) {
     pdb = bio3d::read.pdb(pdbfile)
     id = basename(tools::file_path_sans_ext(id))
  }
  if (is.null(pdb)) {
     cat("get cif file", id, "\n")
     pdb = bio3d::read.cif(id)
  }
  if (is.null(pdb)) {
    message(paste("target",id,"not found"))
    next
  } 

  tchain=as.character(D$tchain[i])# a target chain can be a list of chains (e.g. H,L)
  lchain=as.character(D$pchain[i])
  lchain=substring(lchain,1,1) #only one ligand chain
  outfile=paste(BSBANK,"/",id,tchain,":",lchain,".pdb",sep="")
  
  #pdb=bio3d::trim.pdb(pdb, string="protein")
  pdb=bio3d::trim.pdb(pdb, string="noh")
  chains=unique(pdb$atom$chain)
  tchain = unlist(strsplit(tchain,split=","))
  tchain=intersect(tchain, chains)
  if (length(tchain)==0) {
    next
  }
  
  target.pdb = bio3d::trim.pdb(pdb, chain=c(tchain, lchain) , string="protein")
  cat("write file:", outfile, "\n")
  bio3d::write.pdb(target.pdb, file=outfile)
}


