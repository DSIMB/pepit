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

read.config("pepit.cfg")

D = read.table(infile, colClasses="character",header=TRUE)

# resid.mode for bs given by a list of residue numbers (comma separated)
resid.mode = FALSE
if (ncol(D) == 4) resid.mode = TRUE

for (i in 1:nrow(D)) {
  id=as.character(D[i,1])
  pdb=NULL
  #result = bio3d::get.pdb(id)
  #if (file.exists(result)) pdb=bio3d::read.pdb(result)
  pdb=bio3d::read.pdb(id)
  if (is.null(pdb)) pdb=bio3d::read.cif(id)
  if (is.null(pdb)) {
    message(paste("target file",id,"not found"))
    next
  }
  #unlink(result)
  
  tchain=as.character(D[i,2])# a target chain can be a list of chains (e.g. H,L)
  lchain=as.character(D[i,3])
  lchain=substring(lchain,1,1) #only one ligand chain
  outfile=paste(BSBANK,"/",id,tchain,":",lchain,".dat",sep="")

  #pdb=bio3d::trim.pdb(pdb, string="protein")
  pdb=bio3d::trim.pdb(pdb, string="noh")
  chains=unique(pdb$atom$chain)
  tchain = unlist(strsplit(tchain,split=","))
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
  if (get.pepit("PROTEIN")) { # for peptides (!=sugars)
    ligand.pdb = bio3d::trim.pdb(pdb, chain=lchain, resno=res, string="protein")
  } else {
    ligand.pdb = bio3d::trim.pdb(pdb, chain=lchain, resno=res)
  }
  
  inds = get_binding_sites(target.pdb, ligand.pdb, add=get.pepit("ADD"))
  unlink(outfile)
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

