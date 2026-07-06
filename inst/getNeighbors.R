#!/usr/bin/Rscript
# Rscript getNeighbors.R pdb  residfile
# output : a resid file
library(bio3d)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  cat("usage: getNeighbors.R pdbfile tchain cutoff residues\n")
  q()
}

pdbfile = args[1]
tchain = args[2]
cutoff = args[3]
cutoff = max(0, as.integer(cutoff))
residues = args[4]

library(pepit)

ext = tools::file_ext(pdbfile)
id = tools::file_path_sans_ext(pdbfile)
if (ext == "" | ext == "pdb" | ext == "cif") {
  pdb = NULL
  #if (ext == "" | ext == "pdb") pdb = read.pdb(pdbfile)
  if (ext == "pdb") pdb = read.pdb(pdbfile)
  if (ext == "") {
     pdbfile = get.pdb(pdbfile, format="cif")
     ext = "cif" # patch for error of bio3d
  }
  if(ext == "" | ext == "cif") {
     pdb = read.cif(pdbfile)
  }
  if (is.null(pdb)) {
    message("target file not found")
    q()
  }
}
pdb = bio3d::trim.pdb(pdb, string="noh")

tchain=as.character(tchain)# a target chain can be a list of chains (e.g. H,L)
# decode resi=[chain][resno][insert]->chi,resnoi,insi (i meaning input)
resi = unlist(strsplit(residues,split=","))
nc = nchar(resi)
chi = substring(resi, 1, 1)
resnoi =  as.integer(substring(resi, 2, nc))
insi = rep(NA,length(resi))
insi[is.na(resnoi)] = substring(resi, nc, nc)[is.na(resnoi)]
resnoi[is.na(resnoi)] =  as.integer(substring(resi, 2, nc))[is.na(resnoi)]

chains=unique(pdb$atom$chain)
tchain = unlist(strsplit(tchain,split=","))
tchain=intersect(tchain, chains)
if (length(tchain)==0) {
    message(paste("chain", args[2], " not found"))
    q()
}
  
target.pdb = bio3d::trim.pdb(pdb, chain=tchain, string="protein")

ligand.pdb = bio3d::trim.pdb(pdb, chain=chi, resno=resnoi, string="protein")
  
inds = get_binding_sites(target.pdb, ligand.pdb, cutoff=cutoff, add="calpha")

resno = target.pdb$atom$resno[inds$atom]
chains = target.pdb$atom$chain[inds$atom]
keep = which(!duplicated(resno))

cat("#", pdbfile, tchain, "\n")
cat(paste(chains[keep], resno[keep], sep=""),sep=",")
cat("\n")
