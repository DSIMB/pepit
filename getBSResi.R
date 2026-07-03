#!/usr/bin/Rscript
# Rscript getBSResi.R pdb tchain pchain
# output : a resid file
library(bio3d)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  cat("usage: getBSResi.R pdbfile target_chain [pep_chain]\n")
  q()
}

pdbfile = args[1]
tchain = args[2]

pchain = NULL
if (length(args) >= 3) pchain = args[3]

library(pepit)

read.config("pepit.cfg")

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

chains=unique(pdb$atom$chain)
tchain = unlist(strsplit(tchain,split=","))
tchain = intersect(tchain, chains)
if (length(tchain)==0) {
    message(paste("chain", args[2], " not found"))
    q()
}
target.pdb = bio3d::trim.pdb(pdb, chain=tchain, string="protein")

target.data = encode(target.pdb)

if (!is.null(pchain)) {
   lchain=as.character(pchain)
   lchain=substring(lchain,1,1) #only one ligand chain
   if (get.pepit("PROTEIN")) { # for peptides (!=sugars)
      ligand.pdb = bio3d::trim.pdb(pdb, chain=lchain, string="protein")
   } else {
      ligand.pdb = bio3d::trim.pdb(pdb, chain=lchain)
   }
  inds = get_binding_sites(target.pdb, ligand.pdb, cutoff=get.pepit("CONTACT"), add=get.pepit("ADD"))
  eleno = target.pdb$atom$eleno[inds$atom]
  target.data = target.data[target.data$eleno%in%eleno,]
}

cat("#", pdbfile, tchain, "\n")
ind = !duplicated(target.data$resno)
target.data = target.data[ind,]
ind = which(target.data$hseu <= get.pepit("HSECUTOFF"))
cat(paste(target.data$chain[ind], target.data$resno[ind], sep=""),sep=",")
cat("\n")



