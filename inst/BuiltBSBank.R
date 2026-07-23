#!/usr/bin/Rscript

# 
# Create a bank of binding sites
# inputs:
# a data file containing a table: pdbid [pdb_file] target_chains peptide_chain with header
# header = id [pdb] tchain pchain
# a bank directory for the bs pdb file
# 
# output:
# a set of encode .dat files containing a binding site (chains not renamed) and
# .pdb files containing peptides ligand (chain renamed P). A dat file and a pdb
# file with the same id form a complex

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  cat("usage: BuiltBSBank.R infile.dat Bankdirectory resid\n")
  q()
}
infile = args[1]
BSBANK = args[2]
residout = FALSE

if (length(args) == 3 && args[3] == "resid") {
   residout = TRUE
}

library(pepit)

read.config("pepit.cfg")

D = read.table(infile, colClasses="character",header=TRUE, tryLogical=FALSE)
names(D) = c("id", "tchain", "pchain")

for (i in 1:nrow(D)) {
  id = as.character(D$id[i])
  pdbfile = D$pdb[i]
  if (is.null(D$pdb[i])) pdbfile = ""
  pdb = NULL
  if (!file.exists(pdbfile)) 
     cat("pdb file", pdbfile, "not found\n")
  if (file.exists(pdbfile)) {
     pdb = bio3d::read.pdb(pdbfile)
     #id = basename(tools::file_path_sans_ext(id))
  }
  if (is.null(pdb)) pdb = bio3d::read.cif(id)
  if (is.null(pdb)) {
    message(paste("target",id,"not found"))
    next
  } 

  tchain=as.character(D$tchain[i])# a target chain can be a list of chains (e.g. H,L)
  lchain=as.character(D$pchain[i])
  lchain=substring(lchain,1,1) #only one ligand chain
  outfile=paste(BSBANK,"/",id,tchain,":",lchain,".dat",sep="")
  residfile=paste(BSBANK,"/",id, tchain,":",lchain,".resid",sep="")
  
  #pdb=bio3d::trim.pdb(pdb, string="protein")
  pdb=bio3d::trim.pdb(pdb, string="noh")
  chains=unique(pdb$atom$chain)
  tchain = unlist(strsplit(tchain,split=","))
  tchain=intersect(tchain, chains)
  if (length(tchain)==0) {
    next
  }
  
  target.pdb = bio3d::trim.pdb(pdb, chain=tchain, string="protein")
  if (get.pepit("PROTEIN")) { # for peptides (!=sugars)
    ligand.pdb = bio3d::trim.pdb(pdb, chain=lchain, string="protein")
  } else {
    ligand.pdb = bio3d::trim.pdb(pdb, chain=lchain)
  }
  
  inds = get_binding_sites(target.pdb, ligand.pdb, cutoff=get.pepit("CONTACT"), add=get.pepit("ADD"))
  unlink(outfile)
  unlink(residfile)
  for (ch in tchain) {
    pdb=bio3d::trim.pdb(target.pdb, chain=ch)
    target.data = encode(pdb)
    eleno = pdb$atom$eleno[inds$atom]
    target.data = target.data[target.data$eleno%in%eleno,]
    col = ifelse(file.exists(outfile), FALSE, TRUE)
    cat("creation of ", outfile,"\n")
    write.table(target.data, quote=FALSE, col.names = col, row.names=FALSE, file=outfile, append=TRUE) #
    if (0) {
        pdbtargetfile = paste(BSBANK, "/", id, tchain, ".pdb",sep="")
	target.pdb = bio3d::trim.pdb(target.pdb, chain=ch, eleno=eleno)
        bio3d::write.pdb(target.pdb, file=pdbtargetfile)
    }
  }
  if (residout) {
       target.data = read.table(outfile, header=TRUE, tryLogical=FALSE) # chain F interpreted as FALSE otherwise
       cat("#", pdbfile, tchain, "\n")
       ind = !duplicated(target.data$resno)
       target.data = target.data[ind,]
       ind = which(target.data$hseu <= get.pepit("HSECUTOFF"))
       cat(paste(target.data$chain[ind], target.data$resno[ind], sep=""),sep=",", file=residfile)
  }
  id = substring(id,1,4)
  outfile=paste(tools::file_path_sans_ext(outfile),".pdb",sep="")
  bio3d::write.pdb(ligand.pdb, file=outfile)
}


