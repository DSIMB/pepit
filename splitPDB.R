# split a pdb into pdbs with one chain 
library(bio3d)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  cat("usage: splitPDB pdbfile [chain_list [dest_chain_list]]\n")
  q()
}

pdbfile = args[1]
ext = tools::file_ext(pdbfile)
pdbname = basename(tools::file_path_sans_ext(pdbfile))
if (ext == "") pdb = read.pdb(pdbfile)
if (ext == "pdb") pdb = read.pdb(pdbfile)
if (ext == "cif") pdb = suppressWarnings(read.cif(pdbfile))

chains = unique(pdb$atom$chain)

if (length(args) == 2) {
   charg = args[2]
   charg = unlist(strsplit(charg, split=","))
   if (! charg%in%chains) {
      cat("wrong chains", charg, "\n")
      exit()
   }
   chains = charg
}

dest_chains = NULL
if (length(args) == 3) {
   dest_chains = args[3]
   dest_chains = unlist(strsplit(dest_chains, split=","))
   if (length(dest_chains) != length(chains)) {
      cat("number of destination chains != number of chains\n")
      exit()
   }
}

for (k in 1:length(chains)) {
    chain = chains[k]
    pdb1 = trim.pdb(pdb, chain=chains[k])
    if (!is.null(dest_chains)) {
       chain = dest_chains[k]
       pdb1$atom$chain = chain
    }
    outfile = paste(pdbname,"_",chain,".pdb", sep="")
    write.pdb(pdb1, file=outfile)
}