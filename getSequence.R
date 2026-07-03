library(bio3d)

args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  cat("usage: getPeptide.R pdbfile [chains]\n")
  q()
}

pdbfile = args[1]
ext = tools::file_ext(pdbfile)

chains = NULL
if (length(args)>1) {
   chains = args[2]
   chains = unlist(strsplit(chains,split=","))
}

if (ext == "") pdb = read.pdb(pdbfile, rm.alt=TRUE)
if (ext == "pdb") pdb = read.pdb(pdbfile, rm.alt=TRUE)
if (ext == "cif") pdb = suppressWarnings(read.cif(pdbfile, rm.alt=TRUE))


if (is.null(chains))  {
   chains = unique(pdb$atom$chain)
}

for (chain in chains) {
    pdb1 = trim(pdb, string="calpha", chain=chain)
    cat (">",pdbfile, chain, "\n", append=TRUE)
    cat (aa321(pdb1$atom$resid), "\n", sep="")
    #print(pdb1, printseq=TRUE)
}
