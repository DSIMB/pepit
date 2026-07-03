library(bio3d)
library(pepit)

set.pepit("NBHITS", 10)

args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  cat("usage: getPeptide.R prefix\n")
  q()
}

prefix = args[1] # a prefix for generating  output files
tch = args[2] # a prefix for generating  output files

scorefile = paste(prefix, ".score", sep="")
pepfile = paste(prefix, ".pep", sep="")
#
# sort results for best scores and or p-value (scores here)
#
D = read.table(scorefile, header=TRUE)
if (nrow(D) == 0) {
  message("no hit found")
  q()
}

#
# sorts output of pose files
#

noclash = D$clashes<=get.pepit("MAXCLASHES")
noclash[is.na(noclash)]=TRUE
D = D[noclash,]
o = order(D$alen, decreasing=TRUE)# 
D = D[o, ]
nbhits = min(get.pepit("NBHITS"), nrow(D))
if (nbhits > 0) {
  D = D[1:nbhits,]
}

bsfiles = as.character(D$bs)
bsfiles = basename(bsfiles)
bsfiles = tools::file_path_sans_ext(bsfiles)
ids = substring(bsfiles, 1, 4)
chains = substring(bsfiles, 5)
chains = strsplit(chains, split=":")
pchains = unlist(lapply(chains, "[", 1))
pchains = unlist(lapply(chains, "[", 2))

pdb = read.pdb(prefix)
tch = unlist(strsplit(tch,split=","))
target = trim(pdb, chain=tch, string="calpha")
cat (">", bsfiles[i],"\n", file=pepfile, append=TRUE)
cat (aa321(target$atom$resid), "\n", sep="", file=pepfile, append=TRUE)

for (i in 1:length(bsfiles)) {
    pdb = read.pdb(ids[i])
    pep = trim(pdb, chain=pchains[i], string="calpha")
    cat (">", bsfiles[i],"\n")
    cat (aa321(pep$atom$resid), "\n", sep="")
    cat (">", bsfiles[i],"\n", file=pepfile, append=TRUE)
    cat (aa321(pep$atom$resid), "\n", sep="", file=pepfile, append=TRUE)
}
