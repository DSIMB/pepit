library(pepit)

set.pepit("MAXCLASHES", 10) # clashes allowed
set.pepit("NBHITS", 100)

args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  cat("usage: selectPeptide.R prefix\n")
  q()
}

prefix = args[1] # a prefix for generating  output files

allscorefile = paste(prefix, "all.score", sep="")
scorefile = paste(prefix, ".score", sep="")
alignfile = paste(prefix, ".al", sep="")
residfile = paste(prefix, ".resi", sep="")

#
# sort results for best scores and or p-value (scores here)
#
D = read.table(allscorefile, header=TRUE)
if (nrow(D) == 0) {
  message("no hit found")
  q()
}

if (0) { # p_value needs a larger number of scores of negative bs 
  value = D$alen/D$precision
  o = order(value, decreasing=TRUE)
  value = value[o]
  ind = which(!duplicated(D$bs[o]))
  value = value[ind]
  pval = evd(value)
  D$pval = p_values(D$alen/D$precision, pval$loc, pval$scale)
}
#
# sorts output of pose files
#
noclash = D$clashes<=get.pepit("MAXCLASHES")
D = D[noclash,]
o = order(D$alen, decreasing=TRUE)# 
D = D[o, ]
nbhits = min(get.pepit("NBHITS"), nrow(D))

if (nbhits > 0) {
  D = D[1:nbhits,]
  
  Lresid = readLines(residfile)
  unlink(residfile)
  for (k in 1:nbhits) {
    i = D$index[k]
    cat(">", k, as.character(D[k,-1]), "\n", file=residfile, append=TRUE)
    cat(sort(Lresid[i]),"\n", file=residfile, append=TRUE)
  }
  
  Lalign = readLines(alignfile)
  unlink(alignfile)
  for (k in 1:nbhits) {
    i = D$index[k]-1
    cat(">", k,  as.character(D[k,-1]), "\n", file=alignfile, append=TRUE)
    cat(Lalign[(3*i+2):(3*i+3)], sep="\n", file=alignfile, append=TRUE)
  }
  
  if (get.pepit("POSE")) {
    for (k in 1:nbhits) {
      pepfile = paste(prefix, "-", D[k,"index"], ".pdb", sep="")
      if (!file.exists(pepfile)) next
      newfile = paste(prefix, "_peptide","-", k, ".pdb", sep="")
      command = paste("cp ", pepfile, " ", newfile, sep="")
      cat(command,"\n")
      system(command)
    }
  }
  
  D[,1] = 1:nrow(D)
  write.table(D, quote=FALSE, row=FALSE, file=scorefile)
  
}
