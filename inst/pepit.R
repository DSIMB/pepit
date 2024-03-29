# 
# search against binding site bank
#
# inputs:
# target pdb file
# target chain(s)
# binding site bank or one bs pdb file
# label for generating output files
# outputs:
# score file
# alignment/matching file
# pdb files of posed peptides 
library(pepit)
library(bio3d)

#
# user defined parameters
#
set.pepit("RESIDUES","")
set.pepit("PRECISION", 1)
set.pepit("ADD","calpha")
set.pepit("MINCLIQUE",6)
set.pepit("MINSCORE", 0) # tous les scores
set.pepit("MAXCLASHES", 10) # clashes allowed
set.pepit("NBCLIQUES", 1)
set.pepit("NBHITS", 100)
set.pepit("POSE", TRUE)
set.pepit("TYPES", c("A","C","O","N","c","o","b","n","a")) # all atoms
#set.pepit("TYPES", c("A","C","O","N")) # backbone atoms
set.pepit("PVALUE", TRUE)
set.pepit("MODE",1)
set.pepit("HSE", 10)
#
criteria = "score"
#
#

args = commandArgs(trailingOnly=TRUE)

if (length(args)<4) {
  cat("usage: pepit.R target chain bs_bank out_prefix\n")
  q()
}

tfile = args[1] # a pdb id or a pdb file (can be a .gz file)
tchain = args[2] # a chain or a chain list A,B,C or * for all chains
bank = args[3] # a directory with binding sites or a pdb file
prefix = args[4] # a prefix for generating  output files

deltadist = get.pepit("PRECISION")
maxdist = get.pepit("MAXDELTADIST")
allscorefile = paste(prefix, "all.score", sep="")
scorefile = paste(prefix, ".score", sep="")
alignfile = paste(prefix, ".al", sep="")
residfile = paste(prefix, ".resi", sep="")

if (tools::file_ext(tfile)=="" | tools::file_ext(tfile)==".pdb" | tools::file_ext(tfile)==".cif") {
  #pdb = read.pdb(tfile)
  pdb = NULL
  pdb = read.pdb(tfile)
  if (is.null(pdb)) pdb = read.cif(tfile)
  if (is.null(pdb)) {
    message("target file not found")
    q()
  }
  pdb = bio3d::trim.pdb(pdb, string="protein")
  pdb = bio3d::trim.pdb(pdb, string="noh")

  chainlist = unique(pdb$atom$chain)
  tchain = unlist(strsplit(tchain,split=""))
  if (tchain[1]=="*") tchain = chainlist
  tchain = intersect(tchain, chainlist)

  pdb = bio3d::trim.pdb(pdb,chain=tchain)
  if(nrow(pdb$atom)==0) {
    message("wrong target chain(s)")
    q()
  }

  cat("encode target...\n")
  target.data = encode(pdb)
}

if (tools::file_ext(tfile)==".dat") {
  target.data = read.table(tfile, header=TRUE)
  chainlist = unique(target.data$chain)
  tchain = unlist(strsplit(tchain,split=""))
  if (tchain[1]=="*") tchain = chainlist
  tchain = intersect(tchain, chainlist)
  
  target.data = target.data[target.data$chain %in% tchain,]
  if(nrow(target.dat)==0) {
    message("wrong target chain(s)")
    q()
  }
}

###
unlink(scorefile)
unlink(alignfile)
unlink(residfile)
###
resi=get.pepit("RESIDUES")
resi=unlist(strsplit(resi,split=","))
resi=as.integer(resi)
###



if (file.exists(bank) & dir.exists(bank)) { # if bank is directory of bs files
  bslist = dir(bank, pattern=".dat")
  bslist = paste(bank,"/",bslist,sep="")
} else if (file.exists(bank)) { # if bank is a bs file
  bslist = bank
} else {
  message("unknown bank")
  q()
}


if (length(resi)>0) target.data = target.data[target.data$resid%in%resi,]

X = target.data[,c("x","y","z")]
X = as.matrix(X)
XProp = target.data
N = nrow(X)

cat("index bs target precision bslen alen rmsd coverage meandist score clashes\n", file=allscorefile)
count = 0
for (bsfile in bslist) {
  cat("bsfile =", bsfile,"\n")
  bs.data = read.table(bsfile, header=TRUE)
  chains = unique(bs.data$chain)

  Y = bs.data[,c("x","y","z")]
  Y = as.matrix(Y)
  YProp = bs.data
  
  result = cliques(X, XProp, Y, YProp, deltadist=1.0, mindist=0.0, maxdist=get.pepit("MAXDELTADIST"), types=get.pepit("TYPES"))
  if (length(result)==0) {
    cat("no clique\n")
    #count = count+1
    #cat(count, bsfile, tfile, deltadist, nrow(Y), 0, 100, 0, maxdist+1, 0, 0, "\n", file=allscorefile, append=TRUE)
    next
  }
  #clusters = result
  clusters = extend_cliques(X, XProp, Y, YProp, result, deltadist=get.pepit("PRECISION"))
  clusters = remove_redundant_clusters (clusters)
  
  scores = score_clusters(clusters, X, Y, deltadist=deltadist)
  # test
  # scores$score = sqrt(pmax(0,scores$score))
  # scores$score = scores$score/scores$len
  # fin test
  o = order(scores[[criteria]], decreasing=TRUE)
  nbcliques = min(length(clusters), get.pepit("NBCLIQUES"))
  o = o[1:nbcliques]
  scores = lapply(scores, "[", o)
  clusters = clusters[o]

  for (i in 1:length(clusters)) {
      if(scores$score[i] >= get.pepit("MINSCORE")) {
          I = (clusters[[i]]-1)%%N+1 # target indices
          J = (clusters[[i]]-1)%/%N+1# bs indices
 	        count = count + 1
          # patch 
          bs.data$insert[is.na(bs.data$insert)]=""
          target.data$insert[is.na(target.data$insert)]=""
          # fin patch 
          # output matching (= non seq. alignment) in .al file
          bs_res = paste(bs.data$resno[J],":",bs.data$insert[J],":",bs.data$elety[J],":",bs.data$chain[J], sep="")
          target_res = paste(target.data$resno[I],":",target.data$insert[I],":",target.data$elety[I],":",target.data$chain[I], sep="")
          cat(">",count, bsfile, tfile, deltadist, scores$len[i], scores$alen[i], scores$rmsd[i], scores$coverage[i], scores$distorsion[i], scores$score[i], scores$pval[i],"\n", file=alignfile, append=TRUE)
          cat(bs_res, "\n", file=alignfile, append=TRUE)
          cat(target_res, "\n", file=alignfile, append=TRUE)
          # output of BS residues
          resno = target.data$resno[I]
          ins = target.data$insert[I]
          ins[is.na(ins)] = ""
          bsres = unique(paste(resno, ins, sep=""))
          cat(sort(bsres), "\n", file=residfile, append=TRUE)
          nbclashes = 0
          if (get.pepit("POSE")) {
             pos = min(gregexpr(":", bsfile)[[1]])
             pepchain = substring(bsfile, pos+1, pos+1)
             #bsid = substring(basename(bsfile),1,4)
             bsid = tools::file_path_sans_ext(bsfile)
             pepfile = paste(bsid, ".pdb", sep="")
             cat("pepfile=", pepfile, "\n")
             if (file.exists(pepfile)) {
                peptide = bio3d::read.pdb(pepfile)
             # output binding site posed on target in a .pdb file 
                ligand.moved = superpose_sites2(clusters[[i]], bs.data, target.data, peptide)
                nbclashes = clashes(ligand.moved, pdb)
                cat("nbclashes=", nbclashes,"\n")
                if (nbclashes <= get.pepit("MAXCLASHES")) {
              	  pepfile = paste(prefix, "-", count, ".pdb",sep="")
              	  write.pdb(pdb=ligand.moved, file = pepfile)
                }
             }
          }
          cat(count, bsfile, tfile, deltadist, scores$len[i], scores$alen[i], scores$rmsd[i], scores$coverage[i], scores$distorsion[i], scores$score[i], nbclashes, "\n", file=allscorefile, append=TRUE)
        }
  }
}
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
      command = paste("mv ", pepfile, " ", newfile, sep="")
      cat(command,"\n")
      system(command)
    }
    unlink(paste(prefix,"-*.pdb"))
  }
  
  D[,1] = 1:nrow(D)
  write.table(D, quote=FALSE, row=FALSE, file=scorefile)
  
}

