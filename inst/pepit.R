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
set.pepit("CONTACT",10.0)
set.pepit("PRECISION", 1.0)
set.pepit("ADD","calpha")
set.pepit("MINCLIQUE",6)
set.pepit("MINSCORE", 0) # tous les scores
set.pepit("MAXCLASHES", 0)
set.pepit("NBCLIQUES", 5)
set.pepit("NBHITS", 100)
set.pepit("POSE", FALSE)
set.pepit("TYPES", c("A","C","O","N"))
set.pepit("PVALUE", TRUE)
#
criteria = "score"
#
#

args = commandArgs(trailingOnly=TRUE)

if (length(args)<4) {
  cat("usage: pepit.R target target_chain bs_bank out_prefix\n")
  q()
}

tfile = args[1] # a pdb id or a pdb file (can be a .gz file)
tchain = args[2] # a chain or a chain list A,B,C or * for all chains
bank = args[3] # a directory with binding sites or a pdb file
prefix = args[4] # a prefix for generating  output files

deltadist = get.pepit("PRECISION")

scorefile = paste(prefix, ".score", sep="")
alignfile = paste(prefix, ".al", sep="")

pdb = read.pdb(tfile)
###
resi=get.pepit("RESIDUES")
resi=unlist(strsplit(resi,split=","))
resi=as.integer(resi)
###

chainlist = unique(pdb$atom$chain)
tchain = unlist(strsplit(tchain,split=""))
if (tchain[1]=="*") tchain = chainlist
tchain = intersect(tchain, chainlist)

pdb = bio3d::trim.pdb(pdb,chain=tchain)

if(nrow(pdb$atom)==0) {
  message("unknown target")
  q()
}

pdb = type_atoms(pdb)
target = access_atoms(pdb, chains=tchain, add=get.pepit("ADD"), select=FALSE)
target = trim_non_acc_atoms(target,minacc=get.pepit("ACC"))
if (length(resi)>0) target=bio3d::trim.pdb(target, resno=resi)

X = target$atom[,c("x","y","z")]
X = as.matrix(X)
XProp = target$atom$elesy
N = nrow(X)

if (file.exists(bank) & dir.exists(bank)) { # if bank is directory of bs files
  bslist = dir(bank)
  bslist = paste(bank,"/",bslist,sep="")
} else if (file.exists(bank)) { # if bank is a bs file
  bslist = bank
} else {
  message("unknown bank")
  q()
}


if (!file.exists(scorefile)) cat("index bs target precision bslen alen rmsd coverage meandist score\n", file=scorefile)
count = 0
for (bsfile in bslist) {
  cat("bsfile =", bsfile,"\n")
  bsnlig = read.pdb(bsfile)
  chains = unique(bsnlig$atom$chain)
  pos = min(gregexpr(":", bsfile)[[1]])
  pepchain = substring(bsfile, pos+1, pos+1)
  chains = setdiff(chains, pepchain)
  bs = bio3d::trim.pdb(bsnlig, chain=chains) # to get rid of peptide (chain P)
  
  Y = bs$atom[,c("x","y","z")]
  Y = as.matrix(Y)
  YProp = bs$atom$elesy
  
  result = cliques(X, XProp, Y, YProp, deltadist=1.0, mindist=0.0, maxdist=get.pepit("MAXDELTADIST"), types=get.pepit("TYPES"))
  if (length(result)==0) {
    cat("no clique\n")
    count = count+1
    cat(count, bsfile, tfile, deltadist, nrow(Y), 0, 100, 0, maxdist+1, 0, "\n", file=scorefile, append=TRUE)
    next
  }
  clusters = extend_cliques(X, XProp, Y, YProp, result, deltadist=deltadist)
  clusters = remove_redundant_clusters (clusters)
  
  scores = score_clusters(clusters, X, Y, deltadist=deltadist)
  # test
  scores$score = sqrt(pmax(0,scores$score))
  scores$score = scores$score/scores$len
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
          cat(count, bsfile, tfile, deltadist, scores$len[i], scores$alen[i], scores$rmsd[i], scores$coverage[i], scores$distorsion[i], scores$score[i], "\n", file=scorefile, append=TRUE)
          if (get.pepit("POSE")) {
	     # patch 
	     bs$atom$insert[is.na(bs$atom$insert)]=""
	     target$atom$insert[is.na(target$atom$insert)]=""
	     # fin patch 
             # output matching (= non seq. alignment) in .al file
             bs_res = paste(bs$atom$resno[J],":",bs$atom$insert[J],":",bs$atom$elety[J],":",bs$atom$chain[J], sep="")
             target_res = paste(target$atom$resno[I],":",target$atom$insert[I],":",target$atom$elety[I],":",target$atom$chain[I], sep="")
             cat(">",count, bsfile, tfile, deltadist, scores$len[i], scores$alen[i], scores$rmsd[i], scores$coverage[i], scores$distorsion[i], scores$score[i], scores$pval[i],"\n", file=alignfile, append=TRUE)
             cat(bs_res, "\n", file=alignfile, append=TRUE)
             cat(target_res, "\n", file=alignfile, append=TRUE)
             # output binding site posed on target in a .pdb file 
             ligand.moved = superpose_sites(clusters[[i]], bsnlig, target, pepchain=pepchain)
             nbclashes = clashes(ligand.moved, pdb)
             cat("nbclashes=", nbclashes,"\n")
             if (nbclashes <= get.pepit("MAXCLASHES")) {
              	pepfile = paste(prefix, "-", count, ".pdb",sep="")
              	write.pdb(pdb=ligand.moved, file = pepfile)
             }
          }
        }
      }
  }

#
# sorts output of pose files
#
  if (get.pepit("POSE")) {
     indexfile = paste(prefix, ".index", sep="")
     D = read.table(scorefile, header=TRUE)
     if (nrow(D)>0) {
     	o = order(D$score, decreasing=TRUE)
     	D = D[o, ]
     	D[,1] = 1:nrow(D)
     	nbhits = min(get.pepit("NBHITS"), nrow(D))
     	D = D[1:nbhits,]
     	write.table( D, quote=FALSE, row=FALSE, file=scorefile)
     }

     count = 1
     for (k in 1:nrow(D)) {
     	 pdbfile = paste(prefix, "-", D[k,"index"], ".pdb", sep="")
    	 if (!file.exists(pdbfile)) next
    	 outfile = paste(prefix, "_sorted","-", count, ".pdb", sep="")
    	 command = paste("mv ", pdbfile, " ", outfile, sep="")
    	 print(command)
    	 system(command)
	     cat (count,  k, D$bs[k], D$bslen[k], D$alen[k], D$rmsd[k], D$score[k], "\n", file=indexfile, append=TRUE)
	     count = count + 1
     }
   }