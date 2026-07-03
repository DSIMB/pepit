# search against a binding site bank
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
criteria = "alen"
#
#

args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("usage: pepit.R target chain bs_bank out_prefix [resid_file]\n")
}

tfile = args[1] # a pdb id or a pdb file (can be a .gz file)
tchain = args[2] # a chain or a chain list A,B,C or * for all chains
bank = args[3] # a directory with binding sites or a pdb file
prefix = args[4] # a prefix for generating  output files

#
read.config("pepit.cfg")
#

resi=c()
if (length(args)>4) {
   residfile = args[5]
   if (file.exists(residfile)) {
      resi = scan(residfile, what="", comment=">")
   } else {
      stop("no resid file\n")
   }
# decode resi=[chain][resno][insert]->chi,resnoi,insi (i meaning input)
  resi = unlist(strsplit(resi,split=","))
  nc = nchar(resi)
  chi = substring(resi, 1, 1)
  resnoi =  as.integer(substring(resi, 2, nc))
  insi = rep(NA,length(resi))
  insi[is.na(resnoi)] = substring(resi, nc, nc)[is.na(resnoi)]
  resnoi[is.na(resnoi)] =  as.integer(substring(resi, 2, nc))[is.na(resnoi)]
}
#
deltadist = get.pepit("PRECISION")
maxdist = get.pepit("MAXDELTADIST")
scorefile = paste(prefix, "scores.dat", sep="")
alignfile = paste(prefix, "alignments.dat", sep="")

ext = tools::file_ext(tfile)
id = tools::file_path_sans_ext(tfile)
if (ext == "" | ext == "pdb" | ext == "cif") {
  pdb = NULL
  #if (ext == "" | ext == "pdb") pdb = read.pdb(tfile)
  if (ext == "pdb") pdb = read.pdb(tfile)
  if (ext == "") {
     tfile = get.pdb(tfile, format="cif")
     ext = "cif" # patch for error of bio3d
  }
  if(ext == "" | ext == "cif") {
     pdb = read.cif(tfile)
  }
  if (is.null(pdb)) {
    message("target file not found")
    q()
  }
  pdb = bio3d::trim.pdb(pdb, string="protein")
  pdb = bio3d::trim.pdb(pdb, string="noh")

  chainlist = unique(pdb$atom$chain)
  tchain = unlist(strsplit(tchain,split=","))
  if (tchain[1]=="*") tchain = chainlist
  tchain = intersect(tchain, chainlist)
  
  pdb = bio3d::trim.pdb(pdb,chain=tchain)
  if(nrow(pdb$atom)==0) {
    message("wrong target chain(s)")
    q()
  }

  message("encode target...")
  target.data = encode(pdb)
  datfile = paste(tools::file_path_sans_ext(tfile),".dat", sep="")
  write.table(target.data, quote=FALSE, col.names = TRUE, row.names=FALSE, file=datfile, append=FALSE) # append ?
}

if (ext == "dat") {
  target.data = read.table(tfile, header=TRUE, tryLogical=FALSE)
  chainlist = unique(target.data$chain)
  tchain = unlist(strsplit(tchain,split=","))
  if (tchain[1]=="*") tchain = chainlist
  tchain = intersect(tchain, chainlist)
  
  target.data = target.data[target.data$chain %in% tchain,]
  if(nrow(target.data)==0) {
    message("wrong target chain(s)")
    q()
  }
}

###
unlink(scorefile)
unlink(alignfile)
###
###
if (file.exists(bank) & dir.exists(bank)) { # if bank is directory of bs files
  bslist = dir(bank, pattern="\\.dat")
  bslist = paste(bank,"/",bslist,sep="")
} else if (file.exists(bank)) { # if bank is a bs file
  bslist = bank
} else {
  message("unknown bank")
  q()
}
# filter with HSE = select Solvent Exposed atoms
hse = get.pepit("HSECUTOFF")
target.data = target.data[target.data$hseu<= hse | target.data$hsed <= hse,]
#
X = target.data[,c("x","y","z")]
X = as.matrix(X)
N = nrow(X)
XProp = target.data

Xbs = X
XbsProp = XProp
Nbs = N

if (length(resi)>0) {
   Ibs = c()
   for (k in 1:length(resi)) {
      Ibs = c(Ibs, which(target.data$resno%in%resnoi[k] & target.data$chain%in%chi[k] & target.data$insert%in%insi[k]))
   }
   Xbs = target.data[Ibs ,c("x","y","z")]
   XbsProp = target.data[Ibs,]
   Xbs = as.matrix(Xbs)
   Nbs = nrow(Xbs)
}
target.data$insert[is.na(target.data$insert)]=""
#
#
#
minalen = get.pepit("MINALEN")
nbhits = get.pepit("NBHITS")

allclusters = list()
allscores = data.frame()

for (bsfile in bslist) {
  message("bsfile =", bsfile)
  bs.data = read.table(bsfile, header=TRUE, tryLogical=FALSE)
  # filter with HSE
  bs.data = bs.data[bs.data$hseu <= hse | bs.data$hsed <= hse,]
  chains = unique(bs.data$chain)
  bs.data$insert[is.na(bs.data$insert)]=""
  Y = bs.data[,c("x","y","z")]
  Y = as.matrix(Y)
  YProp = bs.data
  
  clusters = cliques(Xbs, XbsProp, Y, YProp, deltadist=1.0, mindist=0.0, maxdist=get.pepit("MAXDELTADIST"), types=get.pepit("TYPES"))
  if (length(clusters)==0) {
    message("no matches")
    next
  }
  nbcliques = min(length(clusters), get.pepit("NBCLIQUES"))
  clusters = clusters[1:nbcliques] # clusters are sorted by the function cliques()
                    
  # renumbering of links in clusters to consider all target residues
  if (length(resi)>0) {
     for (k in 1:length(clusters)) {
     	 I = (clusters[[k]]-1)%%Nbs + 1
	 J = (clusters[[k]]-1)%/%Nbs + 1
	 I = Ibs[I]
 	 clusters[[k]] = (J-1)*N + I
     }
  }
  #
  if (get.pepit("EXTEND")) {
     message("extend cliques...", get.pepit("EXTEND"), length(clusters), "clusters\n", append=TRUE)
     clusters = extend_cliques(X, XProp, Y, YProp, clusters, deltadist=get.pepit("PRECISION"))
  }
  #clusters = remove_redundant_clusters (clusters)

  bsid = tools::file_path_sans_ext(bsfile)
  pepfile = paste(bsid, ".pdb", sep="")
  peptide = bio3d::read.pdb(pepfile)
  keep = select_clusters(clusters, target.data, bs.data, pdb, peptide, minalen)
  if (length(keep) > 0) {
      clusters = clusters[keep]
      allclusters = c(allclusters, clusters)
      scores = score_clusters(clusters, X, XProp, Y, YProp, deltadist=deltadist, bsfile=bsfile)
      allscores = rbind( allscores, scores)
   }
}
message("End of computation of all clusters")

if (nrow(allscores) == 0) {
   stop("no match found")
}
#colnames(allscores) = c("len", "alen", "coverage", "rmsd", "distorsion", "score", "lddt", "bsfile")

# sort allclusters; keep first NBHITS and equivalent last items
  message("Sort all clusters", length(allclusters))
  nbclusters = min(length(allclusters), get.pepit("NBHITS"))
  o = order(allscores[[criteria]], decreasing=TRUE)
  o = o[1:nbclusters]
  allscores =  allscores[o,]
  allclusters = allclusters[o]

  seqfile = paste(prefix, "sequences.faa",sep="")

# outputs
cat("index bs target precision bslen alen rmsd coverage meandist score lddt nres peplen\n", file=scorefile)

  for (i in 1:nrow(allscores)) {
       bsfile = allscores$bsfile[i]

bs.data = read.table(bsfile, header=TRUE, tryLogical=FALSE)
       # filter with HSE
       bs.data = bs.data[bs.data$hseu <= hse | bs.data$hsed <= hse,]
       chains = unique(bs.data$chain)
       bs.data$insert[is.na(bs.data$insert)]=""
       Y = bs.data[,c("x","y","z")]
       Y = as.matrix(Y)
       YProp = bs.data
       bsid = tools::file_path_sans_ext(bsfile)
       pepfile = paste(bsid, ".pdb", sep="")
       peptide = bio3d::read.pdb(pepfile)
       peplen = sum(peptide$calpha)

# output peptide pdbs
        I = (allclusters[[i]]-1)%%N+1 # target indices
        J = (allclusters[[i]]-1)%/%N+1# bs indices
        if (get.pepit("POSE")) {
            # output binding site posed on target in a .pdb file 
             ligand.moved = superpose_sites(allclusters[[i]], bs.data, target.data, peptide)
    	     #pdb1 = bio3d::trim.pdb(peptide, resno=target.data$resno[I], chain=target.data$chain[I], insert=target.data$insert[I])
   	     targetbs.pdb = bio3d::trim.pdb(pdb, resno=target.data$resno[I], chain=target.data$chain[I], insert=target.data$insert[I])
	     pepresno = NULL
	     if (get.pepit("PEPTRIM")) {
		   message("reducing peptide size")
		   pepsel = bio3d::binding.site(ligand.moved, targetbs.pdb, cutoff=get.pepit("CONTACT"), byres=TRUE)
		   if (0) {
		      peprm = bio3d::binding.site(ligand.moved, pdb, cutoff=1.4, byres=TRUE)
		      pepsel$resno = setdiff(pepsel$resno, peprm$resno)
		   }
		   if (!is.null(pepsel)) {
		      pepresno = seq(min(pepsel$resno), max(pepsel$resno))
		      peplen = length(pepresno)
		      ligand.moved = bio3d::trim.pdb(ligand.moved, resno=pepresno)
		   }
	     }
	     ligand.moved$atom$chain="P" # required for futher analysis (peptide is always chain P)
             pepfile = paste(prefix, "peptide-", i, ".pdb",sep="")
             write.pdb(pdb=ligand.moved, file = pepfile)
 	     if (get.pepit("PEPSEQ")) {
	     	seq = paste(aa321(ligand.moved$atom$resid[ligand.moved$calpha]),collapse="")
		cat(">",i, bsfile, "\n", file=seqfile, append=TRUE)
		cat(seq,"\n", file=seqfile, append=TRUE)
	     }
	}
# output matching (= non seq. alignment) in .al file
	if (get.pepit("ALIGNFILE")) {
             	      bs_res = paste(bs.data$resno[J],":",bs.data$insert[J],":",bs.data$elety[J],":",bs.data$chain[J], sep="")
             	      target_res = paste(target.data$resno[I],":",target.data$insert[I],":",target.data$elety[I],":",target.data$chain[I], sep="")
             	      cat(">",i, bsfile, tfile, deltadist, allscores$len[i], allscores$alen[i], allscores$rmsd[i], allscores$coverage[i], allscores$distorsion[i], allscores$score[i], allscores$pval[i],"\n", file=alignfile, append=TRUE)
             	      cat(bs_res, "\n", file=alignfile, append=TRUE)
             	      cat(target_res, "\n", file=alignfile, append=TRUE)
	 }
	 message("output score file")
# output score file 
        cat(i, bsfile, tfile, deltadist, allscores$len[i], allscores$alen[i], allscores$rmsd[i], allscores$coverage[i], allscores$distorsion[i], allscores$score[i], allscores$lddt[i], allscores$nres[i], peplen, "\n", file=scorefile, append=TRUE)
    }

