# output a complex with first chain A and second chain P if any if 4 arguments
# output a pdb with the selected chain if 2 arguments
library(bio3d)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  cat("usage: mergePDB pdbfile1 chain1 [pdbfile2 chain2]\n")
  q()
}

pdbfile1 = args[1]
chain1 = args[2]
ext1 = tools::file_ext(pdbfile1)

if (ext1 == "") pdb = read.pdb(pdbfile1)
if (ext1 == "pdb") pdb = read.pdb(pdbfile1)
if (ext1 == "cif") pdb = suppressWarnings(read.cif(pdbfile1))
pdb1 = trim.pdb(pdb, chain=chain1)
pdb1 = bio3d::trim.pdb(pdb1, string="protein")
pdb1 = bio3d::trim.pdb(pdb1, string="noh") 
pdb1 = clean.pdb(pdb1, fix.chain = TRUE, fix.aa = TRUE)
pdball = pdb1

pdb2 = NULL
if (length(args) == 4) {
   pdbfile2 = args[3]
   chain2 = args[4]
   ext2 = tools::file_ext(pdbfile2)
   if (pdbfile2 != pdbfile1) {
      if (ext2 == "pdb") pdb = read.pdb(pdbfile2)
      if (ext2 == "cif") pdb = suppressWarnings(read.cif(pdbfile2))
      pdb = bio3d::trim.pdb(pdb, string="protein")
      pdb = bio3d::trim.pdb(pdb, string="noh") 
   }
   pdb2 = trim.pdb(pdb, chain=chain2)
   clean.pdb(pdb2, fix.chain = TRUE, fix.aa = TRUE)
   pdb1$atom$chain = "A"
   pdb2$atom$chain = "P"
   pdball = cat.pdb(pdb1, pdb2, renumber=FALSE, rechain=FALSE)
}

pdball = clean.pdb(pdball, fix.chain = FALSE, fix.aa = TRUE)

write.pdb(pdball, file="output.pdb")
