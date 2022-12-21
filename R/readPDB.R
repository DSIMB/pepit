#' read a PDB file into a list
#'
#' @param pdbfile  pdb file
#' @param chain  vector of selected chains
#' @param access boolean to select atoms with b-factor > 0
#' @param resi vector of residu numbers of selected residus
#'
#' @import bio3d
#' @return list of atoms (PepiT format)
#' @export
read_pdb <- function(pdbfile, chain="",  access=FALSE, resi=integer(0)) {
	pdb=suppressWarnings(bio3d::read.pdb(pdbfile))
	cat("read", nrow(pdb$atom),"atoms\n")
	if (length(resi)>0) pdb=bio3d::trim.pdb(pdb, resno=resi)
	if (access) {
	   cat("trim to accessible atoms", access, nrow(pdb$atom),"-->")
	   pdb=bio3d::trim.pdb(pdb, eleno=pdb$atom$eleno[pdb$atom$b>0])
	   cat(nrow(pdb$atom),"\n")
	}
	pdb=bio3d::trim.pdb(pdb, string="protein")
	if (chain[1]!="") {
	   pdb=bio3d::trim.pdb(pdb, string="protein", chain=chain)
	   cat("trim to ",chain, nrow(pdb$atom),"atoms\n")
	}
	x=pdb$atom$x
	y=pdb$atom$y
	z=pdb$atom$z
	access=pdb$atom$b
	type=gsub(" ","",pdb$atom$elesy)
	n=length(type)
	resname=aa321(pdb$atom$resid)
	ch=pdb$atom$chain
	ch[is.na(ch)]=" "
	pdb$atom$chain=ch
	# modif du 4/7/18 (sinon non bijection atomnum et resnum:chain:atomname)
	pdb$atom$insert[is.na(pdb$atom$insert)]=""
	resnum=paste(pdb$atom$resno,":",pdb$atom$insert,sep="")
	#
  list(coord=cbind(x,y,z), type=type, access=access, num=pdb$atom$eleno, resno=pdb$atom$resno, resnum=resnum, insert=pdb$atom$insert, resname=resname, natoms=n, aname=pdb$atom$elety, chains=ch)
}
