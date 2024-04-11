#' Atom typing
#'
#' @param pdb a pdb structure returned by bio3d::read.pdb
#' @description The function produces a data frame with pdb ATOM information
#' @description and add columns with HSE values, physico-chemcal encoding of residues and an atom type column:
#' @description alpha carbons changed to A, beta carbons changed to B,
#' @description other backbone carbons remains C,
#' @description aromatic carbons changed to a, other side-chain atoms changed to
#' @description to lower case (o,n,c)
#' @description ligand chain elements are not modified
#' @description 
#' @return a data frame
#' @import bio3d
#' @export
#'
#' @examples
encode = function(pdb, cutoff=get.pepit("HSECUTOFF"),  pc.code=get.pepit("PCCODE") , chain=NULL) {
  
    if (!is.null(chain)) pdb = bio3d::trim(pdb, chain=chain)
    # HSE
    hse = HSEB(pdb, cutoff=cutoff)

    # atom type
    pdb = type_atoms(pdb, ligchain=NULL)

    # physico-chemical code
    resid = pdb$atom$resid
    resid = aa321(resid)
    pc = pc.code[resid]

    eleno = pdb$atom$eleno
    elety = pdb$atom$elesy # modified by type_atoms
    resid = pdb$atom$resid
    chain = pdb$atom$chain
    resno = pdb$atom$resno
    insert = pdb$atom$insert
    x = pdb$atom$x
    y = pdb$atom$y
    z = pdb$atom$z

    hseu = integer(nrow(pdb$atom))
    hsed = integer(nrow(pdb$atom))
    r = resno[pdb$calpha]
    for (i in 1:length(r)) {
	    sel = atom.select(pdb, resno=r[i])
	    hseu[sel$atom] = hse$hsebu[i]
	    hsed[sel$atom] = hse$hsebd[i]
    }

    data.frame(eleno, elety, resid, chain, resno, insert, x, y, z, pc, hseu, hsed)
}
