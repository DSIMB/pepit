#' Annotate surface atoms
#'
#' @param pdb a pdb structure returned by bio3d::read.pdb
#' @param inds a vector of atom indices (of pdb$atom) to be processed
#' @param add "calpha", "cab" or "" add alpha carbons of residues with accessible atoms
#' @description it returns indices of alpha carbons of residues of given atoms, (or alpha and cbeta atoms)
#' @description protein chains are processed independantly
#' @return atom indices , list of class select
#' @export
#'
#' @examples
get_carbons=function(pdb, inds, add="calpha") {
  #ind=inds$atom
  ind=c()
  res=unique(pdb$atom[inds$atom,c("resno", "chain")])
  chains=unique(res$chain)
  for (i in 1:length(chains)) {
       resno=res$resno[res$chain==chains[i]]
       chind=chain_find_carbons(pdb, chains[i], resno, add)
       ind=c(ind,chind$atom)
  }
  ind=unique(ind)
  ind=sort(ind)
  as.select(ind)
}

chain_find_carbons=function(pdb, chain, resno, add) {
    if (add=="calpha") sel=bio3d::atom.select(pdb, elety="CA", chain=chain, resno=resno)
    if (add=="cab") sel=bio3d::atom.select(pdb, elety=c("CA","CB"), chain=chain, resno=resno)
    if (add=="cbeta") sel=bio3d::atom.select(pdb, string="cbeta", chain=chain, resno=resno) #!backbone+cbeta
    sel
}
