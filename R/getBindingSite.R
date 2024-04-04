
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


#' Annotate surface atoms
#'
#' @param pdb a pdb structure returned by bio3d::read.pdb
#' @param inds a vector of atom indices (of pdb$atom) to be processed
#' @param add "calpha", "cab" or "" add alpha carbons of residues with accessible atoms or calpha and cbeta or backbones plus cbeta atom
#' @description it returns indices of supplementary atoms  of residues of given atoms
#' @description protein chains are processed independantly
#' @return atom indices relative to pdb$atom
#' @export
#'
#' @examples
get_atoms=function(pdb, resno, chain, add="calpha") {
  ind = c()
  for (i in 1:length(resno)) {
      if (add=="calpha") sel=bio3d::atom.select(pdb, elety="CA", chain=chain[i], resno=resno[i])
      if (add=="cab") sel=bio3d::atom.select(pdb, elety=c("CA","CB"), chain=chain[i], resno=resno[i])
      if (add=="cbeta") sel=bio3d::atom.select(pdb, string="cbeta", chain=chain[i], resno=resno[i]) #!backbone+cbeta
      if (add=="residue") sel=bio3d::atom.select(pdb, chain=chain[i], resno=resno[i]) #!backbone+cbeta
      ind = c(ind,sel$atom)
  }
  ind = unique(ind)
  ind = sort(ind)
}

#' create pdb file with binding site atoms marked with b field>0 + ligand
#'
#' @param target pdb structure as read by bio3d
#' @param ligand pdb structure as read by bio3d
#' @param cutoff distance cutoff for defining contact atoms
#' @param add type of atoms to be added to the binding sites
#'
#' @import bio3d
#' @return inds selected indices from target 
#' @export

get_binding_sites = function(target, ligand, cutoff=get.pepit("CONTACT"), add=get.pepit("ADD")) {
    add_values=c("", "calpha", "cbeta", "cab")
    add=intersect(add, add_values)
    if (length(add)==0) {
      message("add what?")
      return(NULL)
    }

    selbs=bio3d::binding.site(target, ligand, cutoff=cutoff, byres=FALSE)

    if (length(selbs$inds$atom)==0) {
      message("no binding sites")
      return(NULL)
    }
    # atoms in contact
    inds=selbs$inds

    # add calpha atoms of selected atoms
    if (length(add)>0) {
      ind=get_atoms(target, selbs$resno, selbs$chain, add)
      ind=c(ind,inds$atom)
      ind=unique(ind)
      inds=as.select(ind)
    }
    
    inds
}

