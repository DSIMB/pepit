
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

#' create pdb file with binding site atoms marked with b field>0 + ligand
#'
#' @param pdb  pdb file as read by bio3d
#' @param target_chains  vector of selected chains as target
#' @param ligand_chains  vector of selected chains as ligand
#' @param cutoff distance cutoff for defining contact atoms
#' @param add type of carbon atoms to be added to the binding sites
#' @param add.ligand add ligand to the output pdb
#' @param access if >0 keep only surface atoms on returned pdb
#'
#' @import bio3d
#' @return pdb bio3d structure
#' @export

get_binding_sites = function(pdb, target_chains, ligand_chains, cutoff=get.pepit("CONTACT"), add=get.pepit("ADD"), add.ligand=TRUE, access=-1.0, pepchain=NULL) {
    chains=unique(pdb$atom$chain)
    target_chains=intersect(target_chains, chains)
    ligand_chains=intersect(ligand_chains, setdiff(chains, target_chains))
    if (length(target_chains)==0 | length(ligand_chains)==0) {
      return(NULL)
    }

    add_values=c("", "calpha", "cbeta", "cab")
    add=intersect(add, add_values)
    if (length(add)==0) {
      message("add what?")
      return(NULL)
    }

    ligand.pdb=bio3d::trim.pdb(pdb, chain=ligand_chains, string="protein")
    target.pdb=bio3d::trim(pdb, string="protein", chain=target_chains)
    target.pdb$atom[,"o"]=0 # test
    selbs=bio3d::binding.site(target.pdb, ligand.pdb, cutoff=cutoff, byres=FALSE)

    if (length(selbs$inds$atom)==0) {
      message("no binding sites")
      return(NULL)
    }
    # atoms in contact
    inds=selbs$inds
    target.pdb$atom[inds$atom,"o"]=1 #test
    # if access>0 select atoms in contact and accessible
    if (access>0) {
      ind_surf=which(target.pdb$atom$b>=access)
      ind_surf=intersect(ind_surf, inds$atom)
      inds=list(atom=ind_surf, xyz=atom2xyz(ind_surf))
    }

    # add calpha atoms of selected atoms
    if (length(add)>0) {
      ind=get_carbons(target.pdb, inds, add)$atom
      target.pdb$atom[ind,"b"]=99
      ind=c(ind,inds$atom)
      ind=unique(ind)
      #ind=sort(ind)
      inds=as.select(ind)
    }
    
    pdb=bio3d::trim.pdb(target.pdb, inds=inds)
# add the ligand to the patch file for script SiteAlign.R
    if (add.ligand) {
   #ligand.pdb$atom$b[]=99
      ligand.pdb$atom$b=99 # utile ?
      if (!is.null(pepchain)) ligand.pdb$atom$chain=pepchain
      pdb=suppressWarnings(bio3d::cat.pdb(pdb, ligand.pdb, rechain=FALSE, renumber=FALSE))
    }

  pdb
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

get_binding_sites_2 = function(target, ligand, cutoff=get.pepit("CONTACT"), add=get.pepit("ADD")) {
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

