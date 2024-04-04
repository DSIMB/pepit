#' Annotate surface atoms
#'
#' @param pdb a pdb structure returned by bio3d::read.pdb
#' @param chains vector of chains to be processed
#' @param minacc accessibility threshold. b field set to 0 if SASA is less than ACC.
#' @param probeRadius probe radius in Angstroms for the FreeSASA application
#' @param add type of carbon atoms to be added
#' @param select boolean for keeping atoms in output pdb
#' @description It replaces the b field by the solvent accessible surface area
#' @description protein chains are processed independantly
#' @description if chains!=NULL, only b field of chains in chains are set to acc. 
#' @description b field of the other chains are set to 0
#' @description it implies that binding sites between chains in contact are accessible
#' @return modified pdb structure (b field with SASA values)
#' @import bio3d
#' @import vanddraabe
#' @export
#'
#' @examples
access_atoms=function(pdb, chains=NULL, minacc=get.pepit("ACC"), probeRadius = get.pepit("radius"), add=get.pepit("ADD"), select=FALSE) {
  if (is.null(chains)) chains=unique(pdb$atom$chain)
  add_values=c("", "calpha", "cbeta", "cab")
  if (!add%in%add_values) {
    message("add what?")
    return(NULL)
  }
  
  pdb$atom$b[]=0 # otherwise all chain not in chains remains acc.
  pdb=chain_surf_annotate(pdb, chains[1], minacc, probeRadius)

  if (length(chains)>1) {
    for (i in 2:length(chains)) {
      pdb=chain_surf_annotate(pdb, chains[i], minacc, probeRadius)
      # surf.pdb=bio3d::cat.pdb(surf.pdb, tmp.pdb, renumber=FALSE, rechain=FALSE)
      # renumber=FALSE doesn't work when chain break
    }
  }
  
  ind=which(pdb$atom$b>=minacc)
  inds=as.select(ind)
  # add calpha atoms of selected atoms
  if (length(add)>0) {
    ind=get_carbons(pdb, inds, add)$atom
    pdb$atom[ind,"b"]=99
    ind=c(ind,inds$atom)
    ind=unique(ind)
    #ind=sort(ind)
    inds=as.select(ind)
  }

  if (select) {
    pdb=bio3d::trim.pdb(pdb, inds=inds)
  }
  #pdb$atom$b[pdb$calpha]=99
  return(pdb)
}

SASA=function(pdb, probeRadius = get.pepit("radius")) {
    result=vanddraabe::FreeSASA.diff(atoms.oi = pdb$atom, probeRadius = probeRadius)
    acc=result$SASA.prot
    x=strsplit(result$uniq.atom.ids,split="_")
    for (i in 1:length(acc)) {
        resno=as.integer(x[[i]][2])
        chain=x[[i]][3]
        eleno=as.integer(x[[i]][5])
	      ind=bio3d::atom.select(pdb, resno=resno, chain=chain, eleno=eleno)
	      #ind=atom.select(pdb,  eleno=eleno)
	      pdb$atom[ind$atom,"b"]=acc[i]
    }
    pdb
}

chain_surf_annotate=function(fullpdb, chain, minacc, probeRadius) {
    selchain=bio3d::atom.select(fullpdb, chain=chain)
    pdb=bio3d::trim.pdb(fullpdb, chain=chain)
    modified.pdb=SASA(pdb, probeRadius)# b-factors replaced by acc
    modified.pdb$atom$b[modified.pdb$atom$b<minacc]=0
    fullpdb$atom[selchain$atom,]=modified.pdb$atom
    fullpdb
}

#' select surface atoms
#'
#' @param pdb a pdb structure returned by bio3d::read.pdb
#' @param chains vector of chains to be processed
#' @param minacc double accessibility threshold
#' @description it selects accessible atoms (access>=minacc) and return a pdb structure
#' @return pdb structure with non accessible atoms removed
#' @import bio3d
#' @export
#'
#' @examples
trim_non_acc_atoms=function(pdb, chains=NULL, minacc=get.pepit("ACC")) {
  if (is.null(chains)) chains=unique(pdb$atom$chain)
  if (minacc>=0) {
    ind=which(pdb$atom$b>=minacc & pdb$atom$chain%in%chains)
    inds=as.select(ind)
    pdb=bio3d::trim.pdb(pdb, inds=inds)
  }
  pdb
}

