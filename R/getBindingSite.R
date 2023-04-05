
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

