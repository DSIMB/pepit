
#' create pdb file with binding site atoms marked with b field>0 + ligand
#'
#' @param pdb  pdb file as read by bio3d
#' @param target_chains  vector of selected chains as target
#' @param ligand_chains  vector of selected chains as ligand
#' @param cutoff distance cutoff for defining contact atoms
#' @param add type of carbon atoms to be added to the binding sites
#' @param add.ligand add ligand to the output pdb
#' @param surface if >0 keep only surface atoms on returned pdb
#'
#' @import bio3d
#' @return pdb bio3d structure
#' @export

get_binding_sites = function(pdb, target_chains, ligand_chains, cutoff=5, add="calpha", add.ligand=TRUE, surface=-1.0) {
    chains=unique(pdb$atom$chain)
    target_chains=intersect(target_chains, chains)
    ligand_chains=intersect(ligand_chains, setdiff(chains, target_chains))
    if (length(target_chains)==0 | length(ligand_chains)==0) {
      return(NULL)
    }

    add_values=c("", "calpha", "cbeta", "cab","backbone", "residue", "deep")
    add=intersect(add, add_values)
    if (length(add)==0) {
      message("add what?")
      return(NULL)
    }

    ligand.pdb=bio3d::trim.pdb(pdb, chain=ligand_chains, string="protein")
    target.pdb=bio3d::trim(pdb, string="protein", chain=target_chains)
    selbs=bio3d::binding.site(target.pdb, ligand.pdb, cutoff=cutoff, byres=FALSE)

    if (length(selbs$inds$atom)==0) {
      message("no binding sites")
      return(NULL)
    }
    inds=selbs$inds
    if (length(add)>0) {
      ind=get_carbons(target.pdb, inds, add)$atom
      ind=c(ind,inds$atom)
      ind=unique(ind)
      ind=sort(ind)
      inds=as.select(ind)
    }

    if (surface>0) {
      ind_surf=which(target.pdb$atom$b>=surface)
      ind_surf=intersect(ind_surf, inds$atom)
      inds=list(atom=ind_surf, xyz=atom2xyz(ind_surf))
    }
    pdb=bio3d::trim.pdb(target.pdb, inds=inds)
    print(pdb)
# add the ligand to the patch file for script SiteAlign.R
    if (add.ligand) {
   #ligand.pdb$atom$b[]=99
      ligand.pdb$atom$b=99 # utile ?
      pdb=suppressWarnings(bio3d::cat.pdb(pdb, ligand.pdb, rechain=FALSE, renumber=FALSE))
    }

  pdb
}

