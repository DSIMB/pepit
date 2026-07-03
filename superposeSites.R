#' superpose query (binding site + ligand) onto the target. 
#' The ligand is kept while the binding site is removed.
#' 
#' @param cluster indices of mapping edges
#' @param bs data frame containing bs features (mobile)
#' @param target data frame containing target features (fixed)
#' @param  peptide pdb file containing peptide (or ligand)
#' @return a pdb structure containing the moved ligand 
#' @export
#'
#' @examples

# target : eleno, elety, resid, chain, resno, insert, x, y, z, pc, hseu, hsed

superpose_sites = function(cluster, bs, target, peptide) {
  N = nrow(target)
  I = (cluster-1)%%N+1 # target indices
  J = (cluster-1)%/%N+1# query indices
  target.ind = as.select(I)$xyz
  bs.ind = as.select(J)$xyz
  target.xyz = as.vector(rbind(target$x, target$y, target$z))
  bs.xyz = as.vector(rbind(bs$x, bs$y, bs$z))
  mobile.xyz = c(bs.xyz, peptide$xyz)
  moved.xyz = fit.xyz(fixed=target.xyz, mobile=mobile.xyz, fixed.inds=target.ind, mobile.inds=bs.ind, verbose=FALSE)
  
  hitrmsd = bio3d::rmsd(moved.xyz, b=target.xyz, a.inds=bs.ind, b.inds=target.ind)
  #cat("-----------> bs rmsd=", hitrmsd,"\n")

  pep.ind = (length(bs.xyz)+1):(length(mobile.xyz))
  peprmsd=bio3d::rmsd(peptide$xyz, b=moved.xyz, b.inds=pep.ind)
  #cat("-----------> pep rmsd=", peprmsd,"\n")
  
  peptide$xyz = moved.xyz[(length(bs.xyz)+1):(length(mobile.xyz))]
  peptide
}

#' adjust peptide to aligned binding site
trim_peptide = function(cluster, bs, peptide) {
}

#'compute the number of ligand atoms too close to the target atoms
#' 
#' @param ligand pdb structure
#' @param target pdb structure
#' @return number of clashes
#' @export
#'
#' @examples
clashes = function(ligand, target, cutoff=1.4) {
  result = bio3d::binding.site(ligand, target, cutoff=cutoff, byres=FALSE)
  nclash = length(result$inds$atom)  
  nclash
}