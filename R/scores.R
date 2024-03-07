#' Compute BC score between X and Y
#'
#' @param X matrix Nx3
#' @param Y matrix Nx3
#' @export
BCscore<-function(X,Y) {
  X=scale(X,scale=FALSE)
  Y=scale(Y,scale=FALSE)
  if (det(t(X)%*%X)<1 | det(t(Y)%*%Y)<1) return(0)
  det(t(X)%*%Y)/(sqrt(det(t(X)%*%X))*sqrt(det(t(Y)%*%Y)))
}

#' Compute rmsd between X and Y
#'
#' @param X matrix Nx3
#' @param Y matrix Nx3
#' @export
rmsd<-function(X,Y) {
  bio3d::rmsd(as.vector(t(X)),as.vector(t(Y)), fit=TRUE)
}

#' evaluate cliques and clusters
#'
#' @param clusters list of mappings : a mapping is a set of correspondence identifiers
#' @param X atom coordinate matrix from target (Nx3)
#' @param Y atom coordinate matrix from query (a binding site) (Mx3)
#' @param score_function score optimized by maximal matching procedure
#' @import bio3d
#' @import igraph
#' @import Rcpp
#' @return list of mapping criteria :
#' @export
#'
#' @examples
score_clusters = function(clusters, X, Y, score_function=get.pepit("SCORE"), deltadist=1.0) {
    N = nrow(X)
    M = nrow(Y)
    results = list()
    for (i in 1:length(clusters)) {
        C = clusters[[i]]
        nC = length(C)
        score = matching_score(C, X, Y, score_function, deltadist)
        dist = deltadist*(1-1/nC**2*score)
        #score = sqrt(score)
        coverage = nC/M
        I = (clusters[[i]]-1)%%N+1 # target indices
        J = (clusters[[i]]-1)%/%N+1# query indices
        rms = rmsd(X[I,],Y[J,])
        results$len = append(results$len,M)
        results$alen = append(results$alen,nC)
        results$coverage = append(results$coverage, coverage)
        results$rmsd = append(results$rmsd, rms)
        results$distorsion = append(results$distorsion, dist)
        results$score = append(results$score, score)
    }
    results
}

#' evaluate pvalues from all scores
#'
#' @param score vector of scores
#' @return vector of p values
#' @export
#'
#' @examples
p_values=function(score, loc, scale) {
  m = mean(score)
  s = sd(score)
  scale = sqrt(6)/pi*s
  loc = m + scale*digamma(1) # digamma(1)=-Euler
  sc = (score - loc)/scale
  1 - exp(-exp(-sc))
}

#' evaluate evd
#'
#' @param score a
#' @return parameters of extreme distribution loc and scale
#' @export
#'
#' @examples
evd=function(score) {
  m = mean(score)
  s = sd(score)
  scale = sqrt(6)/pi*s
  loc = m + scale*digamma(1) # digamma(1)=-Euler
  list(scale=scale, loc=loc)
}



