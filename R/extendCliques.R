matching_score=function(C, X, Y,  score_function, thresh=0) {
  N=nrow(X)
  CI=(C-1)%%N+1 # target indices
  CJ=(C-1)%/%N+1# query indices
  W=score_function(X, CI, CI, Y, CJ, CJ, thresh)
  sum(W)
}

filter_greedy=function(C, X, Y, score_function, thresh) {
  N=nrow(X)
  CI=(C-1)%%N+1 # target indices
  CJ=(C-1)%/%N+1# query indices
  W=score_function(X, CI, CI, Y, CJ, CJ, thresh)
  #maxs=sum(W)
  ind=which.min(W)
  while (W[ind]<=0) {
    CI=CI[-ind];CJ=CJ[-ind];C=C[-ind]
    W=score_function(X, CI, CI, Y, CJ, CJ,thresh)
    ind=which.min(W)
  }
  C
}

max_bipartite_score=function(C, K, X, Y, gp, score_function, thresh) {
  C0=C
  N=nrow(X)
  CI=(C-1)%%N+1 # target indices
  CJ=(C-1)%/%N+1# query indices
  KI=(K-1)%%N+1
  KJ=(K-1)%/%N+1
  W=score_function(X, KI, CI, Y, KJ, CJ, thresh)
  W[W<=0]=-1
  #
  igraph::E(gp)$weight=W
  igraph::V(gp)$type=igraph::bipartite_mapping(gp)$type
  result=igraph::max_bipartite_match(gp)
  x=result$matching[igraph::V(gp)$type]
  x=x[!is.na(x)]
  jnd=as.integer(substring(x,first=2))
  ind=as.integer(substring(names(x),first=2))
  #W[W<=0]=0 # inutile
  list(C=(jnd-1)*N+ind, weight=W, score=result$matching_weight)
}

max_bipartite=function(C, K, X, Y, thresh, verbose=FALSE, score_function=get.pepit("SCORE")) {
  Nref=min(nrow(X), nrow(Y))
  N=nrow(X)    

  KI=(K-1)%%N+1
  KJ=(K-1)%/%N+1
  P=cbind(paste("q",KJ,sep=""),paste("t",KI,sep=""))
  gp=igraph::graph.edgelist(P, directed=FALSE)
  stop=FALSE
  niter=0
  Cprev=C
  C0=C
  stop=FALSE
  niter=0
  score_init=matching_score(C,X, Y, score_function, thresh)
  if (verbose) cat("C0=",length(C),score_init, score_init/length(C),"\n")
  while(!stop & niter<20) {
    result=max_bipartite_score(Cprev, K, X, Y,  gp, score_function, thresh) #f(C1=result$C,C0=Cprev)
    nC=length(result$C)
    if (result$score<=0) stop=TRUE
    if(result$score>0) C=filter_greedy(result$C,X,Y,score_function, thresh) # C1bar
    if (verbose) cat("C0=",length(Cprev),"C1=",length(result$C),"C1b=",length(C),"C0.C1=",length(intersect(result$C,Cprev)),"C0.C1b=",length(intersect(C,Cprev)),"f(C1,C0)=",result$score,"\n")
    if (all(Cprev%in%C) & all(C%in%Cprev)) {stop=TRUE}
    if (verbose) cat("f(C1)=", matching_score(result$C,X, Y,score_function,thresh),"f(C1b)=", matching_score(C,X,Y,score_function,thresh),"\n")
    Cprev=C
    niter=niter+1
  }
  nC=length(C)
  score=matching_score(C,X, Y, score_function,thresh)
  dist=thresh*(1-1/nC**2*score)
  if (stop & verbose) cat ("convergence\n")
  C=union(C0,C) # hope that clique comes first
  list(C=C, score=score/nC, dist=dist, coverage=nC/Nref)
}


max_bipartite_simple_greedy=function(C, K, X, Y, thresh, verbose=FALSE, score_function=get.pepit("SCORE")) {
  Nref=min(nrow(X), nrow(Y))
  N=nrow(X)
  CI=(C-1)%%N+1 # target indices
  CJ=(C-1)%/%N+1# query indices
  nC=length(C)
  K=setdiff(K,C)
  KI=(K-1)%%N+1
  KJ=(K-1)%/%N+1
  stop=FALSE
  score_init=matching_score(C,X, Y, score_function,thresh)
  #cat("C=",length(C),score_init, score_init/length(C),"\n")
  W=score_function(X, KI, CI, Y, KJ, CJ, thresh)
  W1=score_function(X, CI, CI, Y, CJ, CJ, thresh)

  while(!stop) {
    stop=TRUE
    ind=which.max(W)
    w=W[ind]
    nC=length(C)
    score=matching_score(C, X, Y,  score_function, thresh)
    #cat("score =", w, score,  score/nC, "\n")
    if (w>0) {
      #if (w>(score/nC-1)/2) #!si wii=1
      CI=c(CI,KI[ind])
      CJ=c(CJ,KJ[ind])
      C=c(C,K[ind])
      nC=length(C)
      KI=KI[-ind]
      KJ=KJ[-ind]
      K=K[-ind]
      if (length(K)>0) {
        W=score_function(X,  KI, CI, Y, KJ, CJ, thresh)
        stop=FALSE
      }
    }
  }
  score=matching_score(C, X, Y, score_function, thresh)
  dist=thresh*(1-1/nC**2*score)
  if (stop & verbose) cat ("convergence\n")
  list(C=C, score=score/nC, dist=dist, coverage=nC/Nref)
}


#' extends cliques
#'
#' @param X atom coordinate matrix from target (Nx3)
#' @param XProp atom label vector from target (N)
#' @param Y atom coordinate matrix from query (a binding site) (Mx3)
#' @param YProp atom label vector from query  (M)
#' @param cliques list of mappings that are lists of corresponding edge indices 
#' @param deltadist precision parameter
#' @import bio3d
#' @import igraph
#' @import Rcpp
#' @return list of clusters : a cluster is a list of X indices (cluster$I) and of Y indices (cluster$J)
#' @export
#'
#' @examples
extend_cliques=function(X, XProp, Y, YProp, cliques, deltadist, score_function=get.pepit("SCORE"), verbose=FALSE) {
  N=nrow(X)
  nbclique=length(cliques)
  clusters=list()

  V=vertex(XProp, YProp, mode=get.pepit("MODE"), size=0, hse=100) # hse=100 => +buried atoms 
  for (ic in 1:nbclique) {
    cat("enlarging clique ", ic, "\n")
    C=cliques[[ic]]
    nbefore= length(C)
    #
    # calcul de K=ensemble de links à tester pour enrichir clique choisie
    # 1. ensemble des voisins de C
    #
    #K=selectLinks(C, V, X, Y, deltadist, get.pepit("NBNEI"))
    #K=union(K,C)
    K=(V[,1]-1)*N+V[,2]
    result=max_bipartite(C, K, X, Y, thresh=deltadist, verbose, score_function)
    #result=max_bipartite_simple_greedy(C, K, X, Y, deltadist, verbose=FALSE, score_function)
    clusters[[ic]]=result$C
    #clen[ic]=length(C)
  }
  return(clusters)
}

#
# remove redundant clusters
#
#' extends cliques
#'
#' @param clusters list of cluster
#' @param common number of pairs of atoms in common to merge cluster
#' @return list of clusters 
#' @export
#'
#' @examples
remove_redundant_clusters = function(clusters, common=get.pepit("INTERCLIQUE")) {
  o=order(unlist(lapply(clusters, length)), decreasing=TRUE)
  clusters=clusters[o]
  n=length(clusters)
  i=1
  while (i<= n-1) {
    j=i+1
    while (j<= n) {
      lmin=min(length(clusters[[i]]),length(clusters[[j]]),common)
      #lmin=common
      if (length(base::intersect(clusters[[i]],clusters[[j]]))>=lmin) {
        clusters[[j]]=c()
        n=n-1
        j=j-1
      }
      j=j+1
    }
    i=i+1
  }
  return(clusters)
}
