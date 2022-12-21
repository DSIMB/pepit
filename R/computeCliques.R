
extract_cliques<-function(gp) {
  #Cl=igraph::largest_cliques(gp)
  size=igraph::clique.number(gp)
  cat("extract_cliques-->", igraph::clique.number(gp),"\n")
  Cl=igraph::maximal.cliques(gp,min=size-1)
  #Cl=cliques(gp, min=6,max=6)
  if (length(Cl)==0) return(list(clique=NULL,size=0))
  ##
  for (i in 1:length(Cl)) Cl[[i]]=as.integer(names(Cl[[i]]))
  ##
  list(clique=Cl, size=max(mapply(length,Cl)))
}

union_clique <- function(C1,C2,N) {
  I1=(C1-1)%%N+1
  I2=(C2-1)%%N+1
  J1=(C1-1)%/%N+1
  J2=(C2-1)%/%N+1
  I=which(!I2%in%I1 & !J2%in%J1)
  base::union(C1,C2[I])
}
#
# sans BC
#
merge_cliques<-function(clusters,N,common=patchsearch_MERGECLIQUE) {
  cat("------> merge_cliques", length(clusters),"\n")
  stop=FALSE
  while (!stop) {
    stop=TRUE
    n=length(clusters)
    i=1
    while (i<= n-1) {
      j=i+1
      while (j<= n) {
        if (length(base::intersect(clusters[[i]],clusters[[j]]))>=common) {
          #clusters[[i]]=union(clusters[[i]],clusters[[j]]) # may create forks
          C=union_clique(clusters[[i]],clusters[[j]],N)
          clusters[[i]]=C
          clusters[[j]]=c()
          n=n-1
          stop=FALSE
        }
        j=j+1
      }
      i=i+1
    }
  }
  cat("------> merge_cliques", length(clusters),"\n")
  return(clusters)
}

clique_patchsearch=function(g, X, Y, minclique=patchsearch_MINCLIQUE, bcmin=patchsearch_BCMIN) {
  cat("clique_patchsearch::MINCLIQUE", patchsearch_MINCLIQUE, "\n")
  N=nrow(X);M=nrow(Y)
  cl=extract_cliques(g)
  ##
  if (cl$size<minclique) {
    cat("cl$size<minclique\n")
    return(NULL)
  }
  nbclique=length(cl$clique)
  size=bc=double(nbclique)
  for (i in 1:nbclique) {
    C=cl$clique[[i]]
    size[i]=length(C)
    I=(C-1)%%N+1 # target indices
    J=(C-1)%/%N+1# query indices
    bc[i]=BCscore(X[I,], Y[J,])
  }
  cl$clique=cl$clique[bc>=bcmin]
  #
  if (length(cl$clique)==0) return(NULL)
  bc=bc[bc>=bcmin]
  o=order(bc,decreasing=TRUE)
  cl$bc=bc[o]
  cl$clique=cl$clique[o]
  cl$size=size[o]
  return(cl)
}

#
# sort clusters
#
sort_clusters=function(clusters,N, X, Y, score_function=patchsearch_SCORE, thresh=0) {
    nbclique=length(clusters)
    cat("nb of clique clusters=",nbclique,"\n")
    bc=alen0=double(nbclique)
    for (i in 1:nbclique) {
    	C=clusters[[i]]
   	  size=length(C)
   	  I=(C-1)%%N+1 # target indices
   	  J=(C-1)%/%N+1# query indices
   	  s=sign(BCscore(X[I,], Y[J,]))
	    bc[i]=s*matching_score(C, X, Y, score_function, thresh) 
	    alen0[i]=length(C)
    }
    clusters=clusters[bc>=0]
    bc=bc[bc>=0]
    o=order(bc,decreasing=TRUE)
    clusters=clusters[o]
    clusters
}

#
# deltadist inutile ?
#
build_graph=function(types, X, XProp, Y, YProp, mindist, maxdist, deltadist, verbose=FALSE) {
  N=nrow(X)
  J=which(YProp%in%types)
  I=which(XProp%in%types)
  cat ("build_graph::",length(I),length(J),"\n")
  Mmotif=length(J)
  Nmotif=length(I)
  if (Mmotif<patchsearch_MINCLIQUEATOMS | Nmotif<patchsearch_MINCLIQUEATOMS) {
    types=c("a","b","A","C","O","N")
    J=which(YProp%in%types)
    I=which(XProp%in%types)
   cat ("+build_graph::",length(I),length(J),"\n")
  }
  V=vertex(XProp[I], YProp[J])
  if (verbose) cat("vertices:",length(I), length(J),length(V)/2,"\n")
  if (all(V==0) | length(V)<=1) return(NULL)
  V=cbind(J[V[,1]],I[V[,2]]) # atom ids 1..M, 1..N
  #nV=nrow(V)
  E=buildGraph(X, Y, V, as.double(deltadist), as.double(mindist), as.double(maxdist))
  while (nrow(E)>patchsearch_MAXEDGES){
    cat("edges:",nrow(E),"\n")
    maxdist=maxdist/2.;
    E=buildGraph(X, Y, V, deltadist, mindist, maxdist)
  }
  cat("edges:",nrow(E),"\n")
  if (nrow(E)==0) {
    message("no edge...")
    return(igraph::make_empty_graph())
  }
  P=cbind((E[,3]-1)*N+E[,1],(E[,4]-1)*N+E[,2])
  mode(P)="character" #vertex label and not id to avoid huge graphs with non connected nodes
  gp=igraph::graph.edgelist(P, directed=FALSE)
  if (verbose) cat("graph:vertices:",igraph::vcount(gp), "graph.edges:", igraph::ecount(gp),"\n")
  return(gp)
}
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

matching_score=function(C, X, Y, score_function, thresh=0) {
    N=nrow(X)
    CI=(C-1)%%N+1 # target indices
    CJ=(C-1)%/%N+1# query indices
    W=score_function(X, CI, CI, Y,  CJ, CJ,thresh)
    sum(W)
}

#' cliques
#'
#' @param X atom coordinate matrix from target (Nx3)
#' @param XProp atom label vector from target (N)
#' @param Y atom coordinate matrix from query (a binding site) (Mx3)
#' @param YProp atom label vector from query  (M)
#' @param precision parameter ?
#' @param deltadist precision parameter
#' @param mindist minimum atom distance considered
#' @param maxdist maximal atom distance considered
#' @param types atom types from which to build graph
#' @import bio3d
#' @import igraph
#' @import Rcpp
#' @return list of clusters : a cluster is a list of X indices (cluster$I) and of Y indices (cluster$J)
#' @export
#'
#' @examples
cliques=function(X, XProp, Y, YProp, deltadist, mindist, maxdist, types, verbose=FALSE) {
    clusters=list()
    N=nrow(X)
  #
  # construction du graphe produit
  #
  graph=build_graph(types, X, XProp, Y, YProp, mindist, maxdist, deltadist, verbose)
  if (length(igraph::E(graph))==0) {# from make_empty_graph()
    if (verbose) cat("empty graph...\n")
    clusters[[1]]=-1
    return(clusters)
  }
  if (verbose) cat("clique search...\n")
  #
  # recherche des cliques
  #
  cliques=clique_patchsearch(graph, X, Y, patchsearch_MINCLIQUE, patchsearch_BCMIN)
  if (verbose) cat("number of cliques=",length(cliques$clique),"\n")
  if (is.null(cliques)) {
    message("no clique found")
    clusters[[1]]=-2
    return(clusters)
  }
  cliques=cliques$clique
  #
  # merging cliques
  #
  if (verbose) cat("clique merge...\n")
  clusters=merge_cliques(cliques,N, patchsearch_INTERCLIQUE)
  nbclique=length(clusters)
  cat("nbclique", nbclique, "\n")
  clusters=sort_clusters(clusters,N,X,Y, patchsearch_SCORE)
  nbclique=length(clusters)
  cat("nbclique", nbclique, "\n")
  if (nbclique==0) {
    message("no valid clique after merging")
    clusters[[1]]=-3
    return(clusters)
  }
  return(clusters)
}

#' convert list of atom pairs (edges of bipartite mapping) to a list of atom pairs
#'
#' @param clusters : vector list of atom pairs
#' @return list of mappings : a mapping is a list of X indices (mapping$I) and of corresponding Y indices (mapping$J)
#' @export
#' @examples
get_mapping=function(clusters, N) {
  mapping=list()
  for (i in 1:length(clusters)) {
    I=(clusters[[i]]-1)%%N+1 # target indices
    J=(clusters[[i]]-1)%/%N+1# query indices
    mapping[[i]]=list(I=I,J=J)
  }
  return(mapping)
}
