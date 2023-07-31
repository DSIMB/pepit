#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double distloc(NumericMatrix X, int ind1, int ind2) {
  int i;
  double x,y,sum=0;
  for (i=0;i<3;i++) {
    x=X(ind1-1,i);
    y=X(ind2-1,i);
    sum+=(x-y)*(x-y);
  }
  return sqrt(sum);
}

// [[Rcpp::export]]
int gaploc(IntegerVector ResNo, int ind1, int ind2) {
  int x,y;
  x=ResNo(ind1-1);
  y=ResNo(ind2-1);
  return abs(x-y);
}

//' Constructs vertices of correspondence graph
//'
//' @param XProp size N vector of atom types
//' @param YProp size M vector of atom types
//' @param V output matrix Nvx2 of vertices
//' @export
// [[Rcpp::export]]
IntegerMatrix vertex(CharacterVector XProp,  CharacterVector YProp) {
  int i,ip,v=0;
  int N=XProp.size();
  int M=YProp.size();
  IntegerMatrix V(N*M,2);
  for(i=0; i<M; i++) {
    for (ip=0; ip<N; ip++) {
      if (strcmp(XProp(ip),YProp(i))==0) {
        V(v,0)=i+1;
        V(v,1)=ip+1;
        v++;
      }
    }
  }
  if (v<=1) {
        fprintf(stderr,"no product graph...\n");
        return V(Range(0,0),_);
  }
  fprintf(stdout, "---> vertices of correspondence graph = mapping edges: %d\n", V.nrow());
  return V(Range(0,v-1),_);
}

//' Constructs vertices of correspondence graph
//'
//' @param XProp size N vector of atom types 
//' @param XResname size N vector of residue names 
//' @param YProp size M vector of atom types 
//' @param YResname size M vector of residue names 
//' @param V output, matrix Nvx2 of vertices
//' @export
// [[Rcpp::export]]
IntegerMatrix vertex_ho(StringVector XProp, StringVector XResname, StringVector YProp, StringVector YResname) {
  int i,ip,v=0;
  int N=XProp.size();
  int M=YProp.size();
  IntegerMatrix V(N*M,2);
  for(i=0; i<M; i++) {
    for (ip=0; ip<N; ip++) {
      if (strcmp(XProp(ip),YProp(i))==0 && strcmp(XResname(ip),YResname(i))==0) {
        //Rcpp::Rcout << XProp(ip) << "," << YProp(i) << ","<< XRes(ip) << ","<< YRes(i) << "\n";
        V(v,0)=i+1;
        V(v,1)=ip+1;
        v++;
      }
    }
  }
  if (v<=1) {
    fprintf(stderr,"no product graph...\n");
    return V(Range(0,0),_);
  }
  return V(Range(0,v-1),_);
}

//' Constructs vertices of correspondence graph based on residues
//'
//' @param XResno size N vector of residue number
//' @param XResname size N vector of residue names 
//' @param YResno size M vector of residue number
//' @param YResname size M vector of residue names 
//' @param V output, matrix Nvx2 of vertices
//' @export
// [[Rcpp::export]]
IntegerMatrix vertex_ho2(IntegerVector XResno, StringVector XResname, IntegerVector YResno, StringVector YResname) {
  int i,ip,v=0;
  int N=XResno.size();
  int M=YResno.size();
  IntegerMatrix V(N*M,2);
  for(i=0; i<M; i++) {
    for (ip=0; ip<N; ip++) {
      if (strcmp(XResname(ip),YResname(i))==0) {
        V(v,0)=i+1;
        V(v,1)=ip+1;
        v++;
      }
    }
  }
  if (v<=1) {
    fprintf(stderr,"no product graph...\n");
    return V(Range(0,0),_);
  }
  return V(Range(0,v-1),_);
}

//' Constructs correspondence graph
//'
//' @param X matrix Nx3
//' @param Y matrix Nx3
//' @param V matrix Nvx2
//' @param deltadist double
//' @param mindist double
//' @param maxdist double
//' @export
// [[Rcpp::export]]
IntegerMatrix buildGraph(NumericMatrix X, NumericMatrix Y, IntegerMatrix V, double deltadist, double mindist, double maxdist) {
  int i,j, nv, e=0;
  double d,d1,d2;
  std::vector<int> Etmp;
  
  nv=V.nrow();
  for(i=0; i<nv-1; i++) {
    for (j=i+1; j<nv; j++) {
      //printf("%d %d %d %d\n",V(i,0),V(j,0),V(i,1),V(j,1));
      d1=distloc(Y,V(i,0),V(j,0)); // dist between atoms i and j in Y
      d2=distloc(X,V(i,1),V(j,1)); // dist between atoms i and j in X
      d=fabs(d1-d2);
      if (d<=deltadist &&  V(i,0)!=V(j,0) && V(i,1)!=V(j,1) && d1<=maxdist && d2<=maxdist) {
//     if (d<=deltadist && d1>=mindist && d2>=mindist && d1<=maxdist && d2<=maxdist) {
        Etmp.push_back(V(i,1));
        Etmp.push_back(V(j,1));
        Etmp.push_back(V(i,0));
        Etmp.push_back(V(j,0));
        e+=1;
      }
    }
  }
  int n=Etmp.size();
  IntegerMatrix E(e,4);
  for (i=0,j=0; i<n-3; i+=4,j++) {
    E(j,0)=Etmp[i];
    E(j,1)=Etmp[i+1];
    E(j,2)=Etmp[i+2];
    E(j,3)=Etmp[i+3];
  }
  Etmp.clear();
  return E;
}

//' Constructs correspondence graph
//'
//' @param X matrix Nx3
//' @param Y matrix Nx3
//' @param V matrix Nvx2
//' @param deltadist double
//' @param mindist double
//' @param maxdist double
//' @export
// [[Rcpp::export]]
IntegerMatrix buildGraph_ho(NumericMatrix X, IntegerVector XResno, NumericMatrix Y, IntegerVector YResno, IntegerMatrix V, double deltadist, double mindist, double maxdist, int maxgap) {
  int i,j, nv, e=0, gap;
  double d,d1,d2;
  std::vector<int> Etmp;
  
  nv=V.nrow();
  
  for(i=0; i<nv-1; i++) {
    for (j=i+1; j<nv; j++) {
      //printf("%d %d %d %d\n",V(i,0),V(j,0),V(i,1),V(j,1));
      d1=distloc(Y,V(i,0),V(j,0)); // dist between atoms i and j in Y
      d2=distloc(X,V(i,1),V(j,1)); // dist between atoms i and j in X
      d=fabs(d1-d2);
      //gap=gaploc(XResno,V(i,1),V(j,1));
      gap=gaploc(YResno,V(i,0),V(j,0));
      if (d<=deltadist &&  V(i,0)!=V(j,0) && V(i,1)!=V(j,1) && gap<=maxgap) {
//     if (d<=deltadist && d1>=mindist && d2>=mindist && d1<=maxdist && d2<=maxdist && gap<=maxgap) {
        Etmp.push_back(V(i,1));
        Etmp.push_back(V(j,1));
        Etmp.push_back(V(i,0));
        Etmp.push_back(V(j,0));
        e+=1;
      }
    }
  }
  int n=Etmp.size();
  IntegerMatrix E(e,4);
  for (i=0,j=0; i<n-3; i+=4,j++) {
    E(j,0)=Etmp[i];
    E(j,1)=Etmp[i+1];
    E(j,2)=Etmp[i+2];
    E(j,3)=Etmp[i+3];
  }
  Etmp.clear();
  return E;
}


//' Constructs correspondence graph
//'
//' @param X matrix Nx3
//' @param Y matrix Nx3
//' @param V matrix Nvx2
//' @param deltadist double
//' @param mindist double
//' @param maxdist double
//' @export
// [[Rcpp::export]]
IntegerMatrix buildGraph_ho2(NumericMatrix X, IntegerVector XResno, IntegerVector XResInd, NumericMatrix Y, IntegerVector YResno, IntegerVector YResInd, IntegerMatrix V, double deltadist, int maxgap) {
  int i,j,k,l,nv, e=0, gap;
  double d,d1,d2;
  std::vector<int> Etmp;
  
  nv=V.nrow();

  for(i=0; i<nv-1; i++) {
    for (j=i+1; j<nv; j++) {
      //printf("%d %d %d %d\n",V(i,0),V(j,0),V(i,1),V(j,1));
      d=0.0;
      for (k=0; k<YResInd(V(i,0))-YResInd(V(i,0)-1); k++) {
        for (l=0; l<XResInd(V(i,1))-XResInd(V(i,1)-1); l++) {
          d1=distloc(Y, YResInd(V(i,0)-1)+k, YResInd(V(j,0)-1)+l); // dist between residue Vi, atom k and Vj,k in Y
          d2=distloc(X, XResInd(V(i,1)-1)+k, XResInd(V(j,1)-1)+l); // dist between residue Vi, atom k and Vj,k in X
          d=fmax(d,fabs(d1-d2));
        }
      }
      //gap=gaploc(XResno,V(i,1),V(j,1));
      gap=gaploc(YResno, YResInd(V(i,0)-1), YResInd(V(j,0)-1));
      if (d<=deltadist &&  V(i,0)!=V(j,0) && V(i,1)!=V(j,1) && gap<=maxgap) {
        Etmp.push_back(V(i,1));
        Etmp.push_back(V(j,1));
        Etmp.push_back(V(i,0));
        Etmp.push_back(V(j,0));
        e+=1;
      }
    }
  }
  int n=Etmp.size();
  IntegerMatrix E(e,4);
  for (i=0,j=0; i<n-3; i+=4,j++) {
    E(j,0)=Etmp[i];
    E(j,1)=Etmp[i+1];
    E(j,2)=Etmp[i+2];
    E(j,3)=Etmp[i+3];
  }
  Etmp.clear();
  return E;
}



//' @export
// [[Rcpp::export]]
NumericVector mapping_dist_sum2(NumericMatrix X, IntegerVector I, IntegerVector CI, NumericMatrix Y, IntegerVector J, IntegerVector CJ, double thresh) {
  double  score, sc, d1, d2, deltadist;
  int i,j,n;

  n=I.size();
  //printf("mapping_dist_sum %d %d\n", I.size(), CI.size());
  if (thresh<=0.) thresh=15.; //patchsearch_NEIDIST
  NumericVector scores(n);
  for (i=0; i<I.size();i++) {
    score=0.;
    for (j=0; j<CI.size();j++) {
       if ((CJ(j)==J(i) && CI(j)!=I(i)) || (CI(j)==I(i) && CJ(j)!=J(i))) {
        score=0.;
        break;
      } else {
	        d1=distloc(Y,CJ(j),J(i));
	        d2=distloc(X,CI(j),I(i));
	        deltadist=fabs(d1-d2);
	        sc=1-deltadist/thresh;
      }
      //if (deltadist>thresh) sc=0.;
      score=score+sc;
    }
    //printf("scores[%d]=%lf\n", i, score);
    scores[i]=score;
  }
  return scores;
}

// all nodes in V connected to at least nbnei in C and not to far in target (X) sequence for hot spots application 
// [[Rcpp::export]]
IntegerVector selectLinks(IntegerVector C, IntegerMatrix V, NumericMatrix X, IntegerVector XResno, NumericMatrix Y, IntegerVector YResno, double thresh, int nbnei, int maxgap) {
  int i, j, inC, gap, count=0;
  int nC=C.size();
  IntegerVector K;
  IntegerVector I(nC), J(nC);
  int N=X.nrow();
  
  for (int i=0; i < nC; i++) {
    I(i)=(C(i)-1)%N+1;
    J(i)=(int)(C(i)-1)/N+1;
  }
  
  double deltadist, d1, d2;
  for (j=0; j<V.nrow(); j++) {
    count=0;
    inC=0;
    for (i=0; i<nC; i++) {
      if (J(i)!=V(j,0) || I(i)!=V(j,1)) {
        d1=distloc(Y,J(i),V(j,0));
        d2=distloc(X,I(i),V(j,1));
        deltadist=fabs(d1-d2);
        gap=gaploc(XResno, I(i), V(j,1));
          //deltadist=fabs(d1-d2)/std::min(d1,d2);
        if (deltadist<=thresh && gap<=maxgap) count++;
      }
      else inC=1; // j is in C
    }
    if (count>=nbnei && !inC) K.push_back((V(j,0)-1)*N+V(j,1));
  }
  return K;
}
