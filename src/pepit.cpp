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
//' @param XProp data frame  N nrows
//' @param YProp data frame  M nrows
//' @param mode 0 atom types comparison, 1 all features' comparison
//' @param size estimated initial size of V
//' @param V output matrix Nvx2 of vertices
//' @description data.frame(eleno, elety, resid, chain, resno, insert, x, y, z, pc, hseu, hsed)
//' @export
// [[Rcpp::export]]
IntegerMatrix vertex(DataFrame XProp,  DataFrame YProp, int mode, int size, int hse) {
  int i,ip,v = 0;
  int N = XProp.nrows();
  int M = YProp.nrows();
  if (size == 0) size = N*M;
  
  CharacterVector Xelety = XProp["elety"];
  CharacterVector Yelety = YProp["elety"];
  
  IntegerMatrix V(size,2);
  if (mode == 0) {
    for(i=0; i<M; i++) {
      for (ip=0; ip<N; ip++) {
        if (strcmp(Xelety(ip), Yelety(i))==0) {
          V(v,0)=i+1;
          V(v,1)=ip+1;
          v++;
        }
      }
    }
  }
  if (mode != 0) {
    IntegerVector XHSEu = XProp["hseu"];
    IntegerVector XHSEd = XProp["hsed"];
    IntegerVector Xpc = XProp["pc"];
    IntegerVector YHSEu = YProp["hseu"];
    IntegerVector YHSEd = YProp["hsed"];
    IntegerVector Ypc = YProp["pc"];
    for(i=0; i<M; i++) {
      for (ip=0; ip<N; ip++) {
        if (strcmp(Xelety(ip), Yelety(i)) == 0) {
          if (Xpc(ip) == Ypc(i)) {
            if ((XHSEu(ip) <= hse && YHSEu(i) <= hse) || (XHSEd(ip) <= hse && YHSEd(i) <= hse)) {
              V(v,0)=i+1;
              V(v,1)=ip+1;
              v++;
            }
          }
        }
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


//' @export
// [[Rcpp::export]]
NumericVector lDDT(NumericMatrix X, IntegerVector I, NumericMatrix Y, IntegerVector J) {
  double d1, d2, deltadist;
  int i, j, n, total, count0, count1, count2, count4;
  NumericVector out(1);
  
  n=I.size();
  NumericVector scores(n);
  total = 0;
  count0 = count1 = count2 = count4 = 0;
  for (i=0; i<I.size();i++) {
    for (j=0; j<I.size();j++) {
      d2=distloc(X,I(j),I(i));
      if (i != j && d2 < 15.0) {
	total++;
	d1=distloc(Y,J(j),J(i));
	deltadist=fabs(d1-d2);
	if (deltadist < 0.5) count0++;
	if (deltadist < 1.0) count1++;
	if (deltadist < 2.0) count2++;
	if (deltadist < 4) count4++;
      }
    }
  }

  if (total == 0)
    return 0.;

  out(0) = 0.25 * (double)(count0 + count1 + count2 + count4) / total;
  return out;
}


//' @export
// [[Rcpp::export]]
NumericVector mapping_dist_sum2(NumericMatrix X, IntegerVector I, IntegerVector CI, NumericMatrix Y, IntegerVector J, IntegerVector CJ, double thresh) {
  double  score, sc, d1, d2, deltadist;
  int i,j,n;

  //printf("X=%ld, Y=%ld, CI=%ld, CJ=%ld\n", X.size(), Y.size(), I.size(), J.size());
  n=I.size();
  if (thresh<=0.) {
    thresh=15.; //patchsearch_NEIDIST
    printf("thresh=%lf\n", thresh);
  }
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
	        sc=1.0 - deltadist/thresh;
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
