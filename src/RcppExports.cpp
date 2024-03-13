// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// distloc
double distloc(NumericMatrix X, int ind1, int ind2);
RcppExport SEXP _pepit_distloc(SEXP XSEXP, SEXP ind1SEXP, SEXP ind2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type ind1(ind1SEXP);
    Rcpp::traits::input_parameter< int >::type ind2(ind2SEXP);
    rcpp_result_gen = Rcpp::wrap(distloc(X, ind1, ind2));
    return rcpp_result_gen;
END_RCPP
}
// gaploc
int gaploc(IntegerVector ResNo, int ind1, int ind2);
RcppExport SEXP _pepit_gaploc(SEXP ResNoSEXP, SEXP ind1SEXP, SEXP ind2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ResNo(ResNoSEXP);
    Rcpp::traits::input_parameter< int >::type ind1(ind1SEXP);
    Rcpp::traits::input_parameter< int >::type ind2(ind2SEXP);
    rcpp_result_gen = Rcpp::wrap(gaploc(ResNo, ind1, ind2));
    return rcpp_result_gen;
END_RCPP
}
// vertex
IntegerMatrix vertex(DataFrame XProp, DataFrame YProp, int mode, int size, int hse);
RcppExport SEXP _pepit_vertex(SEXP XPropSEXP, SEXP YPropSEXP, SEXP modeSEXP, SEXP sizeSEXP, SEXP hseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type XProp(XPropSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type YProp(YPropSEXP);
    Rcpp::traits::input_parameter< int >::type mode(modeSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< int >::type hse(hseSEXP);
    rcpp_result_gen = Rcpp::wrap(vertex(XProp, YProp, mode, size, hse));
    return rcpp_result_gen;
END_RCPP
}
// vertex2
IntegerMatrix vertex2(CharacterVector XProp, IntegerMatrix XHSE, CharacterVector YProp, IntegerMatrix YHSE);
RcppExport SEXP _pepit_vertex2(SEXP XPropSEXP, SEXP XHSESEXP, SEXP YPropSEXP, SEXP YHSESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type XProp(XPropSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type XHSE(XHSESEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type YProp(YPropSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type YHSE(YHSESEXP);
    rcpp_result_gen = Rcpp::wrap(vertex2(XProp, XHSE, YProp, YHSE));
    return rcpp_result_gen;
END_RCPP
}
// vertex_ho
IntegerMatrix vertex_ho(StringVector XProp, StringVector XResname, StringVector YProp, StringVector YResname);
RcppExport SEXP _pepit_vertex_ho(SEXP XPropSEXP, SEXP XResnameSEXP, SEXP YPropSEXP, SEXP YResnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type XProp(XPropSEXP);
    Rcpp::traits::input_parameter< StringVector >::type XResname(XResnameSEXP);
    Rcpp::traits::input_parameter< StringVector >::type YProp(YPropSEXP);
    Rcpp::traits::input_parameter< StringVector >::type YResname(YResnameSEXP);
    rcpp_result_gen = Rcpp::wrap(vertex_ho(XProp, XResname, YProp, YResname));
    return rcpp_result_gen;
END_RCPP
}
// vertex_ho2
IntegerMatrix vertex_ho2(IntegerVector XResno, StringVector XResname, IntegerVector YResno, StringVector YResname);
RcppExport SEXP _pepit_vertex_ho2(SEXP XResnoSEXP, SEXP XResnameSEXP, SEXP YResnoSEXP, SEXP YResnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type XResno(XResnoSEXP);
    Rcpp::traits::input_parameter< StringVector >::type XResname(XResnameSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type YResno(YResnoSEXP);
    Rcpp::traits::input_parameter< StringVector >::type YResname(YResnameSEXP);
    rcpp_result_gen = Rcpp::wrap(vertex_ho2(XResno, XResname, YResno, YResname));
    return rcpp_result_gen;
END_RCPP
}
// buildGraph
IntegerMatrix buildGraph(NumericMatrix X, NumericMatrix Y, IntegerMatrix V, double deltadist, double mindist, double maxdist);
RcppExport SEXP _pepit_buildGraph(SEXP XSEXP, SEXP YSEXP, SEXP VSEXP, SEXP deltadistSEXP, SEXP mindistSEXP, SEXP maxdistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type V(VSEXP);
    Rcpp::traits::input_parameter< double >::type deltadist(deltadistSEXP);
    Rcpp::traits::input_parameter< double >::type mindist(mindistSEXP);
    Rcpp::traits::input_parameter< double >::type maxdist(maxdistSEXP);
    rcpp_result_gen = Rcpp::wrap(buildGraph(X, Y, V, deltadist, mindist, maxdist));
    return rcpp_result_gen;
END_RCPP
}
// buildGraph_ho
IntegerMatrix buildGraph_ho(NumericMatrix X, IntegerVector XResno, NumericMatrix Y, IntegerVector YResno, IntegerMatrix V, double deltadist, double mindist, double maxdist, int maxgap);
RcppExport SEXP _pepit_buildGraph_ho(SEXP XSEXP, SEXP XResnoSEXP, SEXP YSEXP, SEXP YResnoSEXP, SEXP VSEXP, SEXP deltadistSEXP, SEXP mindistSEXP, SEXP maxdistSEXP, SEXP maxgapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type XResno(XResnoSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type YResno(YResnoSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type V(VSEXP);
    Rcpp::traits::input_parameter< double >::type deltadist(deltadistSEXP);
    Rcpp::traits::input_parameter< double >::type mindist(mindistSEXP);
    Rcpp::traits::input_parameter< double >::type maxdist(maxdistSEXP);
    Rcpp::traits::input_parameter< int >::type maxgap(maxgapSEXP);
    rcpp_result_gen = Rcpp::wrap(buildGraph_ho(X, XResno, Y, YResno, V, deltadist, mindist, maxdist, maxgap));
    return rcpp_result_gen;
END_RCPP
}
// buildGraph_ho2
IntegerMatrix buildGraph_ho2(NumericMatrix X, IntegerVector XResno, IntegerVector XResInd, NumericMatrix Y, IntegerVector YResno, IntegerVector YResInd, IntegerMatrix V, double deltadist, int maxgap);
RcppExport SEXP _pepit_buildGraph_ho2(SEXP XSEXP, SEXP XResnoSEXP, SEXP XResIndSEXP, SEXP YSEXP, SEXP YResnoSEXP, SEXP YResIndSEXP, SEXP VSEXP, SEXP deltadistSEXP, SEXP maxgapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type XResno(XResnoSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type XResInd(XResIndSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type YResno(YResnoSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type YResInd(YResIndSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type V(VSEXP);
    Rcpp::traits::input_parameter< double >::type deltadist(deltadistSEXP);
    Rcpp::traits::input_parameter< int >::type maxgap(maxgapSEXP);
    rcpp_result_gen = Rcpp::wrap(buildGraph_ho2(X, XResno, XResInd, Y, YResno, YResInd, V, deltadist, maxgap));
    return rcpp_result_gen;
END_RCPP
}
// mapping_dist_sum2
NumericVector mapping_dist_sum2(NumericMatrix X, IntegerVector I, IntegerVector CI, NumericMatrix Y, IntegerVector J, IntegerVector CJ, double thresh);
RcppExport SEXP _pepit_mapping_dist_sum2(SEXP XSEXP, SEXP ISEXP, SEXP CISEXP, SEXP YSEXP, SEXP JSEXP, SEXP CJSEXP, SEXP threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type I(ISEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type CI(CISEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type J(JSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type CJ(CJSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    rcpp_result_gen = Rcpp::wrap(mapping_dist_sum2(X, I, CI, Y, J, CJ, thresh));
    return rcpp_result_gen;
END_RCPP
}
// selectLinks
IntegerVector selectLinks(IntegerVector C, IntegerMatrix V, NumericMatrix X, IntegerVector XResno, NumericMatrix Y, IntegerVector YResno, double thresh, int nbnei, int maxgap);
RcppExport SEXP _pepit_selectLinks(SEXP CSEXP, SEXP VSEXP, SEXP XSEXP, SEXP XResnoSEXP, SEXP YSEXP, SEXP YResnoSEXP, SEXP threshSEXP, SEXP nbneiSEXP, SEXP maxgapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type V(VSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type XResno(XResnoSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type YResno(YResnoSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< int >::type nbnei(nbneiSEXP);
    Rcpp::traits::input_parameter< int >::type maxgap(maxgapSEXP);
    rcpp_result_gen = Rcpp::wrap(selectLinks(C, V, X, XResno, Y, YResno, thresh, nbnei, maxgap));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pepit_distloc", (DL_FUNC) &_pepit_distloc, 3},
    {"_pepit_gaploc", (DL_FUNC) &_pepit_gaploc, 3},
    {"_pepit_vertex", (DL_FUNC) &_pepit_vertex, 5},
    {"_pepit_vertex2", (DL_FUNC) &_pepit_vertex2, 4},
    {"_pepit_vertex_ho", (DL_FUNC) &_pepit_vertex_ho, 4},
    {"_pepit_vertex_ho2", (DL_FUNC) &_pepit_vertex_ho2, 4},
    {"_pepit_buildGraph", (DL_FUNC) &_pepit_buildGraph, 6},
    {"_pepit_buildGraph_ho", (DL_FUNC) &_pepit_buildGraph_ho, 9},
    {"_pepit_buildGraph_ho2", (DL_FUNC) &_pepit_buildGraph_ho2, 9},
    {"_pepit_mapping_dist_sum2", (DL_FUNC) &_pepit_mapping_dist_sum2, 7},
    {"_pepit_selectLinks", (DL_FUNC) &_pepit_selectLinks, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_pepit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
