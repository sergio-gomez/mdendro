// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcppLinkage
Rcpp::List rcppLinkage(const Rcpp::NumericVector& prox, bool isDistance, int digits, std::string method, double methodPar, bool isWeighted, bool isVariable);
RcppExport SEXP _mdendro_rcppLinkage(SEXP proxSEXP, SEXP isDistanceSEXP, SEXP digitsSEXP, SEXP methodSEXP, SEXP methodParSEXP, SEXP isWeightedSEXP, SEXP isVariableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type prox(proxSEXP);
    Rcpp::traits::input_parameter< bool >::type isDistance(isDistanceSEXP);
    Rcpp::traits::input_parameter< int >::type digits(digitsSEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type methodPar(methodParSEXP);
    Rcpp::traits::input_parameter< bool >::type isWeighted(isWeightedSEXP);
    Rcpp::traits::input_parameter< bool >::type isVariable(isVariableSEXP);
    rcpp_result_gen = Rcpp::wrap(rcppLinkage(prox, isDistance, digits, method, methodPar, isWeighted, isVariable));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mdendro_rcppLinkage", (DL_FUNC) &_mdendro_rcppLinkage, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_mdendro(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}