#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// knl_mnq
RcppExport SEXP knl_mnq(SEXP, SEXP);
RcppExport SEXP knl_mn2(SEXP, SEXP);
RcppExport SEXP egn_mnq(SEXP, SEXP);
RcppExport SEXP vcm_dv1(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] =
{
    {"knl_mnq", (DL_FUNC) &knl_mnq, 2},
    // {"knl_mn2", (DL_FUNC) &knl_mn2, 2},
    {"egn_mnq", (DL_FUNC) &egn_mnq, 2},
    {"vcm_dv1", (DL_FUNC) &vcm_dv1, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_knn(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
