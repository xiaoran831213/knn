#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP vcm_dv1(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] =
{
    {"vcm_dv1", (DL_FUNC) &vcm_dv1, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_knn(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
