#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// knl_mnq
RcppExport SEXP knl_mnq(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] =
{
    {"knl_mnq", (DL_FUNC) &knl_mnq, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_knn(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
