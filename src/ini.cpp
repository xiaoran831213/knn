// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// MINQUE by Cpp
RcppExport SEXP mnq(SEXP _K, SEXP _p);

static const R_CallMethodDef CallEntries[] = {
    {"_cpp_mnq", (DL_FUNC) &mnq, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_knn(DllInfo *dll)
{
    printf("R_init_mylib\n");
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

RcppExport void R_unload_mylib(DllInfo *info)
{
    printf("R_unload_mylib\n");
    /* Release resources. */
}
