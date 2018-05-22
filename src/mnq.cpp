// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Cpp implementation of kernel MINQUE (kmq)
// only explicit kernels, no fixed effect, no implicit product kernels
// from random effect either.
SEXP mnq(SEXP _K, SEXP _p)
{
  List K=as<List>(_K);
  dvec p=as<dvec>(_p);

  dmat V = as<dmat>(K[0]);
  int  N = V.n_rows;
  for(int i=1; i<N; i++)
    {
      V += as<dmat>(K[i]);
    }
  dmat I(N, N, fill::eye);

  Rcout << "V=" << V << '\n';
  Rcout << "p=" << p << '\n';
  Rcout << "N=" << N << '\n';

  return(wrap(N));
}

// mnq
// SEXP mnq(List K, NumericVector x, NumericVector p);
// RcppExport SEXP _knn_mnq(SEXP KSEXP, SEXP xSEXP, SEXP pSEXP) {
// BEGIN_RCPP
//     Rcpp::RObject rcpp_result_gen;
//     Rcpp::RNGScope rcpp_rngScope_gen;
//     Rcpp::traits::input_parameter< List >::type K(KSEXP);
//     Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
//     Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
//     rcpp_result_gen = Rcpp::wrap(mnq(K, x, p));
//     return rcpp_result_gen;
// END_RCPP
// }

