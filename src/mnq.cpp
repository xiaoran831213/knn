// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <vector>
using namespace Rcpp;
using namespace arma;

// Cpp implementation of kernel MINQUE (kmq)
// only explicit kernels, no fixed effect, no implicit product kernels
// from random effect either.
RcppExport SEXP mnq(SEXP _K, SEXP _C)
{
    List K=as<List>(_K);	// kernels / MCVs (matrics of covariance)
    dmat C=as<dmat>(_C);	// contrast matrix for variance components

    dmat V = as<dmat>(K[0]);	// sum of kernels
    const int L = K.length();	// number of kernels
    for(int i=1; i<L; i++)
    {
	V += as<dmat>(K[i]);
    }

    int  N = V.n_rows;		// sample size
    dmat I(N, N, fill::eye);	// identity kernel (matrix of white noise)

    printf("N=%d, L=%d\n", N, L);

    /* Inverse of V -- the sum of kernels:*/
    dmat W = inv_sympd(V);

    /* Cache the result: B_i = inv(V) K_i inv(V) for i = 1 .. L */
    List B(L);
    printf("B.length()=%d\n", B.length());
    for(int i = 0; i < L; i++)
    {
    	B[i] = W * as<dmat>(K[i]) * W;
    }
    printf("Done with B\n");

    /* The L-L matrix S */
    dmat S(L, L);
    for(int i = 0; i < L; i++)
    {
	for(int j=i + 1; j < L; j++)
	{
	    S(i, j) = accu(as<dmat>(B[i]) % as<dmat>(K[j]));
	    S(j, i) = S(i, j);
	}
	S(i, i) = accu(as<dmat>(B[i]) % as<dmat>(K[i]));
    }
    dmat T = pinv(S);		// Moore-Penrose pseudo-inverse of S
    dmat R = T * C;		// lambdas

    List A(R.n_cols);
    for(int j = 0; j < R.n_cols; j++)
    {
	dmat x(N, N, fill::zeros);
	for(int i = 0; i < R.n_rows; i++) // assertion: R.n_rows==L
	    x += as<dmat>(B[i]); // * R(i, j);
	A[j] = x;
    }
    
    List ret;
    ret["B"] = B;
    ret["K"] = K;
    ret["N"] = N;
    ret["L"] = L;
    ret["C"] = C;
    ret["V"] = V;
    ret["W"] = W;
    ret["S"] = S;
    ret["T"] = T;
    ret["R"] = R;
    ret["A"] = A;

    return(wrap(ret));
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

