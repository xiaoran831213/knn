// #include <Rcpp.h>
#include "RcppArmadillo.h"
#include <vector>
#include <numeric>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

//' Gradient over variance components
//'
//'
//' @param W a L-Q matrix, Q sets of variance components;
//' @param K a list of L kernel matrices;
//' @param Y a N-Q matrix of Q responses;
RcppExport SEXP vcm_dv1(SEXP _W, SEXP _K, SEXP _Y)
{
    dmat W = exp(as<dmat>(_W));	// weights - vcs
    dmat Y = as<dmat>(_Y);	// responses
    int M = Y.n_cols;		// number of responses
    int L = W.n_rows;		// number of kernels

    /* kernels */
    std::vector<dmat> K(L);
    for(int i = 0; i < L; i++)
	K[i] = as<dmat>(as<List>(_K)[i]);

    /* covariance */
    std::vector<dmat> C(M);
    for(int j = 0; j < M; j++)
    {
	C[j] = W(0, j) * K[0];
	for(int i = 1; i < L; i++)
	    C[j] += W(i, j) * K[i];
    }

    /* derivative */
    dmat G = zeros(W.n_rows, W.n_cols);
    for(int j = 0; j < M; j++)
    {
	dmat A = inv_sympd(C[j]);
	dvec b = A * Y.col(j);
	dmat T = b * b.t() - A;

	for(int i = 0; i < L; i++)
	    G(i, j) = -.5 * accu(T % K[i]) * W(i, j);
    }

    return(Rcpp::wrap(G));
}
