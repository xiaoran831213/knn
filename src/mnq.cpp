// #include <Rcpp.h>
#include "RcppArmadillo.h"
#include <vector>
#include <numeric>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

//' Kernel MINQUE
//'
//' Only explicit kernels without fixed effect, no implicit product kernels
//' from random effect either.
//' 
//' @param V a list of k kernel matrices matched to k variance components.
//' @param P a m-k matrix, each row is a contrast of k variance components.
//' @export
RcppExport SEXP knl_mnq(SEXP _y, SEXP _V)
{
    // // [[Rcpp::export]]
    dmat y = as<dmat>(_y);
    int N = y.n_rows;		// sample size

    /* list of kernels */
    List _l = as<List>(_V);	// list of kernels
    int k = _l.length();	// number of kernels
    std::vector<dmat> V(k);
    for(int i = 0; i < k; i++)
	V[i] = as<dmat>(_l[i]);

    // Z: sum of kernels (Rao. 1971)
    dmat Z = std::accumulate(V.begin() + 1, V.end(), V.front());

    /* R = Z^{-1} (I - P_v), (Rao. 1971) */
    /* in case of no fixed effect X, project P_v = 0, R = Z^{-1} */
    dmat R = inv_sympd(Z);

    /* W_i = R V_i R, (Rao. 1971) */
    std::vector<dmat> B(k);
    for(int i = 0; i < k; i++)
	B[i] = R * V[i] * R;

    /* The k-k matrix S, S_ij = Tr(B_i V_j) */
    dmat S(k, k);
    for(int i = 0; i < k; i++)
    {
	for(int j = i + 1; j < k; j++)
	{
	    S(i, j) = accu(B[i] % V[j]);
	    S(j, i) = S(i, j);
	}
	S(i, i) = accu(B[i] % V[i]);
    }

    /* Lambda_i = P_i' S^{-1} (Rao. 1917), i = 1 .. nrow(P).
       To solve each VC, the contrast matrix must be I(k), so
       Lambda = S^{-1} */
    dmat L = pinv(S);
    List A(L.n_rows);		      // k A-matrices
    dvec f(L.n_rows);		      // k contrasts -> k VC(s)
    // dmat y1(y);
    for(int i = 0; i < L.n_rows; i++) // L.n_rows == P.n_rows
    {
	// A = sum_{i=1}^k lamda[i, ] (R V[i, ] R)
	dmat a(N, N, fill::zeros);
	for(int j = 0; j < L.n_cols; j++) // L.n_cols == k
	    a += B[j] * L(i, j);
	A[i] = a;

	// estimate linear function of variance components
	// modified MINQUE, project A to PSD space
	dvec d;
	dmat v;
	eig_sym(d, v, (a + a.t())/2);
	d = d % (d > 0);
	// A_hat_i = v I(d) v.t()
	// sum(d * crossprod(v, y)^2) in R
	// f[i] = as_scalar(y.t() * v *  diagmat(d) * v.t() * y);
	f[i] = as_scalar(square(y.t() * v) * d);
    }

    List ret;
    ret["S"] = S;		// S_ij = Tr[B_i V_j]
    ret["L"] = L;		// lambda for each contrast
    ret["A"] = A;		// A matrices
    ret["f"] = f;		// contrasts
    return(Rcpp::wrap(ret));
}
