// #include <Rcpp.h>
#include "RcppArmadillo.h"
#include <vector>
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
// [[Rcpp::export]]
List knl_mnq(arma::dmat y, List V, arma::dmat P, arma::dvec psd)
{
    int k = V.length();		   // number of kernels
    int N = as<dmat>(V[0]).n_rows; // sample size

    // Z: sum of kernels (Rao. 1971)
    dmat Z = as<dmat>(V[0]);
    for(int i = 1; i < k; i++)
	Z += as<dmat>(V[i]);

    /* R = Z^{-1} (I - P_v), (Rao. 1971) */
    /* in case of no fixed effect X, project P_v = 0, R = Z^{-1} */
    dmat R = inv_sympd(Z);

    /* W_i = R V_i R, (Rao. 1971) */
    std::vector<dmat> B(k);
    for(int i = 0; i < k; i++)
    	B[i] = R * as<dmat>(V[i]) * R;

    /* The k-k matrix S, S_ij = Tr(B_i V_j) */
    dmat S(k, k);
    for(int i = 0; i < k; i++)
    {
	for(int j = i + 1; j < k; j++)
	{
	    S(i, j) = accu(B[i] % as<dmat>(V[j]));
	    S(j, i) = S(i, j);
	}
	S(i, i) = accu(B[i] % as<dmat>(V[i]));
    }

    /* Lambda_i = P_i' S^{-1} (Rao. 1917), i = 1 .. nrow(P). */
    dmat s = pinv(S);
    dmat L = P * s;

    List A(L.n_rows);		      // k A-matrices
    dvec f(L.n_rows);		      // k linear contrast outcome
    // dmat y1(y);
    for(int i = 0; i < L.n_rows; i++) // L.n_rows == P.n_rows
    {
	// A = sum_{i=1}^k lamda[i, ] (R V[i, ] R)
    	dmat a(N, N, fill::zeros);
    	for(int j = 0; j < L.n_cols; j++) // L.n_cols == k
    	    a += B[j] * L(i, j);
	A[i] = a;

	// estimate linear function of variance components
	if(as_scalar(psd) != 0)		// PSD MINQUE (modifiec MINQUE)
	{
	    // printf("PSD\n");
	    dvec d;
	    dmat v;
	    eig_sym(d, v, (a + a.t())/2);
	    d = d % (d > 0);
	    // sum(d * crossprod(v, y)^2)
	    // f[i] = as_scalar(y.t() * v *  diagmat(d) * v.t() * y);
	    f[i] = as_scalar(square(y.t() * v) * d);
	}
	else			// RAW MINQUE
	{
	    // printf("RAW\n");
	    f[i] = as_scalar(y.t() * a * y);
	}
    }

    List ret;
    ret["S"] = S;		// S_ij = Tr[B_i V_j]
    ret["s"] = s;
    // ret["L"] = L;		// lambda for each contrast
    // ret["A"] = A;		// A matrices
    ret["f"] = f;		// contrasts
    return(ret);
}
