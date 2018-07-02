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
    dvec W(L.n_rows);		      // k contrasts -> k VC(s)
    dvec d(N);			      // eigen values of A
    dmat v(N, N);		      // eigen vectors of A
    dmat u(N, N);		      // eigen vectors of A
    for(int i = 0; i < k; i++)
    {
        // A = sum_{i=1}^k lamda[i, ] (R V[i, ] R)
	dmat a = B[0] * L(i, 0);
	for(int j = 1; j < L.n_cols; j++) // L.n_cols == k
	    a += B[j] * L(i, j);

	// estimate the i th. variance component
	double w = as_scalar(y.t() * a * y);
	if(w < 0) 		// modified MINQUE required
	{
	    // project A to its nearest in the PSD space
	    eig_sym(d, v, (a + a.t())/2);
	    for(int j = v.n_rows - 1; j >=0; j--)
	    {
		if(d[j] > 0)
		    u.col(j) = v.col(j) * d[j];
		else
		    u.col(j).zeros();
	    }
	    a = u * v.t();
	    w = as_scalar(y.t() * a * y);
	}
	A[i] = a;
	W[i] = w;
    }

    /* Standard Error of VC(s), assuming y is normal */
    dmat C = V[0] * W[0];	// marginal cov of y
    dvec e(k);
    for(int i = 1; i < k; i++)
	C += V[i] * W[i];

    // se(s2_i) = 2Tr(A_i C_y A_i C_y)
    for(int i = 0; i < k; i++)
	e[i] = 2.0 * accu(pow(as<dmat>(A[i]) * C, 2));

    List ret;
    ret["S"] = S;		// S_ij = Tr[B_i V_j]
    ret["L"] = L;		// lambda for each contrast
    ret["A"] = A;		// A matrices
    ret["C"] = C;
    ret["vcs"] = W;		// variance components
    ret["se2"] = e;		// standard errors
    return(Rcpp::wrap(ret));
}
