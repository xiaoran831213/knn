#include "RcppArmadillo.h"
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
RcppExport SEXP knl_mnq(SEXP _y, SEXP _V, SEXP _psd)
{
    // // [[Rcpp::export]]
    int psd = as<int>(_psd);
    
    fmat y = as<fmat>(_y);
    int N = y.n_rows;		// sample size

    /* list of kernels */
    List _l = as<List>(_V);	// list of kernels
    int k = _l.length();	// number of kernels
    fmat V[k];
    for(int i = 0; i < k; i++)
	V[i] = as<fmat>(_l[i]);

    // Z: sum of kernels (Rao. 1971)
    fmat Z(N, N, fill::zeros);
    for(int i = 0; i < k; i++)
	Z += V[i];

    /* R = Z^{-1} (I - P_v), (Rao. 1971) */
    /* in case of no fixed effect X, project P_v = 0, R = Z^{-1} */
    fmat R = pinv(Z);

    /* W_i = R V_i R, (Rao. 1971) */
    fmat B[k];
    for(int i = 0; i < k; i++)
	B[i] = R * V[i] * R;

    /* The k-k matrix S, S_ij = Tr(B_i V_j) */
    fmat S(k, k);
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
    // fmat L = pinv(S);
    fmat L = inv_sympd(S);
    // List A(L.n_rows);	// k A-matrices
    fmat A[L.n_rows];
    fvec W(L.n_rows);		// k contrasts -> k VC(s)
    fvec d(N);			// eigen values of A
    fmat v(N, N);		// eigen vectors of A
    fmat u(N, N);		// eigen vectors of A
    for(int i = 0; i < k; i++)
    {
        // A_i = sum_{j=1}^k lamda[i, j] (R V[i] R)
	fmat a = B[0] * L(i, 0);
	for(int j = 1; j < L.n_cols; j++) // L.n_cols == k
	    a += B[j] * L(i, j);

	// estimate the i th. variance component
	float w = as_scalar(y.t() * a * y);
	if(w < 0.0f && psd)	// modified MINQUE required
	{
	    // project A to its nearest in the PSD space
	    eig_sym(d, v, (a + a.t())/2.0f);
	    for(int j = v.n_rows - 1; j >= 0; j--)
	    {
		if(d[j] > 0.0f)
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
    // fmat C = V[0] * W[0];	// marginal cov of y
    // fvec e(k);
    // for(int i = 1; i < k; i++)
    //  	C += V[i] * W[i];

    /* se(s2_i) = 2Tr(A_i C_y A_i C_y) */
    // for(int i = 0; i < k; i++)
    //  	e[i] = 2.0f * accu(square(A[i] * C));

    List ret;
    // ret["S"] = S;		// S_ij = Tr[B_i V_j]
    // ret["L"] = L;		// lambda for each contrast
    // ret["A"] = A;		// A matrices
    // ret["C"] = C;
    ret["vcs"] = W;		// variance components
    ret["se2"] = 1;		// standard errors
    return(Rcpp::wrap(ret));
}

