// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <vector>
using namespace Rcpp;
using namespace arma;

//' Kernel MINQUE
//'
//' Only explicit kernels, no fixed effect, no implicit product kernels
//' from random effect either.
//' Multiply a number by two
//' 
//' @param v r-list, of the kernels
//' @param P matrix, of each column is a contrast of VCs
//' @export
// [[Rcpp::export]]
List mnq(List v, arma::dmat p)
{
    int k = v.length();		   // number of kernels
    int N = as<dmat>(v[0]).n_rows; // sample size
    dmat V = as<dmat>(v[0]);

    // V: sum of kernels (Rao. 1971)
    for(int i = 1; i < k; i++)
	V += as<dmat>(v[i]);

    /* R = V^{-1} (I - P_v), (Rao. 1971)*/
    /* in case of no fixed effect X, project P_v = 0, R = V^{-1} */
    dmat R = inv_sympd(V);

    /* W_i = R v_i R, (Rao. 1971) */
    List W(k);
    for(int i = 0; i < k; i++)
    	W[i] = R * as<dmat>(v[i]) * R;

    printf("done with W=R v_i R\n");
    /* The k-k matrix S, S_ij = Tr(W_i v_j) */
    dmat S(k, k);
    for(int i = 0; i < k; i++)
    {
	for(int j = i + 1; j < k; j++)
	{
	    S(i, j) = accu(as<dmat>(W[i]) % as<dmat>(v[j]));
	    S(j, i) = S(i, j);
	}
	S(i, i) = accu(as<dmat>(W[i]) % as<dmat>(v[i]));
    }
    /* Lambda_i = p_i' S^{-1} (Rao. 1917), i = 1 .. nrow(P). */
    dmat L = p * pinv(S);

    // List lsA(lmx.n_cols);
    // for(int j = 0; j < lmx.n_cols; j++)
    // {
    // 	dmat x(N, N, fill::zeros);
    // 	for(int i = 0; i < lmx.n_rows; i++) // assert: lmx.n_rows==k
    // 	    x += as<dmat>(W[i]) * lmx(i, j);
    // 	lsA[j] = x;
    // }
    
    List ret;
    ret["dim"] = ivec{k, N};
    ret["v"] = v;
    ret["p"] = p;
    ret["V"] = V;		// sum of v_i, i = 1 .. k
    ret["R"] = R;		// V^{-1} (I - P_v)
    ret["W"] = W;
    ret["S"] = S;		// S_ij = Tr[W_i V_j]
    ret["L"] = L;		// lambda for each contrast
    return(ret);
}
