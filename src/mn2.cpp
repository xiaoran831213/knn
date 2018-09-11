#include <Rcpp.h>
#include <RcppEigen.h>
using Rcpp::as;
using Rcpp::List;
using Eigen::MatrixXd;
using Eigen::MatrixXf;
using Eigen::VectorXf;
using std::cout;

//' Eigen3 MINQUE
//'
//' Only explicit kernels without fixed effect, no implicit product kernels
//' from random effect either.
//' 
//' @param V a list of k kernel matrices matched to k variance components.
//' @param P a m-k matrix, each row is a contrast of k variance components.
RcppExport SEXP egn_mnq(SEXP _y, SEXP _V)
{
    const VectorXf y(as<VectorXf>(_y));
    const int      N(y.rows());	    // sample size

    cout << "y=" << y << "\n";
    cout << "N=" << N << "\n";

    /* list of kernels */
    const List __V(as<List>(_V));
    const int  K(__V.length());
    MatrixXf   V[K];
    for(int i = 0; i < K; i++)
	V[i] = as<MatrixXf>(__V[i]);

    // Z: sum of kernels (Rao. 1971)
    MatrixXf Z(V[0]);
    for(int i = 1; i < K; i++)
	Z += V[i];
    cout << std::setw(3) << "Z=" << Z << "\n";

    /* R = Z^{-1} (I - P_v), (Rao. 1971) */
    /* in case of no fixed effect X, project P_v = 0, R = Z^{-1} */
    // MatrixXf R = inv_sympd(Z);

    /* W_i = R V_i R, (Rao. 1971) */
    // fmat B[k];
    // for(int i = 0; i < k; i++)
    // 	B[i] = R * V[i] * R;

    /* The k-k matrix S, S_ij = Tr(B_i V_j) */
    // fmat S(k, k);
    // for(int i = 0; i < k; i++)
    // {
    // 	for(int j = i + 1; j < k; j++)
    // 	{
    // 	    S(i, j) = accu(B[i] % V[j]);
    // 	    S(j, i) = S(i, j);
    // 	}
    // 	S(i, i) = accu(B[i] % V[i]);
    // }

    return _y;
}
