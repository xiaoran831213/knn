#include <Rcpp.h>
#include <RcppEigen.h>
using Rcpp::as;
using Rcpp::List;
using Eigen::MatrixXf;
using Eigen::ArrayXXf;
using Eigen::VectorXf;
using Eigen::Map;
using Eigen::LLT;
using Eigen::BDCSVD;
using Eigen::JacobiSVD;
using Eigen::DecompositionOptions::ComputeThinU;
using Eigen::DecompositionOptions::ComputeThinV;
using Eigen::DecompositionOptions::ComputeFullU;
using Eigen::DecompositionOptions::ComputeFullV;
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

    /* list of kernels */
    const List __V(as<List>(_V));
    const int  K(__V.length());
    MatrixXf   V[K];
    for(int i = 0; i < K; i++)
	V[i] = as<MatrixXf>(__V[i]);

    // sum of kernels (Rao. 1971)
    MatrixXf Z(V[0]);
    for(int i = 1; i < K; i++)
	Z += V[i];

    /* R = Z^{-1} (I - P_v), (Rao. 1971) */
    /* in case of no fixed effect X, project P_v = 0, R = Z^{-1} */
    LLT<MatrixXf> L = Z.llt();

    /* W_i = R V_i R, (Rao. 1971) */
    MatrixXf D[K];
    MatrixXf B[K];
    for(int i = 0; i < K; i++)
    {
	D[i] = L.solve(V[i]);	// R V_i
	B[i] = L.solve(D[i]);	// R V_i R
    }

    /* The k-k matrix S, S_ij = Tr[(R V_i) (R V_j)] */
    MatrixXf S(K, K);
    for(int i = 0; i < K; i++)
    {
  	for(int j = i + 1; j < K; j++)
  	{
  	    S(i, j) = D[i].cwiseProduct(D[j]).sum();
  	    S(j, i) = S(i, j);
  	}
  	S(i, i) = D[i].array().pow(2).sum();
    }

    /* Lambda_i = P_i' S^{-1} (Rao. 1917), i = 1 .. nrow(P). To solve
       each VC, the contrast matrix must be I(k), so Lambda = S^{-1} */
    MatrixXf s = S.llt().solve(MatrixXf::Identity(K, K));
    VectorXf W(K);
    MatrixXf A[K];
    // fvec d(N);	              // eigen values of A
    // fmat v(N, N);		      // eigen vectors of A
    // fmat u(N, N);		      // eigen vectors of A
    for(int i = 0; i < K; i++)
    {
	// A_i = sum_{j=1}^k lamda[i, j] (R V[i] R)
	ArrayXXf a = B[0].array() * s(i, 0);
	for(int j = 1; j < K; j++)
	    a += B[j].array() * s(i, j);

	// estimate the i th. variance component
	float w = y.dot(a.matrix() * y);
	if(w < 0.0f) 		// modified MINQUE required
	{
	    // project A to its nearest in the PSD space
	    // eig_sym(d, v, (a + a.t())/2.0f);
	    JacobiSVD<MatrixXf> svd(((a + a.transpose()) / 2.0f).matrix(), ComputeFullU|ComputeFullV);
	    ArrayXXf            d(svd.singularValues());
	    cout << "d=" << d << std::endl;
	    d = d * (d > d(0) * svd.threshold()).cast<float>();
	    // cout << "U=" << svd.matrixU() << std::endl;
	    cout << "d=" << d << std::endl;
	    // cout << "V=" << svd.matrixV() << std::endl;
	    a = svd.matrixU() * d.matrix().asDiagonal() * svd.matrixV().transpose();
	    w = y.dot(a.matrix() * y);
	    // for (int l = 0; l < N; l++)
	    // {
	    // 	if(d(l) < thd)
	    // 	    d(l) = 0.0f;
	    // }
	    // for(int j = v.n_rows - 1; j >= 0; j--)
	    // {
	    //     if(d[j] > 0.0f)
	    // 	       u.col(j) = v.col(j) * d[j];
	    // 	   else
	    // 	       u.col(j).zeros();
	    // }
	    // a = u * v.t();
	    // w = as_scalar(y.t() * a * y);
	}
	A[i] = a.matrix();
	W(i) = w;
    }
  
    List ret;
    ret["S"] = S;
    ret["vcs"] = W;
    return ret;
}
