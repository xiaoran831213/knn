#include <RcppEigen.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace Eigen;
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
    // MatrixXf A[K];
    for(int i = 0; i < K; i++)
    {
	// A_i = sum_{j=1}^k lamda[i, j] (R V[i] R)
	ArrayXXf a = B[0].array() * s(i, 0);
	for(int j = 1; j < K; j++)
	    a += B[j].array() * s(i, j);

	// estimate the i th. variance component
	float w = y.dot(a.matrix() * y);
	// if(w < 0.0f) 		// modified MINQUE required
	// {
	//     // project A to its nearest in the PSD space
	//     SelfAdjointEigenSolver<MatrixXf> egn(((a + a.transpose()) / 2.0f).matrix());
	//     ArrayXXf            d(egn.eigenvalues());
	//     MatrixXf            v(egn.eigenvectors());
	//     cout << "d.egn=" << d << std::endl;
	//     d = d * (d > d.maxCoeff() * 1.49e-08f).cast<float>();
	//     for(int i = 0; i < K; i++)
	// 	v.col(i) *= d(i);
	//     a = v * egn.eigenvectors().transpose();
	//     w = y.dot(a.matrix() * y);
	// }
	// A[i] = a.matrix();
	W(i) = w;
    }
  
    List ret;
    ret["vcs"] = W;
    return ret;
}
