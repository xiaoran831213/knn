#include <Rcpp.h>
#include <RcppEigen.h>
using Rcpp::as;
using Eigen::MatrixXd;

//' Eigen3 MINQUE
//'
//' Only explicit kernels without fixed effect, no implicit product kernels
//' from random effect either.
//' 
//' @param V a list of k kernel matrices matched to k variance components.
//' @param P a m-k matrix, each row is a contrast of k variance components.
SEXP egn_mnq(SEXP _y, SEXP _V)
{
    return _y;
}
