#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#if !defined(EIGEN_USE_MKL) // Don't use R's Lapack.h if MKL is enabled
#include <R_ext/Lapack.h>
#endif

#include <numeric>
#include <math.h>

// [[Rcpp::plugins(cpp11)]]

using Eigen::Lower;
using Eigen::MatrixXd;

using Rcpp::as;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::wrap;

//' @export
// [[Rcpp::export]]
MatrixXd ibs(NumericMatrix X) {
  MatrixXd G(as<MatrixXd>(X));
  G *= 0.5;
  int r(G.rows()), c(G.cols());
  
  MatrixXd GGt(MatrixXd(r, r).setZero().selfadjointView<Lower>().rankUpdate(G));
  MatrixXd G2(MatrixXd::Ones(r, r) - G);
  MatrixXd G2G2t(MatrixXd(r, r).setZero().selfadjointView<Lower>().rankUpdate(G2));
  
  // Calculate identity-by-state
  return (1/double(c))*(GGt + G2G2t);
}
