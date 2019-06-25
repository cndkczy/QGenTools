// Calculates the realized relationship method using the method proposed by
// Astle and Balding (2009) doi: 10.1214/09-STS307

#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#if !defined(EIGEN_USE_MKL) // Don't use R's Lapack.h if MKL is enabled
#include <R_ext/Lapack.h>
#endif

#include <numeric>
#include <math.h>

// [[Rcpp::plugins(cpp11)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;

using Rcpp::as;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::wrap;

//' @export
// [[Rcpp::export]]
NumericMatrix realized_ab(NumericMatrix X) {
  MatrixXd G(as<MatrixXd>(X));
  int r(G.rows()), c(G.cols());
  
  // Calculate allele frequencies
  VectorXd freqs(G.colwise().sum()/double(r));
  VectorXd sd((0.5*freqs.array()*(2.0 - freqs.array())).sqrt());
  
  // Center and standardize the marker matrix
  for (size_t i = 0; i < c; ++i)
    G.col(i) = (G.col(i).array() - freqs(i))/sd(i);
  
  // Form the cross product
  MatrixXd GGt(G*G.transpose());
  
  // Realized AB matrix
  NumericMatrix AB(wrap((1/double(c))*GGt));
  rownames(AB) = rownames(X);
  colnames(AB) = rownames(X);
  
  return AB;
}
