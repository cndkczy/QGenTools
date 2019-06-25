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
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using Rcpp::as;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::wrap;

//' @export
// [[Rcpp::export]]
double hudson_fst(NumericMatrix X1, NumericMatrix X2) {
  // Map the genotype matrices to Eigen matrices
  const Map<MatrixXd> G1(as<Map<MatrixXd> >(X1));
  const Map<MatrixXd> G2(as<Map<MatrixXd> >(X2));
  
  int n1(G1.rows()), n2(G2.rows());
  int m1(G1.cols()), m2(G2.cols());
  
  // Calculate allele frequencies
  VectorXd p1(G1.colwise().sum()/double(2*n1));
  VectorXd p2(G2.colwise().sum()/double(2*n2));
  VectorXd q1((p1.array() - 1.0).abs());
  VectorXd q2((p2.array() - 1.0).abs());
  
  // Calculate numerator terms
  VectorXd nterm1(((p1 - p2).array()*(p1 - p2).array()).matrix());
  VectorXd nterm2((p1.array()*q1.array()).matrix()/double(n1 - 1));
  VectorXd nterm3((p2.array()*q2.array()).matrix()/double(n2 - 1));
  
  // Calculate denominator terms
  VectorXd dterm1((p1.array()*q2.array()).matrix());
  VectorXd dterm2((p2.array()*q1.array()).matrix());
  
  // Compute the average numerator and denominator
  double numer((nterm1 - nterm2 - nterm3).sum()/double(m1));
  double denom((dterm1 + dterm2).sum()/double(m2));
  
  return numer/denom;
}
