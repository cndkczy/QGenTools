#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#if !defined(EIGEN_USE_MKL) // Don't use R's Lapack.h if MKL is enabled
#include <R_ext/Lapack.h>
#endif

#include <numeric>
#include <math.h>

// [[Rcpp::plugins(cpp11)]]

using Eigen::BDCSVD;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using Rcpp::_;
using Rcpp::as;
using Rcpp::NumericMatrix;
using Rcpp::wrap;

//' Matrix pseudoinverse
//' 
//' Compute the Moore-Penrose pseudoinverse using the SVD approximation.
//' 
//' To use the same tolerance as \code{MASS::ginv}, set the tolerance parameter 
//' to \code{sqrt(.Machine$double.eps)}.
//' 
//' @param V A numeric matrix
//' @param tol A relative tolerance for detecting zero singular values
//' 
//' @value The Moore-Penrose pseudoinverse for \code{V}.
//' 
//' @export
// [[Rcpp::export]]
MatrixXd pinv(MatrixXd V, double tol = 1.e-6) {
  BDCSVD<MatrixXd> decomp(V, Eigen::ComputeFullU | Eigen::ComputeFullV);
  
  // Set very small singular values to zero
  VectorXd lambdas(decomp.singularValues());
  double tolAdj(lambdas(0)*tol);
  for (size_t i = 0; i < lambdas.size(); ++i)
    if (lambdas(i) > tolAdj)
      lambdas(i) = 1.0/lambdas(i);
    else lambdas(i) = 0.0;
    
    return decomp.matrixV()*lambdas.asDiagonal()*decomp.matrixU().transpose();
}
