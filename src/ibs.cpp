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

//' Identity-by-state between genotypes.
//' 
//' Calculates the proportion of shared alleles between pairs of genotypes.
//' 
//' @param X A numeric matrix of allelic dosages where rows represent different
//'   genotypes and columns represent different markers. If \code{X} has rownames, 
//'   these will be preserved in the output.
//' 
//' @details Each element of \code{X} should be in the interval \eqn{[0,2]}. 
//'   The matrix is pre-multiplied by one-half to correspond to equation 1 of 
//'   Bustos-Korts \emph{et al.} (2016).
//' 
//' @return A numeric matrix containing the proportion of shared alleles between
//'   each pair of genotypes.
//' 
//' @references Bustos-Korts \emph{et al.} (2016) Improvement of predictive 
//'   ability by uniform coverage of the target genetic space. \emph{G3}, 
//'   \strong{6}: 3733-47.
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix ibs(NumericMatrix X) {
  MatrixXd G(as<MatrixXd>(X));
  G *= 0.5;
  int r(G.rows()), c(G.cols());
  
  MatrixXd GGt(MatrixXd(r, r).setZero().selfadjointView<Lower>().rankUpdate(G));
  MatrixXd G2(MatrixXd::Ones(r, r) - G);
  MatrixXd G2G2t(MatrixXd(r, r).setZero().selfadjointView<Lower>().rankUpdate(G2));
  
  // Calculate identity-by-state
  NumericMatrix IBS(wrap((1/double(c))*(GGt + G2G2t)));
  rownames(IBS) = rownames(X);
  colnames(IBS) = colnames(X);
  
  return IBS;
}
