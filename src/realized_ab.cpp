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

//' Realized genetic relationships between genotypes.
//' 
//' Calculates pairwise genetic relationships based on marker scores.
//' 
//' @param X A numeric matrix of allelic dosages. The rows should correspond to 
//'   genotypes and the columns should correspond to markers. If \code{X} has 
//'   rownames, these will be preserved in the output.
//' 
//' @details Each element of \code{X} should be in the interval \eqn{[0,2]}.
//' 
//' @return The matrix of realized genetic relationships (variances and covariances) 
//'   between each pair of genotypes.
//' 
//' @references Astle and Balding. (2009) Population structure and cryptic 
//'   relatedness in genetic association studies. \emph{Statistical Science}, 
//'   \strong{24}\emph{(4)}: 451-71.
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix realized_ab(NumericMatrix X) {
  MatrixXd G(as<MatrixXd>(X));
  int r(G.rows()), c(G.cols());
  
  // Calculate allele frequencies
  VectorXd freqs(G.colwise().sum()/(2*double(r)));
  double sd(0);
  for (size_t i = 0; i < c; ++i)
    sd = sd + freqs(i)*(1 - freqs(i));
  sd = 2*(sd/double(c));
  
  // Center the marker matrix
  for (size_t i = 0; i < c; ++i)
    G.col(i) = (G.col(i).array() - 2*freqs(i));
  
  // Form the cross product
  MatrixXd GGt(G*G.transpose());
  
  // Realized AB matrix
  NumericMatrix AB(wrap(GGt/sd/double(c)));
  rownames(AB) = rownames(X);
  colnames(AB) = rownames(X);
  
  return AB;
}
