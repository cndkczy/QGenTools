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

using Rcpp::_;
using Rcpp::as;
using Rcpp::CharacterVector;
using Rcpp::DataFrame;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::wrap;

//' Correlation between loci.
//' 
//' Calculates the correlation (\eqn{r^2}) between marker loci within a population.
//' 
//' @param X An \eqn{n} x \eqn{m} numeric matrix of allelic dosages. The matrix must 
//'   have column names for labelling correlations in the output.
//' 
//' @details Elements of \code{X} should be in the interval \eqn{[0,2]}. This 
//'   computes the traditional \eqn{r^2} measure of linkage disequilibrium 
//'   between pairs of marker loci.
//' 
//' @return A dataframe with \eqn{m(m-1)/2} rows containing pairwise correlations
//'   between all possible pairs of markers in \code{X}.
//' 
//' @export
// [[Rcpp::export]]
DataFrame r2(NumericMatrix X) {
  MatrixXd G(as<MatrixXd>(X));
  
  // Calculate unadjusted covariances and variances
  MatrixXd centered(G.rowwise() - G.colwise().mean());
  MatrixXd sigmaV((centered.adjoint()*centered)/double(G.rows() - 1));
  
  // Initialize results vectors
  CharacterVector loci = colnames(X);
  int nloci(G.cols()); // Number of loci
  CharacterVector locus1(nloci*(nloci - 1)/2);
  CharacterVector locus2(nloci*(nloci - 1)/2);
  NumericVector R2(nloci*(nloci - 1)/2);
  
  // Compute corrected LD
  int counter(0);
  for (size_t i = 0; i < sigmaV.rows(); ++i) {
    for (size_t j = i + 1; j < sigmaV.cols(); ++j) {
      // Record the loci
      locus1(counter) = loci(i);
      locus2(counter) = loci(j);
      
      // Compute the correlation if the variances are non-zero
      if (sigmaV(i, i) < 1.0e-7 | sigmaV(j, j) < 1.0e-7)
        R2(counter) = 0.0;
      else R2(counter) = std::pow(sigmaV(i, j), 2)/(sigmaV(i, i)*sigmaV(j, j));
      
      // Advance the counter
      counter += 1;
    }
  }
  
  return DataFrame::create(_["Locus1"] = locus1, 
                           _["Locus2"] = locus2, 
                           _["R2"] = R2);
}
