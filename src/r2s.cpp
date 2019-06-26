#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#if !defined(EIGEN_USE_MKL) // Don't use R's Lapack.h if MKL is enabled
#include <R_ext/Lapack.h>
#endif

#include <numeric>
#include <math.h>

// [[Rcpp::depends(FastMath)]]
#include <FastMath.h>

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

//' @export
// [[Rcpp::export]]
DataFrame r2s(NumericMatrix X, NumericMatrix st) {
  MatrixXd G(as<MatrixXd>(X));
  MatrixXd S(as<MatrixXd>(st));
  
  // Concatenate the genotypes and structure matrix
  MatrixXd GS(G.rows(), G.cols() + S.cols());
  GS << G, S;
  
  // Calculate unadjusted covariances and variances
  MatrixXd centered(GS.rowwise() - GS.colwise().mean());
  MatrixXd sigmaXS((centered.adjoint()*centered)/double(GS.rows() - 1));
  
  // Separate blocks for easier calculations
  MatrixXd XX(sigmaXS.topLeftCorner(G.cols(), G.cols()));
  MatrixXd SX(sigmaXS.bottomLeftCorner(S.cols(), G.cols()));
  MatrixXd XS(sigmaXS.topRightCorner(G.cols(), S.cols()));
  MatrixXd Si(as<MatrixXd>(FastMath::pinv(wrap(sigmaXS.bottomRightCorner(S.cols(), S.cols())))));
  
  // Calculate the adjusted covariances
  MatrixXd sigmaV(XX - XS*Si*SX);
  
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
