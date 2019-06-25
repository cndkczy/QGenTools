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
using Rcpp::CharacterVector;
using Rcpp::DataFrame;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::wrap;

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

//' @export
// [[Rcpp::export]]
DataFrame r2vs(NumericMatrix X, NumericMatrix V, NumericMatrix st) {
  MatrixXd G(as<MatrixXd>(X));
  MatrixXd Vi(pinv(as<MatrixXd>(V)));
  MatrixXd S(as<MatrixXd>(st));
  
  // Concatenate the genotypes and structure matrix
  MatrixXd GS(G.rows(), G.cols() + S.cols());
  GS << G, S;
  
  // Calculate the relationship adjustment factor
  MatrixXd ones(MatrixXd::Ones(GS.rows(), 1));
  MatrixXd num(ones*ones.transpose()*Vi);
  double denom((ones.transpose()*Vi*ones)(0, 0));
  MatrixXd factor((1/denom)*num);
  
  // Adjust the observed genotypes and structure for relationships
  MatrixXd adjG(GS - factor*GS);
  
  // Calculate corrected covariances and variances
  MatrixXd sigmaXS(adjG.transpose()*Vi*adjG);
  
  // Separate blocks for easier calculations
  MatrixXd XX(sigmaXS.topLeftCorner(G.cols(), G.cols()));
  MatrixXd SX(sigmaXS.bottomLeftCorner(S.cols(), G.cols()));
  MatrixXd XS(sigmaXS.topRightCorner(G.cols(), S.cols()));
  MatrixXd Si(pinv(sigmaXS.bottomRightCorner(S.cols(), S.cols())));
  
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
