#' Linkage disequilibrium decay distance.
#' 
#' Estimates the linkage disequilibrium (LD) decay distance by fitting a 
#' monotonically decreasing spline to estimates of pairwise LD.
#' 
#' @param ld A data frame with two columns: \code{Bin} contains the distance 
#'   between SNP pairs (or bin widths) and \code{R2} contains the estimated LD.
#' @param threshold The strength at which LD is considered decayed (default: 0.3).
#' @param k The number of knots to use for spline fitting (default: 30)
#' 
#' @details This function estimates the LD decay distance by fitting a monotonically 
#'   decreasing spline to empirical estimates of LD. The units of \code{Bin} may 
#'   vary, depending on application. For example, local LD decay around an 
#'   individual SNP can be calculated where each row of \code{ld} is the measure 
#'   of LD between that SNP and a second SNP with the base-pair distance between 
#'   them. Alternatively, global LD decay could be calculated where LD measures 
#'   are aggregated into bins of 20,000 bp, for example.
#' 
#' @return The distance in units of the \code{Bin} column at which LD is estimated 
#'   to reach \code{threshold}.
#'   
#' @references Bustos-Korts \emph{et al.} (2019) Exome sequences and multi-environment 
#'   field trials elucidate the genetic basis of adaptation in barly. \emph{The Plant Journal}, 
#'   doi: 10.1111/tpj.14414.
#' @references Based on code from <https://gist.github.com/jnpaulson/c47f9bd3246f1121ad3a911e5f707355>.
#'
#' @export
est_distance <- function(ld, threshold = 0.3, k = 30) {
  # Fit an initial smoothing spline
  init_gam <- mgcv::gam(R2 ~ s(Bin, k = k, bs = "cr"), data = ld)
  
  # Estimate the penalized smooth based on the initial fit
  sm <- mgcv::smoothCon(s(Bin, k = k, bs = "cr"), ld, knots = NULL)[[1]]
  mc <- mgcv::mono.con(sm$xp, up = FALSE)
  M <- list(X = sm$X, y = ld$R2, C = matrix(0, 0, 0), Ain = mc$A, 
            bin = mc$b, sp = init_gam$sp, p = -sm$xp, S = sm$S, 
            w = ld$R2*0 + 1, off = 0)
  p <- mgcv::pcls(M)
  
  # If the highest estimated LD is less than the threshold, return the smallest
  # SNP-SNP distance as the lower bound of the decay estimate.
  if (all(is.nan(p))) {
    min(ld$Bin)
  } else {
    x <- seq(min(ld$Bin), max(ld$Bin), 1e3)
    y <- abs(drop(mgcv::Predict.matrix(sm, data.frame(Bin = x)) %*% p) - threshold)
    x[which.min(y)]
  }
}
