#' Monte Carlo Multivariate Subgaussian Stable Distribution
#' 
#' 
#' Computes probabilities of the multivariate subgaussian stable
#' distribution for arbitrary limits, alpha, shape matrices, and
#' location vectors via Monte Carlo (thus the suffix `_mc`).
#' 
#' @param lower the vector of lower limits of length n.
#' @param upper the vector of upper limits of length n.
#' @param alpha default to 1 (Cauchy). Must be 0<alpha<2
#' @param delta location vector.
#' @param Q Shape matrix.  See Nolan (2013).
#' @param which.stable defaults to "libstableR", other option is "stabledist".  Indicates which package 
#' should provide the univariate stable distribution in this production distribution form of a univariate
#' stable and multivariate normal.
#' @param n number of random vectors to be drawn for Monte Carlo calculation.
#' @return a number between 0 and 1, the estimated probability via Monte Carlo
#' @references
#' Nolan JP (2013), \emph{Multivariate elliptically contoured stable distributions:
#' theory and estimation}. Comput Stat (2013) 28:2067–2089
#' DOI 10.1007/s00180-013-0396-7
#' @examples
#' 
#' ## print("mvpd (d=2, alpha=1.71):")
#' U <- c(1,1)
#' L <- -U
#' Q <- matrix(c(10,7.5,7.5,10),2)
#' mvpd::pmvss_mc(L, U, alpha=1.71, Q=Q, n=1e3)
#' mvpd::pmvss   (L, U, alpha=1.71, Q=Q)
#'
#' ## more accuracy = longer runtime
#' mvpd::pmvss_mc(L, U, alpha=1.71, Q=Q, n=1e4)
#' 
#' U <- c(1,1,1)
#' L <- -U
#' Q <- matrix(c(10,7.5,7.5,7.5,10,7.5,7.5,7.5,10),3)
#' ## print("mvpd: (d=3, alpha=1.71):")
#' mvpd::pmvss_mc(L, U, alpha=1.71, Q=Q, n=1e3)
#' 
#' 
#' @export
pmvss_mc <- function(lower=rep(-Inf,d), upper=rep(Inf,d), alpha=1, Q = NULL, delta=rep(0,d), 
                          
                          which.stable=c("libstableR", "stabledist")[1],
                          
                          n=NULL
                          
                          ){
  d <- nrow(Q)
  rmv <- mvpd::rmvss(n, alpha, Q, delta, which.stable)
  mc.result <- sum(rowSums(rmv >= lower & rmv <= upper) == nrow(Q))/n 
  mc.result
}


