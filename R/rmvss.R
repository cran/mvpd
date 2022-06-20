#' Multivariate Subgaussian Stable Random Variates
#'
#'
#' Computes random vectors of the multivariate subgaussian stable
#' distribution for arbitrary alpha, shape matrices, and 
#' location vectors. See Nolan (2013).
#'
#'
#' @param n number of observations
#' @param alpha default to 1 (Cauchy). Must be 0<\code{alpha}<2
#' @param delta location vector.
#' @param Q Shape matrix.  See Nolan (2013).
#' @param which.stable defaults to \code{"libstableR"}, other option is \code{"stabledist"}.  Indicates which package
#' should provide the univariate stable distribution in this production distribution form of a univariate
#' stable and multivariate normal.
#' @return Returns the \code{n} by \code{d} matrix containing multivariate subgaussian stable 
#' random variates where \code{d=nrow(Q)}.
#' @references
#' Nolan JP (2013), \emph{Multivariate elliptically contoured stable distributions:
#' theory and estimation}. Comput Stat (2013) 28:2067–2089
#' DOI 10.1007/s00180-013-0396-7
#' @keywords distribution
#' @importFrom matrixStats rowProds
#' @importFrom stabledist rstable
#' @importFrom libstableR stable_rnd
#' @importFrom mvtnorm rmvnorm GenzBretz
#' @examples
#' ## generate 10 random variates of a bivariate mvss
#' rmvss(n=10, alpha=1.71, Q=matrix(c(10,7.5,7.5,10),2))
#' 
#' ## generate 10 random variates of a trivariate mvss
#' Q <- matrix(c(10,7.5,7.5,7.5,10,7.5,7.5,7.5,10),3)
#' rmvss(n=10, alpha=1.71, Q=Q)
#' 
#'
#' @export
rmvss <- function(n, alpha=1, Q = NULL, delta=rep(0,d),
                  which.stable=c("libstableR", "stabledist")[1]
){
  d <- nrow(Q)

  A <- switch(which.stable,
              "libstableR" =

                libstableR::stable_rnd(n,
                                       pars=c(
                                         alpha = alpha/2,
                                         beta = 1,
                                         sigma = 2*cos(pi*alpha/4)^(2/alpha),
                                         mu = 0),
                                       parametrization = 1L)
              ,
              "stabledist"=

                stabledist::rstable(n,
                                    alpha = alpha/2,
                                    beta = 1,
                                    gamma = 2*cos(pi*alpha/4)^(2/alpha),
                                    delta = 0,
                                    pm = 1)
  )


  G <-   mvtnorm::rmvnorm(n, sigma=Q)


  sqrt(A) * G + matrix(delta,nrow=n,ncol=length(delta),byrow=T)


}
