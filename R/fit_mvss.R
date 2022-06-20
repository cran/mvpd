#' Fit a Multivariate Subgaussian Distribution
#' 
#' 
#' Estimates the parameters (namely, alpha, shape matrix Q, and location vector) 
#' of the multivariate subgaussian
#' distribution for an input matrix X.
#' 
#' 
#' Using the protocols outlined in Nolan (2013), this function uses \code{libstableR}'s univariate 
#' fit functions for each component.  
#' 
#' 
#' @param x a matrix for which the parameters for a \code{d}-dimensional multivariate 
#' subgaussian distribution will be estimated.  The number of columns will be \code{d}.
#' @return A list with parameters from the column-wise univariate fits and 
#' the multivariate alpha and shape matrix estimates (the \code{univ_deltas} are the \code{mult_deltas}):
#' \itemize{
#'   \item \code{univ_alphas} - the alphas from the column-wise univariate fits
#'   \item \code{univ_betas}  - the betas  from the column-wise univariate fits
#'   \item \code{univ_gammas} - the gammas from the column-wise univariate fits
#'   \item \code{univ_deltas} - the deltas from the column-wise univariate fits
#'   \item \code{mult_alpha}  - the mean(univ_alphas); equivalently the multivariate alpha estimate
#'   \item \code{mult_Q_raw}  - the multivariate shape matrix estimate (before applying \code{nearPD()})
#'   \item \code{mult_Q_posdef}   - the nearest positive definite multivariate shape matrix estimate, \code{nearPD(mult_Q_raw)}
#' }
#' @seealso  \code{Rfast::mvnorm.mle}, \code{alphastable::mfitstab.elliptical}
#' @references
#' Nolan JP (2013), \emph{Multivariate elliptically contoured stable distributions:
#' theory and estimation}. Comput Stat (2013) 28:2067–2089
#' DOI 10.1007/s00180-013-0396-7
#' 
#' @keywords distribution
#' @importFrom utils combn
#' @importFrom libstableR stable_fit_mle2d
#' @importFrom Matrix nearPD
#' @examples
#' \donttest{
#' ## create a 4x4 shape matrix symMat
#' S <- matrix(rnorm(4*4, mean=2, sd=4),4); 
#' symMat <- as.matrix(Matrix::nearPD(0.5 * (S + t(S)))$mat)
#' symMat
#' ## generate 10,000 r.v.'s from 4-dimensional mvss
#' X <- mvpd::rmvss(1e4, alpha=1.5, Q=symMat, delta=c(1,2,3,4))
#' ## use fit_mvss to recover the parameters, compare to symMat
#' fmv <- mvpd::fit_mvss(X)
#' fmv
#' symMat
#' ## then use the fitted parameters to calculate a probability:
#' mvpd::pmvss(lower=rep(0,4),
#'             upper=rep(5,4),
#'             alpha=fmv$mult_alpha,
#'             Q=fmv$mult_Q_posdef,
#'             delta=fmv$univ_deltas,
#'             maxpts.pmvnorm = 25000*10)
#' }
#' 
#' @export
fit_mvss <- function(x){
  
  ## x comes in as an nxd matrix.  
  ## the task is to run univariate fit across all the columns
  component_results <- apply(x, 2, function(W){libstableR::stable_fit_mle2d(W, parametrization = 1L)})
  
  ## now need to compute the shape matrix Q
  ## according to Nolan 2013, we subtract off the delta estimates
  xc <- x - matrix(component_results[4,],nrow=nrow(x),ncol=ncol(x),byrow=T)
  
  ## we know that the diagonal of shape matrix Q is just univariate gammas squared.
  ## to calculate the off-diagonal, we have to add up all pairwise combinations
  ## of components, fit a univariate stable distribution, take that estimated gamma
  ## and subtract off the diagonals and then divide by 2.  What's an efficient way
  ## to do this?
  ## https://stackoverflow.com/a/40139629/2727349
  all.pairs <- utils::combn(1:ncol(xc),2)
  xcp <- xc[, all.pairs[1,]] + xc[, all.pairs[2,]]
  
  if(NCOL(xcp)>1){
    ## now apply the same apply() we did to the original x to this xcp (c-centered p-pairwise sums)
    xcp_results <- apply(xcp, 2, function(W){libstableR::stable_fit_mle2d(W, parametrization = 1L)})
    ## initialize Q_est
    Q_est <- diag(component_results[3,]^2)
    ## do a double assignment to get symmetry, loop over the pairwise combos
    for(i in 1:NCOL(all.pairs)){
      Q_est[all.pairs[1,i], all.pairs[2,i]] <- Q_est[all.pairs[2,i], all.pairs[1,i]] <- (xcp_results[3,i]^2 - Q_est[all.pairs[1,i], all.pairs[1,i]] - Q_est[all.pairs[2,i], all.pairs[2,i]])/2
    }
  }
  
  if(NCOL(xcp)==1){
    ## now apply the same apply() we did to the original x to this xcp (c-centered p-pairwise sums)
    xcp_results <- libstableR::stable_fit_mle2d(xcp, parametrization = 1L)
    ## initialize Q_est
    Q_est <- diag(component_results[3,]^2)
    ## do a double assignment to get symmetry, loop over the pairwise combos
    for(i in 1:NCOL(all.pairs)){
      Q_est[all.pairs[1,i], all.pairs[2,i]] <- Q_est[all.pairs[2,i], all.pairs[1,i]] <- (xcp_results[3  ]^2 - Q_est[all.pairs[1,i], all.pairs[1,i]] - Q_est[all.pairs[2,i], all.pairs[2,i]])/2
    }
  }
  
  return(
    list(
      univ_alphas = component_results[1,],
      univ_betas  = component_results[2,],
      univ_gammas = component_results[3,],
      univ_deltas = component_results[4,],
      mult_alpha  = mean(component_results[1,]),
      mult_Q_raw  = Q_est,
      mult_Q_posdef   = as.matrix(Matrix::nearPD(Q_est)$mat)
    )
  )
  
}


