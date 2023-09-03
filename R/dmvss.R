#' Multivariate Subgaussian Stable Density
#' 
#' 
#' Computes the the density function of the multivariate subgaussian stable
#' distribution for arbitrary alpha, shape matrices, and location vectors.
#' See Nolan (2013).
#' 
#' @param x vector of quantiles.
#' @param alpha default to 1 (Cauchy). Must be 0<alpha<2
#' @param delta location vector
#' @param Q Shape matrix.  See Nolan (2013).
#' @param outermost.int select which integration function to use for outermost
#' integral.  Default is "stats::integrate" and one can specify the following options 
#' with the `.si` suffix.  For diagonal Q, one can also specify "cubature::adaptIntegrate"
#' and use the `.ai` suffix options below (currently there is a bug for non-diagnoal Q).
#' @param spherical default is FALSE.  If true, use the spherical transformation.  
#' Results will be identical to spherical = FALSE but may be faster.
#' @param subdivisions.si the maximum number of subintervals.
#' The suffix \code{.si} indicates a \code{stats::integrate()} 
#' option for the outermost semi-infinite integral in the product distribution formulation.
#' @param rel.tol.si relative accuracy requested.  
#' The suffix \code{.si} indicates a \code{stats::integrate()} 
#' option for the outermost semi-infinite integral in the product distribution formulation.
#' @param abs.tol.si absolute accuracy requested. The suffix \code{.si} indicates a \code{stats::integrate()} 
#' option for the outermost semi-infinite integral in the product distribution formulation.
#' @param stop.on.error.si logical. If true (the default) an error stops the function. 
#' If false some errors will give a result with a warning in the message component.
#' The suffix \code{.si} indicates a \code{stats::integrate()} 
#' option for the outermost semi-infinite integral in the product distribution formulation.
#' @param tol.ai The maximum tolerance, default 1e-5.
#' The suffix \code{.ai} indicates a \code{cubature::adaptIntegrate} type 
#' option for the outermost semi-infinite integral in the product distribution formulation.
#' @param fDim.ai The dimension of the integrand, default 1, bears no
#' relation to the dimension of the hypercube
#' The suffix \code{.ai} indicates a \code{cubature::adaptIntegrate} type 
#' option for the outermost semi-infinite integral in the product distribution formulation.
#' @param maxEval.ai The maximum number of function evaluations needed,
#' default 0 implying no limit
#' The suffix \code{.ai} indicates a \code{cubature::adaptIntegrate} type 
#' option for the outermost semi-infinite integral in the product distribution formulation.
#' @param absError.ai The maximum absolute error tolerated
#' The suffix \code{.ai} indicates a \code{cubature::adaptIntegrate} type 
#' option for the outermost semi-infinite integral in the product distribution formulation.
#' @param doChecking.ai A flag to be thorough checking inputs to
#' C routines. A FALSE value results in approximately 9 percent speed
#' gain in our experiments. Your mileage will of course vary. Default
#' value is FALSE.
#' The suffix \code{.ai} indicates a \code{cubature::adaptIntegrate} type 
#' option for the outermost semi-infinite integral in the product distribution formulation.
#' @param which.stable defaults to "libstable4u", other option is "stabledist".  Indicates which package 
#' should provide the univariate stable distribution in this production distribution form of a univariate
#' stable and multivariate normal.
#' @return The object returned depends on what is selected for \code{outermost.int}.  In the case of the default, 
#' \code{stats::integrate}, the value is a list of class "integrate" with components:
#' \itemize{
#' \item{\code{value}}{ the final estimate of the integral.}
#' \item{\code{abs.error}}{ estimate of the modulus of the absolute error.}
#' \item{\code{subdivisions}}{ the number of subintervals produced in the subdivision process.}
#' \item{\code{message}}{ "OK" or a character string giving the error message.}
#' \item{\code{call}}{ the matched call.}
#' }
#' Note: The reported \code{abs.error} is likely an under-estimate as \code{integrate}
#' assumes the integrand was without error, 
#' which is not the case in this application.  
#' @references
#' 
#' Nolan, John P. "Multivariate elliptically contoured stable distributions: theory and estimation." Computational Statistics 28.5 (2013): 2067-2089.
#' 
#' @keywords distribution
#' @importFrom cubature adaptIntegrate
#' @importFrom matrixStats rowProds
#' @importFrom stabledist dstable
#' @importFrom libstable4u stable_pdf
#' @importFrom mvtnorm pmvnorm GenzBretz
#' @importFrom stats integrate dnorm
#' @examples
#' 
#' ## print("mvsubgaussPD (d=2, alpha=1.71):")
#' Q <- matrix(c(10,7.5,7.5,10),2)
#' mvpd::dmvss(x=c(0,1), alpha=1.71, Q=Q)
#' 
#' ## more accuracy = longer runtime
#' mvpd::dmvss(x=c(0,1),alpha=1.71, Q=Q, abs.tol=1e-8)
#' 
#' Q <- matrix(c(10,7.5,7.5,7.5,10,7.5,7.5,7.5,10),3)
#' ## print("mvsubgausPD (d=3, alpha=1.71):")
#' mvpd::dmvss(x=c(0,1,2), alpha=1.71, Q=Q)
#' mvpd::dmvss(x=c(0,1,2), alpha=1.71, Q=Q, spherical=TRUE)
#' 
#' ## How `delta` works: same as centering
#' X <- c(1,1,1)
#' Q <- matrix(c(10,7.5,7.5,7.5,10,7.5,7.5,7.5,10),3)
#' D <- c(0.75, 0.65, -0.35)
#' mvpd::dmvss(X-D, alpha=1.71, Q=Q)
#' mvpd::dmvss(X  , alpha=1.71, Q=Q, delta=D)
#' 
#' 
#' @export
dmvss <- function(x, alpha=1, Q = NULL, delta=rep(0,d), 
                          
                          outermost.int = c("stats::integrate","cubature::adaptIntegrate")[1],
                          
                          spherical = FALSE,
                          
                          subdivisions.si = 100L,
                          rel.tol.si = .Machine$double.eps^0.25, 
                          abs.tol.si = rel.tol.si,
                          stop.on.error.si = TRUE,
                          
                          tol.ai = 1e-05,
                          fDim.ai = 1,
                          maxEval.ai = 0,
                          absError.ai=0,
                          doChecking.ai=FALSE,
                          
                          which.stable=c("libstable4u", "stabledist")[1]
                          
                          ){
  
  f_A <- switch(which.stable,
                "libstable4u" =
                  function(x, alpha) 
                    libstable4u::stable_pdf(x,
                                           pars=c(
                                             alpha = alpha/2,
                                             beta = 1,
                                             sigma = 2*cos(pi*alpha/4)^(2/alpha),
                                             mu = 0),
                                           parametrization = 1L)
                ,
                "stabledist"=
                  function(x, alpha) 
                    stabledist::dstable(x,
                                        alpha = alpha/2,
                                        beta = 1,
                                        gamma = 2*cos(pi*alpha/4)^(2/alpha),
                                        delta = 0,
                                        pm = 1)
  )
  
  
  
  
  f_B <- switch(which.stable,
                "libstable4u" = 
                  function(x, alpha)
                    2*x*
                  libstable4u::stable_pdf(x^2,
                                         pars=c(
                                           alpha = alpha/2,
                                           beta = 1,
                                           sigma = 2*cos(pi*alpha/4)^(2/alpha),
                                           mu = 0),
                                         parametrization = 1L)
                ,
                "stabledist" = 
                  function(x, alpha)
                    2*x*
                  stabledist::dstable(x^2,
                                      alpha = alpha/2,
                                      beta = 1,
                                      gamma = 2*cos(pi*alpha/4)^(2/alpha),
                                      delta = 0,
                                      pm = 1)
  )  
    d <- length(x)
    x <- x - delta
    ifelse(!spherical,
           dens <- function(b,x,Qmat=Q,alp=alpha) f_B(b, alp) * 1/abs(b)^d * mvtnorm::dmvnorm(x/b, sigma=Qmat),
           dens <- function(b,x,Qmat=Q,alp=alpha){
             Linv <- solve(t(chol(Qmat)))
             detLinv <- det(Linv)
             f_B(b, alpha) * 1/abs(b)^d *detLinv *prod(dnorm(t(Linv %*% matrix(x/b,ncol=1))))
           }
    )
    
  switch(outermost.int,
         "cubature::adaptIntegrate"=
           adaptIntegrate_inf_limPD(dens,
                                    lowerLimit=0,
                                    upperLimit=Inf,
                                    x=x,
                                    Q=Q,
                                    tol.ai = tol.ai,
                                    fDim.ai = fDim.ai,
                                    maxEval.ai = maxEval.ai,
                                    absError.ai=absError.ai,
                                    doChecking.ai=doChecking.ai
                                    )
         
         ,
         "stats::integrate" =     
           integrate(Vectorize(dens,vectorize.args = "b"),
                     lower=0,
                     upper=Inf,
                     x=x,
                     Q=Q,
                     subdivisions = subdivisions.si,
                     rel.tol=rel.tol.si,
                     abs.tol=abs.tol.si,
                     stop.on.error=stop.on.error.si)
  )
  
  
  
  
  
}

# use ones in pmvsubgaussPD.R (now they are inside each [dp]mvsubgaussPD)
# f_A <- function(x, alpha) 
#   stabledist::dstable(x,
#                       alpha = alpha/2,
#                       beta = 1,
#                       gamma = 2*cos(pi*alpha/4)^(2/alpha),
#                       delta = 0,
#                       pm = 1)
# f_B <- function(x, alpha)
#   2*x*
#   stabledist::dstable(x^2,
#                       alpha = alpha/2,
#                       beta = 1,
#                       gamma = 2*cos(pi*alpha/4)^(2/alpha),
#                       delta = 0,
#                       pm = 1)
