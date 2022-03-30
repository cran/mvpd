#' Multivariate Subgaussian Stable Distribution
#' 
#' 
#' Computes the probabilities for the multivariate subgaussian stable
#' distribution for arbitrary limits, alpha, shape matrices, and 
#' location vectors.
#' See Nolan (2013).
#' 
#' @param lower the vector of lower limits of length n.
#' @param upper the vector of upper limits of length n.
#' @param alpha default to 1 (Cauchy). Must be 0<alpha<2
#' @param Q Shape matrix.  See Nolan (2013).
#' @param delta location vector.
#' @param maxpts.pmvnorm Defaults to 25000.  Passed to the F_G = pmvnorm() in the integrand of the outermost integral.
#' @param abseps.pmvnorm Defaults to 1e-3. Passed to the F_G = pmvnorm() in the integrand of the outermost integral.
#' @param outermost.int select which integration function to use for outermost
#' integral.  Default is "stats::integrate" and one can specify the following options 
#' with the `.si` suffix.  For diagonal Q, one can also specify "cubature::adaptIntegrate"
#' and use the `.ai` suffix options below (currently there is a bug for non-diagonal Q).
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
#' @param which.stable defaults to "libstableR", other option is "stabledist".  Indicates which package 
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
#' Nolan JP (2013), \emph{Multivariate elliptically contoured stable distributions:
#' theory and estimation}. Comput Stat (2013) 28:2067–2089
#' DOI 10.1007/s00180-013-0396-7
#' 
#' @keywords distribution
#' @importFrom stabledist dstable
#' @importFrom libstableR stable_pdf
#' @importFrom mvtnorm pmvnorm GenzBretz
#' @importFrom stats integrate
#' @examples
#' 
#' ## bivariate
#' U <- c(1,1)
#' L <- -U
#' Q <- matrix(c(10,7.5,7.5,10),2)
#' mvpd::pmvss(L, U, alpha=1.71, Q=Q)
#' 
#' ## trivariate
#' U <- c(1,1,1)
#' L <- -U
#' Q <- matrix(c(10,7.5,7.5,7.5,10,7.5,7.5,7.5,10),3)
#' mvpd::pmvss(L, U, alpha=1.71, Q=Q)
#' 
#' ## How `delta` works: same as centering
#' U <- c(1,1,1)
#' L <- -U
#' Q <- matrix(c(10,7.5,7.5,7.5,10,7.5,7.5,7.5,10),3)
#' D <- c(0.75, 0.65, -0.35)
#' mvpd::pmvss(L-D, U-D, alpha=1.71, Q=Q)
#' mvpd::pmvss(L  , U  , alpha=1.71, Q=Q, delta=D)
#' 
#' 
#' 
#' @export
pmvss <-
  function(lower=rep(-Inf,d), upper=rep(Inf,d), alpha=1, Q = NULL, delta=rep(0,d), 
           
           
           maxpts.pmvnorm = 25000,
           abseps.pmvnorm = 1e-3,
           
           outermost.int = c("stats::integrate","cubature::adaptIntegrate")[1],
           
           subdivisions.si = 100L,
           rel.tol.si = .Machine$double.eps^0.25, 
           abs.tol.si = rel.tol.si,
           stop.on.error.si = TRUE,
           
           tol.ai = 1e-05,
           fDim.ai = 1,
           maxEval.ai = 0,
           absError.ai=0,
           doChecking.ai=FALSE,
           
           which.stable=c("libstableR", "stabledist")[1]
           
  ){
    
    
    f_A <- switch(which.stable,
                  "libstableR" =
                    function(x, alpha) 
                      libstableR::stable_pdf(x,
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
                  "libstableR" = 
                    function(x, alpha)
                      2*x*
                    libstableR::stable_pdf(x^2,
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
    
    
    d <- length(lower)
    lower <- lower - delta
    upper <- upper - delta
    integrand <- function(b,llp, ulp, Qmat=Q, alp=alpha){ 
      F_G <- mvtnorm::pmvnorm(lower=llp/b, 
                              upper=ulp/b, 
                              sigma=Qmat,
                              algorithm=mvtnorm::GenzBretz(maxpts = maxpts.pmvnorm, 
                                                           abseps=abseps.pmvnorm, 
                                                           releps = 0))
      # print(Sys.time())
      if( attr(F_G, "error") > abseps.pmvnorm){ print(paste0("WARNING: EARLY/ATYPICAL COMPLETION of numerical F_G:")) 
                                                print(paste0("abseps.pmvnorm for F_G is ", abseps.pmvnorm, ", ")) 
                                                print(paste0("and the actual error returned is more than that at ", attr(F_G, "error"))) 
                                                print(paste0("Consider increasing `maxpts.pmvnorm` AND/OR enlarging `abseps.pmvnorm`.")  )
                                                }
      
      f_B(b, alp) * F_G[1]
    }
    # integrand <- function(b,llp, ulp, Qmat=Q, alp=alpha) 
    #   f_B(b, alp) * mvtnorm::pmvnorm(lower=llp/b, upper=ulp/b, sigma=Qmat,
    #                                  algorithm=mvtnorm::GenzBretz(maxpts = maxpts.pmvnorm, abseps=abseps.pmvnorm, releps = 0))
    # 
    # integrand2 <- function(b,llp, ulp, Qmat=Q, alp=alpha) 
    #   f_B(b, alp) * mvtnorm::pmvnorm(llp/b, ulp/b)
    
    
    switch(outermost.int,
           "cubature::adaptIntegrate"=
             adaptIntegrate_inf_limPD(integrand,
                                      lowerLimit=0,
                                      upperLimit=Inf,
                                      llp=lower,
                                      ulp=upper,
                                      tol.ai = tol.ai,
                                      fDim.ai = fDim.ai,
                                      maxEval.ai = maxEval.ai,
                                      absError.ai=absError.ai,
                                      doChecking.ai=doChecking.ai
             )
           
           ,
           "stats::integrate" =     
             integrate(Vectorize(integrand,vectorize.args = "b"),
                       lower=0,
                       upper=Inf,
                       llp=lower,
                       ulp=upper,
                       subdivisions = subdivisions.si,
                       rel.tol=rel.tol.si,
                       abs.tol=abs.tol.si,
                       stop.on.error=stop.on.error.si)
    )
    
    
    
    
    
  }