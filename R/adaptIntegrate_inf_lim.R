#' Adaptive multivariate integration over hypercubes (admitting infinite limits)
#'
#' The function performs adaptive multidimensional integration
#' (cubature) of (possibly) vector-valued integrands over hypercubes.
#' It is a wrapper for cubature:::adaptIntegrate, transforming (-)Inf
#' appropriately as described in cubature's help page (http://ab-initio.mit.edu/wiki/index.php/Cubature#Infinite_intervals).
#'
#' @param f The function (integrand) to be integrated
#' @param lowerLimit The lower limit of integration, a vector for
#' hypercubes
#' @param upperLimit The upper limit of integration, a vector for
#' hypercubes
#' @param ... All other arguments passed to the function f
#' @param tol.ai The maximum tolerance, default 1e-5.
#' @param fDim.ai The dimension of the integrand, default 1, bears no
#' relation to the dimension of the hypercube
#' @param maxEval.ai The maximum number of function evaluations needed,
#' default 0 implying no limit
#' @param absError.ai The maximum absolute error tolerated
#' @param doChecking.ai A flag to be thorough checking inputs to
#' C routines. A FALSE value results in approximately 9 percent speed
#' gain in our experiments. Your mileage will of course vary. Default
#' value is FALSE.
#'
#' @keywords multivariate numerical integration
#' @export 
#' 
#' 
#' @examples
#' ## integrate Cauchy Density from -Inf to Inf
#' adaptIntegrate_inf_limPD(function(x) 1/pi * 1/(1+x^2), -Inf, Inf)
#' adaptIntegrate_inf_limPD(function(x, scale) 1/(pi*scale) * 1/(1+(x/scale)^2), -Inf, Inf, scale=4)
#' ## integrate Cauchy Density from -Inf to -3
#' adaptIntegrate_inf_limPD(function(x) 1/pi * 1/(1+x^2), -Inf, -3)$int
#' stats::pcauchy(-3)
#' adaptIntegrate_inf_limPD(function(x, scale) 1/(pi*scale) * 1/(1+(x/scale)^2), -Inf, -3, scale=4)$int
#' stats::pcauchy(-3, scale=4)
#' 
adaptIntegrate_inf_limPD <- function(f, lowerLimit, upperLimit, ..., 
                                   tol.ai = 1e-05, fDim.ai = 1,
                                   maxEval.ai = 0, absError.ai=0, 
                                   doChecking.ai=FALSE
                                   )
{
  ## This function is merely a wrapper that will do an
  ## automated change of variables substitution
  ## to accommodate (-)Inf limits.  See
  ## http://ab-initio.mit.edu/wiki/index.php/Cubature#Infinite_intervals
  ## for more info.

  ## Steps:
  ## 1.  Identify the three classes:
  ##   SL:  semi-infinite bounds with lower -Inf and finite upper bound,
  ##   SU:  semi-infinite bounds with upper  Inf and finite lower bound,
  ##   IB:       infinite bounds with lower -Inf and upper Inf
  ## 2.  Transform the bounds
  ## 3.  Assign the jacobian form
  ## 4.  Substitute new variable parameterization in new function
  ## 5.  Let cubature:::adaptIntegrate() operate on this new (but equivalent) formulation



  ## 0.  assign inputs
  f_orig <- f
  lowerLimit_orig <- lowerLimit
  upperLimit_orig <- upperLimit

  ## get logical info on bounds
  lowerLimit_logical <- is.infinite(lowerLimit_orig)
  upperLimit_logical <- is.infinite(upperLimit_orig)
  SL <-  lowerLimit_logical & !upperLimit_logical
  SU <- !lowerLimit_logical &  upperLimit_logical
  IB <-  lowerLimit_logical &  upperLimit_logical

  jacobian <- function(t, semilower=SL, semiupper=SU, infboth=IB){

    dimt <- dim(t)

    if(!is.null(dimt)){
      semilower_mat <- matrix(rep(SL,dimt[1]), nrow=dimt[1], byrow=T)
      semiupper_mat <- matrix(rep(SU,dimt[1]), nrow=dimt[1], byrow=T)
      infboth_mat <- matrix(rep(IB,dimt[1]), nrow=dimt[1], byrow=T)
    }else{
      semilower_mat <- semilower
      semiupper_mat <- semiupper
      infboth_mat   <- infboth
    }



    ( - 1     / (1 - t  )^2 )^(as.numeric(semilower_mat)) *
      (   1     / (1 - t  )^2 )^(as.numeric(semiupper_mat)) *
      ( (1+t^2) / (1 - t^2)^2 )^(as.numeric(infboth_mat))

  }

  eval_pts_new <- function(t, semilower=SL, semiupper=SU, infboth=IB,
                           ll = lowerLimit_orig, ul = upperLimit_orig){

    dimt <- dim(t)

    if(!is.null(dimt)){
      semilower_mat <- matrix(rep(SL,dimt[1]), nrow=dimt[1], byrow=T)
      semiupper_mat <- matrix(rep(SU,dimt[1]), nrow=dimt[1], byrow=T)
      infboth_mat <- matrix(rep(IB,dimt[1]), nrow=dimt[1], byrow=T)
      ll_mat      <- matrix(rep(ll,dimt[1]), nrow=dimt[1], byrow=T)
      ul_mat      <- matrix(rep(ul,dimt[1]), nrow=dimt[1], byrow=T)
    }else{
      semilower_mat <- semilower
      semiupper_mat <- semiupper
      infboth_mat   <- infboth
      ll_mat      <- ll
      ul_mat      <- ul
    }



    (ul_mat +  - t     / (1 - t  ) )^(as.numeric(semilower_mat)) *
      (ll_mat +    t     / (1 - t  ) )^(as.numeric(semiupper_mat)) *
      (            t     / (1 - t^2) )^(as.numeric(infboth_mat))   *
      t                  ^(as.numeric(!semilower_mat&!semiupper_mat&!infboth_mat))



  }


  f_new <- function(t, f=f_orig, ll = lowerLimit_orig, ul = upperLimit_orig, ...){

    dimt <- dim(t)

    if(!is.null(dimt)){
      f(eval_pts_new(t), ...) * rowProds(jacobian(t))
    }else{
      f(eval_pts_new(t), ...) * prod(jacobian(t))
    }

  }





  new_lowerlimits <- function(ll = lowerLimit_orig, semilower=SL, semiupper=SU, infboth=IB){

    ll[semilower] <-  1
    ll[semiupper] <-  0
    ll[infboth  ] <- -1

    ll

  }


  new_upperlimits <- function(ul = upperLimit_orig, semilower=SL, semiupper=SU, infboth=IB){

    ul[semilower] <-  0
    ul[semiupper] <-  1
    ul[infboth  ] <-  1

    ul

  }


  cubature::adaptIntegrate(f=f_new,
                 lowerLimit=new_lowerlimits(lowerLimit_orig),
                 upperLimit=new_upperlimits(upperLimit_orig),
                 ...,
                 tol = tol.ai, fDim = fDim.ai,
                 maxEval = maxEval.ai, absError=absError.ai, 
                 doChecking=doChecking.ai,
                 vectorInterface = FALSE)

}
