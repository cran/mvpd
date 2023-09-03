#' Multivariate Product Distributions
#'
#' The purpose of this package is to offer density, probability, and 
#' random variate generating (abbreviated as \[d/p/r\], respectively) 
#' functions for 
#' multivariate distributions that can be represented as a product distribution.
#' Specifically, the package will primarily focus on the product of a multivariate 
#' normal distribution and a univariate random variable. 
#' These product distributions are called Scale Mixtures of 
#' Multivariate Normal Distributions, 
#' and for particular choices of the univariate random variable distribution the 
#' resultant product distribution may be a family of interest.  For instance,
#' the square-root of a positive stable random variable multiplied by a 
#' multivariate normal distribution is the multivariate subgaussian stable
#' distribution. Product 
#' distribution theory is applied for implementing their computation.
#' 
#'
#' @section Multivariate subgaussian stable distributions:
#'
#' \code{\link{dmvss}} -- multivariate subgaussian stable distribution density
#' 
#' \code{\link{pmvss}} -- multivariate subgaussian stable distribution probabilities
#' 
#' \code{\link{rmvss}} -- multivariate subgaussian stable distribution random variates
#'
#' \code{\link{pmvss_mc}} -- Monte Carlo version of probabilities, using \code{rmvss}
#' 
#' \code{\link{fit_mvss}} -- Fit a multivariate subgaussian stable distribution (e.g. estimate parameters given data)
#' 
#' @docType package
#' @name mvpd
"_PACKAGE"
#> NULL
