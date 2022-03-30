
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mvpd: An R Package for Multivariate Product Distributions

-   `[dpr]mvss`: mutivariate subgaussian stable distributions

<!-- badges: start -->
<!-- badges: end -->

The goal of `mvpd` is to use product distribution theory to allow the
numerical calculations of specific scale mixtures of the multivariate
normal distribution. The multivariate subgaussian stable distribution is
the product of the square root of a univariate positive stable
distribution and the multivariate normal distribution (see Nolan
(2013)).

<!-- ## Installation -->
<!-- You can install the development version of mvpd from [GitHub](https://github.com/) with: -->
<!-- ``` r -->
<!-- # install.packages("devtools") -->
<!-- devtools::install_github("swihart/mvpd") -->
<!-- ``` -->

## Example

Generate 1000 draws from a random bivariate subgaussian stable
distribution with alpha=1.71 and plot.

``` r
library(mvpd)
set.seed(10)
## basic example code
biv <- rmvss(n=1e3, alpha=1.71, Q=matrix(c(10,7.5,7.5,10),2))
head(biv)
#>            [,1]       [,2]
#> [1,] -0.2260798 -0.6168492
#> [2,] -6.1460819 -4.5603538
#> [3,]  1.4592466  1.6213040
#> [4,] -4.4159078 -2.9252448
#> [5,] -6.7106973 -3.8158068
#> [6,]  5.9107788  5.1332625
plot(biv); abline(h=0,v=0)
```

<img src="man/figures/README-example-1.png" width="100%" />
