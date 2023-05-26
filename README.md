
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spfa

<!-- badges: start -->
<!-- badges: end -->

Estimation, scoring, and plotting functions for the semi-parametric
factor model proposed by Liu & Wang \[(2022)\]
(<https://link.springer.com/article/10.1007/s11336-021-09832-8>) and Liu
& Wang [(2023)](https://arxiv.org/abs/2303.10079). Both the conditional
densities of observed responses given the latent factors and the joint
density of latent factors are estimated non-parametrically. Functional
parameters are approximated by smoothing splines, whose coefficients are
estimated by penalized maximum likelihood using an
expectation-maximization (EM) algorithm. E- and M-steps can be
parallelized on multi-thread computing platforms that support ‘OpenMP’.
Both continuous and unordered categorical response variables are
supported.

## Installation

You can install the development version of spfa from
[GitHub](https://github.com/), with:

``` r
# install.packages("devtools")
devtools::install_github("wwang1370/spfa")
```

## Example

Various examples and help files have been compiled using the knitr
package to generate HTML output, and are available on the package help
file. User contributions are welcome!

``` r
library(spfa)
## basic example code includes fitting an spfa model for response time
#  RT <- spfa::simdata[,1:8]
# spfa_example <- spfa(data = RT, 
#       dimension = rep(0, ncol(RT)), 
#       discrete = rep(F, ncol(RT)))
```

## Citation:

Liu, Y., & Wang, W. (2022). Semiparametric factor analysis for
item-level response time data. *Psychometrika*, *87* (2), 666-692.

Liu, Y., & Wang, W. (2023). What Can We Learn from a Semiparametric
Factor Analysis of Item Responses and Response Time? An Illustration
with the PISA 2015 Data. Retrieved from
[http://arxiv.org/abs/2303.10079](https://arxiv.org/abs/2303.10079)
