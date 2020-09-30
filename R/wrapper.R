# additional wrapper functions
#
# Author: Yang Liu
#
# Last modified: 09/29/2020 */

# xv_risk: cross-validation risk function

xv_risk <- function(
  shortpar,   # parameters (numeric, dim = n_shortpar)
  dat,  # cross validation data (numeric, dim = n_obsn x n_item)
  n_basis,  # number of basis functions (integer)
  n_quad  # number of quadrature points (integer)
)
{
  # initialize
  n_obsn <- nrow(dat)
  n_item <- ncol(shortpar)
  quad <- gl_quad(n_quad, 0, 1)  # quadrature
  capture.output( sco <- spfa_score(shortpar, dat, n_basis, n_quad) )  # EAP score
  # main loop over items
  risk <- 0
  for ( j in seq_len(n_item) )
  {
    cat("Computing risk for item ", j, '\r', sep = '')
    capture.output( fitted <- crossprod( cond_dns(shortpar[, j], 
      quad$x, sco, n_basis, n_quad), quad$x * quad$w) )
    risk <- risk + mean( (dat[, j] - fitted)^2 )  # mean squared error
  }
  cat('\n')
  risk / n_item  # take average
}

 ## 2nd term
 #capture.output( risk <- risk - 2 * marg_lik(shortpar[, j, drop = F], 
 #  dat[, j, drop = F], n_basis, n_quad) / n_obsn )
 ## 1st term
 #for ( q in seq_len(n_quad) )
 #{
 #  capture.output( risk <- risk + marg_lik(shortpar[, j, drop = F], 
 #    matrix(quad$x[q]), n_basis, n_quad)^2 * quad$w[q] )
 #}
