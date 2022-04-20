# plotting functions
#
# -Weimeng: Please write the help document; after that you can delete my
# comments on the arguments 

plotitem.cont <- function(
  param = NULL        # parameter vector
)                     

plotitem.disc <- function(
  param = NULL,             # parameter vector
  n.cat = NULL,             # number of categories (integer)
  n.xlevels = 101,          # number of x levels (integer)
  normal = T,               # normal scale? (logical)
  plot = T,                 # create plot? (logical)
  col = 1:n.cat,            # default color
  lty = 1:n.cat,            # default line type
  ...                       # additional arguments
)
{
  # check validity of arguments
  q <- length(param)
  if (q %% n.cat != 0)
    stop("length(param) is not a multiple of n.cat")
  n.basis <- q / n.cat - 1
  if (length(col) != n.cat)
    stop("length(col) is not equal to n.cat")
  if (length(lty) != n.cat)
    stop("length(lty) is not equal to n.cat")
  n.xlevels <- max(4, n.xlevels)  # n.xlevels >= 4

  # computation
  x <- seq(0, 1, , n.xlevels)  # uniform grid
  if (normal)  # normal grid
  {
    lb <- qnorm(x[2])
    x <- seq(lb, -lb, , n.xlevels)
  }
  Bx <- bspl(x, n.basis, 4, 0, 1)
  parmat <- matrix(par, , nbx + 1L)
  eta <- exp( parmat[, 1L] + tcrossprod(parmat[, -1L], Bx) )
  prob <- t( apply( eta, 2L, function(x) x / sum(x) ) )
  if (!plot) return(prob)

  # plotting
}
