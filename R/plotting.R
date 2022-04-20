# plotting functions
#
# -Weimeng: Please write the help document; after that you can delete my
# comments on the arguments 

plotitem.cont <- function(
  param,                    # parameter vector
  nquad = 21,               # number of quadrature points (integer)
  npoints = 101,            # number of x and y levels (integer)
  xlim = c(-2.5, 2.5),      # limit of x (numeric, length = 2)
  ylim = c(0, 1),           # limit of y (numeric, length = 2), density will be rescaled
  normal = T,               # normal scale? (logical)
  FUN = NULL,               # function to integrate ( function(y) )
  plot = T,                 # create plot? (logical)
  type = "contour",         # type of visualization
  ...                       # additional arguments passed to contour/persp
)                     
{
  # check validity of arguments
  q <- length(param)
  nbasis <- as.integer( 0.5 * (sqrt(1 + 4 * q) - 1) )
  if (nbasis * (nbasis + 1) != q)
    stop("length(param) is cannot be factorized into nbasis * (nbasis + 1)")
  if (length(xlim) != 2)
    stop("length(xlim) must be 2")
  if (length(ylim) != 2)
    stop("length(ylim) must be 2")
  npoints <- max(4, npoints)  # npoints >= 4
  nquad <- max(2, nquad)  # nquad >= 2

  # computation on the [0, 1] scale
  ret <- list(x = NULL, y = NULL)
  ret$y <- seq(0, 1, , npoints)
  By <- bspl(ret$y, nbasis, 4, 0, 1)
  if (normal)  # normal grid
  {
    ret$x <- seq(xlim[1], xlim[2], , npoints)
    Bx <- bspl(pnorm(ret$x), nbasis, 4, 0, 1)
  }
  else  # uniform grid
  {
    xlim <- c(0, 1)
    ret$x <- ret$y
    Bx <- By
  }
  quad <- gl_quad(nquad, 1, 0, 1)  # G-L quadrature
  quad$x <- c(quad$x)
  quad$w <- c(quad$w)
  Bq <- bspl(quad$x, nbasis, 4, 0, 1)
  hy <- drop(By %*% param[1:nbasis]) +
    tcrossprod(By %*% matrix(param[-(1:nbasis)], nbasis, nbasis), Bx)
  hq <- drop(Bq %*% param[1:nbasis]) +
    tcrossprod(Bq %*% matrix(param[-(1:nbasis)], nbasis, nbasis), Bx)
  ret$dns <- t( exp( sweep(hy, 2, log( colSums(exp(hq) * quad$w) ), `-`) ) )

  # rescale to ylim
  dy <- diff(ylim)
  ret$y <- ret$y * dy + ylim[1]  # rescale y
  ret$dns <- ret$dns / dy  # rescale density
  quad$x <- quad$x * dy + ylim[1]  # rescale quadrature nodes
  quad$w <- quad$w * dy  # rescale quadrature weights

  # if user supplies function to integrate
  if ( !is.null(FUN) )
  {
    plot <- F  # automatically turn off plotting
    fy <- sapply(quad$x, FUN)
    # recalculate dns at quad$x
    dns <- t( exp( sweep(hq, 2, log( colSums(exp(hq) * quad$w) ), `-`) ) )
    ret$f <- rowSums( sweep(dns, 2, fy * quad$w, `*`) )
  }
  if (!plot) return(ret)

  # plotting
  if (type == "contour")
  {
    contour(ret$x, ret$y, ret$dns, xlim = xlim, ylim = ylim, 
      xlab = expression( italic(x) ), ylab = expression( italic(y) ), ...)
  }
  else if (type == "persp")
  {
    persp(ret$x, ret$y, ret$dns, xlim = xlim, ylim = ylim, 
      xlab = 'x', ylab = 'y', zlab = "Density", ...)
  }
  else
    stop("type is not yet supported")
}

plotitem.disc <- function(
  param,                    # parameter vector
  ncat,                     # number of categories (integer)
  npoints = 101,            # number of x levels (integer)
  xlim = c(-2.5, 2.5),      # limit of x (numeric, length = 2)
  normal = T,               # normal scale? (logical)
  FUN = NULL,               # function to integrate ( function(y) )
  plot = T,                 # create plot? (logical)
  col = 1:ncat,             # default color
  lty = rep(1, ncat),       # default line type
  ...                       # additional arguments passed to plot
)
{
  # check validity of arguments
  q <- length(param)
  if (q %% ncat != 0)
    stop("length(param) is not a multiple of ncat")
  nbasis <- q / ncat - 1
  if (length(col) != ncat)
    stop("length(col) is not equal to ncat")
  if (length(lty) != ncat)
    stop("length(lty) is not equal to ncat")
  npoints <- max(4, npoints)  # npoints >= 4

  # computation
  ret <- list(x = NULL)
  if (normal)  # normal grid
  {
    ret$x <- seq(xlim[1], xlim[2], , npoints)
    Bx <- bspl(pnorm(ret$x), nbasis, 4, 0, 1)
  }
  else  # uniform grid
  {
    xlim <- c(0, 1)
    ret$x <- seq(0, 1, , npoints)
    Bx <- bspl(ret$x, nbasis, 4, 0, 1)
  }
  parmat <- matrix(param, , nbasis + 1)
  eta <- exp( parmat[, 1] + tcrossprod(parmat[, -1], Bx) )
  ret$prob <- t( apply( eta, 2, function(x) x / sum(x) ) )
  colnames(ret$prob) <- paste0("cat", 1:ncat - 1)
  if ( !is.null(FUN) )  # if user supplies function to integrate
  {
    plot <- F  # automatically turn off plotting
    fy <- sapply(1:ncat - 1, FUN)
    # e.g., if supply function(y) y, then it computes expected score curve
    ret$f <- rowSums( sweep(ret$prob, 2, fy, `*`) )
  }
  if (!plot) return(ret)

  # plotting
  plot(ret$x, ret$prob[, 1], type = 'n', xlim = xlim, ylim = 0:1, 
    xlab = expression( italic(x) ), ylab = "Probability", ...)
  for ( k in seq_len(ncat) )
    lines(ret$x, ret$prob[, k], col = col[k], lty = lty[k], ...)
}
