#' Item level perspective plots or contour plots for spfa models
#' 
#' For continuous response data use \code{plotitem.cont} whereas discrete response data use \code{plotitem.disc}. For joint continuous and discrete data, use \code{plotgroup}.
#' 
#' @param param parameter vector estimated from \code{\link{spfa}} model
#' @param nquad an integer value of number of quadrature points. Default is 21
#' @param npoints an integer value of number of x and y levels in the plot
#' @param xlim the x limits of the plot. Two numerical values indicating the lower and upper limits
#' @param ylim the y limits of the plot. Two numerical values indicating the lower and upper limits of the density. Note y is rescaled to a uniform [0,1] distribution.
#' @param normal a logical value \code{TRUE} or \code{FALSE} indicating which density is used to rescale y.
#' @param FUN a user supplied function to rescale.
#' @param plot a logical value \code{TRUE} or \code{FALSE} indicating whether the plot is visualized. 
#' @param type the type of plot to be visualized. The default is the contour plot \code{\link{contour}}. It can also be changed to "\code{persp}" indicating perspective plots.
#' @param ... additional arguments passed to \code{\link{contour}} and \code{\link{persp}}
#' @seealso \code{\link{contour}} and \code{\link{persp}}
#' @return plots. Item level perspective and contour plot
#' @export
#' @examples
#' 
#' # Contour plot of the first item 
#' 
#' plotitem.cont(spfa::spfa_example$par[[1]])


plotitem.cont <- function(
  param,                    
  nquad = 21,               
  npoints = 101,            
  xlim = c(-2.5, 2.5),      
  ylim = c(0, 1),           
  normal = TRUE,               
  FUN = NULL,              
  plot = TRUE,                
  type = "contour",         
  ...                      
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
    plot <- FALSE  # automatically turn off plotting
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

#' # Item level plot for discrete response data:
#' @rdname plotitem.cont
#' @param ncat an integer value indicating the number of categories for the discrete item.
#' @param col color of the line. 
#' @param lty line type

plotitem.disc <- function(
  param,                    # parameter vector
  ncat,                     # number of categories (integer)
  npoints = 101,            # number of x levels (integer)
  xlim = c(-2.5, 2.5),      # limit of x (numeric, length = 2)
  normal = T,               # normal scale? (logical)
  FUN = NULL,               # function to integrate ( function(y) )
  plot = TRUE,                 # create plot? (logical)
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
    plot <- FALSE  # automatically turn off plotting
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
#' # Item level plot for continuous and discrete response data
#' @param lim limit
#' @rdname plotitem.cont

plotgroup <- function(
  param,                   
  nquad = 21,               
  npoints = 101,           
  lim = c(-2.5, 2.5),      
  normal = TRUE,               
  plot = TRUE,                
  type = "contour",        
  ...                      
)                     
{
  # check validity of arguments
  q <- length(param)
  nbasis <- as.integer( sqrt(q) )
  if (nbasis * nbasis != q)
    stop("length(param) is cannot be factorized into nbasis * nbasis")
  if (length(lim) != 2)
    stop("length(lim) must be 2")
  npoints <- max(4, npoints)  # npoints >= 4
  nquad <- max(2, nquad)  # nquad >= 2

  # computation
  ret <- list(xy = NULL)
  if (normal)  # normal grid
  {
    ret$xy <- seq(lim[1], lim[2], , npoints)
    B <- bspl(pnorm(ret$xy), nbasis, 4, 0, 1)
  }
  else  # uniform grid
  {
    lim <- c(0, 1)
    ret$xy <- seq(0, 1, , npoints)
    B <- bspl(ret$xy, nbasis, 4, 0, 1)
  }
  parmat <- matrix(param, nbasis, nbasis)
  ret$dns <- B %*% parmat %*% t(B)
  if (normal)
    ret$dns <- ret$dns * tcrossprod( dnorm(ret$xy) )
  if (!plot) return(ret)

  # plotting
  if (type == "contour")
  {
    contour(ret$xy, ret$xy, ret$dns, xlim = lim, ylim = lim,
      xlab = expression( italic(x)[1] ), ylab = expression( italic(x)[2] ), ...)
  }
  else if (type == "persp")
  {
    persp(ret$xy, ret$xy, ret$dns, xlim = lim, ylim = lim, 
      xlab = "x1", ylab = "x2", 
      zlab = "Density", ...)
  }
  else
    stop("type is not yet supported")
}
