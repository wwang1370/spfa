#   persp.ecdns: perspective plot for estimated conditional density
# 
# args:
#   f: true conditional density (function, with arguments x and y)
# par: parameters (numeric, dim = p = py + py * px)
#   x: x values for plotting (numeric)
#   y: y values for plotting (numeric)
# ...: additional graphical parameters passed to persp()

#' Get the loglikelihood 
#' 
#' Get the loglikelihood from the fitted model to prepare for the density plots
#'
#' @param par item parameters; numeric, dim = p = py + py * px
#' @param By  basis for y
#' @param Bx  basis for x
#' @param Bq  basis for quadrature points
#' @param qwt quadrature weights
#'
#' @return log-likelihood 

#'
#' @examples
loglik.out <- function(par, By, Bx, Bq, qwt)
{
        nby <- ncol(By)
        hy <- drop(By %*% par[1L:nby]) +
                tcrossprod(By %*% matrix(par[-(1L:nby)], nby, ), Bx)
        hq <- drop(Bq %*% par[1L:nby]) +
                tcrossprod(Bq %*% matrix(par[-(1L:nby)], nby, ), Bx)
        sweep( hy, 2L, log( colSums(exp(hq) * qwt) ) )
}
lik.out.d <- function(par, Bx)  # discrete
{
        nbx <- ncol(Bx)
        Par <- matrix(par, , nbx + 1L)
        eta <- exp( Par[, 1L] + tcrossprod(Par[, -1L], Bx) )
        t( apply( eta, 2L, function(x) x / sum(x) ) )
}
#' perspective plot for estimated conditional density
#'
#' @param par item parameters; numeric, dim = p = py + py * px
#' @param x  x values for plotting (numeric)
#' @param y  x values for plotting (numeric)
#' @param nby number of bases for y
#' @param nbx number of bases for x
#' @param ... additional parameters for perspective plots
#'


#' @seealso \code{\link{persp}}, \code{\link{contour}}, and \code{\link{persp.ecdns}}
#' @examples
persp.ecdns <- function(par, x = seq(0,1,0.02), y = seq(0,1,0.02), nby = 11, nbx = 11, ...)
{
        # set up basis
        nx <- length(x)
        ny <- length(y)
        Bx <- cubic_bspl(x, nbx)
        By <- cubic_bspl(y, nby)
        Bq <- cubic_bspl(rule.le$x, nby)
        # log logistic transform
        z <- t( exp( loglik.out(par, By, Bx, Bq, rule.le$w) ) )
        persp(x, y, z, ticktype = "detailed", ...)
}

#' perspective plot for estimated conditional density
#'
#' @param par item parameters; numeric, dim = p = py + py * px
#' @param x  x values for plotting (numeric)
#' @param y  x values for plotting (numeric)
#' @param nby number of bases for y
#' @param nbx number of bases for x
#' @param ... additional parameters for perspective plots
#'

#' @seealso \code{\link{persp}},\code{\link{contour}}, and \code{\link{persp.ecdns}}
#' @examples
contour.ecdns <- function(par, x = seq(0,1,0.02), y = seq(0,1,0.02), nby = 11, nbx = 11, ...)
{
        # set up basis
        nx <- length(x)
        ny <- length(y)
        Bx <- cubic_bspl(x, nbx)
        By <- cubic_bspl(y, nby)
        Bq <- cubic_bspl(rule.le$x, nby)
        # log logistic transform
        z <- t( exp( loglik.out(par, By, Bx, Bq, rule.le$w) ) )
        contour(x, y, z, ...)
}