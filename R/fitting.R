#' Fitting Semi-parametric Factor Analysis Model 
#' 
#' \code{spfa} fits a unidimensional or two-dimension factor analysis spfa model 
#' using penalized maximum likelihood estimation. A unidimensional spfa model can
#' handle discrete response data (i.e., item responses including binary responses
#' and polytomous responses) or continuous response data (e.g., response time). 
#' A two-dimensional \code{spfa} model can only handle simple structure model 
#' with two latent factors load to continuous and discrete response data, respectively.
#' @param data a matrix that consists of item responses with missing data coded as \code{NA}.
#' @param dimension a vector of integers containing indicators of the latent factor. The default is \code{rep(0, ncol(data))} indicating all item load on the same latent factor. 
#' @param discrete a vector of \code{TRUE} or \code{FALSE} indicating whether the item is discrete. \code{TRUE}: discrete variable. The length of the vector should be the same as the number of columns of the input data. 
#' @param control a list containing technical parameters for estimation. May be: 
#'   \itemize{
#'     \item \code{shortpar} a list containing the starting values of spfa model parameters for each item.
#'     \item  \code{pos} a list containing positivity constraints
#'      \item  \code{lmbd} a vector of the penalty parameter (lambda). Default value is a vector of 1s. 
#'      \item \code{n_basis} number of basis. Default is 11. 
#'      \item  \code{n_quad} number of quadrature points. Default is 21
#'      \item  \code{maxit_em} the maximum number of iterations for the EM cycles. Default is 500.
#'      \item \code{maxit_mstep} the maximum number of iterations for the mstep optimizer.
#'       \item \code{maxit_start} 
#'       \item \code{tol_em} tolerance for the EM convergence. Default is 1e-4.
#'       \item \code{tol_mstep} tolerance for the m-step optimizer. Default is 1e-6.
#'       \item \code{n_thrd} number of cores used for the penalized EM algorithm to run. Default is 1. 
#'   }    
#'
#' @return a list including spfa model parameter estimates and marginal log-likelihood.
#' @references
#' Liu, Y., & Wang, W. (2022). Semiparametric Factor Analysis for Item-Level Response Time Data. \emph{Psychometrika, 87}(2), 666–692. \url{ https://doi.org/10.1007/s11336-021-09832-8}
#' @export spfa
#' 
#' @examples
#'  # load item response time data 
#'  RT <- spfa::simdata[,1:8]
#'  
#'  # Fit a unidimensional spfa model with continuous data (Response time)
#' 
#' 
#' \dontrun{
#' spfa_example <- spfa(data = RT, 
#'        dimension = rep(0, ncol(RT)), 
#'        discrete = rep(F, ncol(RT)))
#'        }
#'  
#'  # In the spfa pacakge, the output of spfa_example can be directly extracted. See eample code below:
#'  
#'  spfa:::spfa_example$shortpar
#'  
#'  # Visualize the result for item 1 as an example 
#'  
#'  plotitem.cont(spfa:::spfa_example$par[[1]])

spfa <- function(
        data,                            
        dimension = rep( 0, ncol(data) ), 
        discrete = rep(F, ncol(data) ),  
        control = list()               
)
{
    # check validity of arguments
    m <- ncol(data)
    if (length(dimension) != m)
        stop("length(dimension) is not equal to ncol(data).")
    uniq.dim <- unique(dimension)
    d <- length(uniq.dim)
    if (d != 1 && d != 2)
        stop("can only handle one- or two-dimensional models for now.")
    dim <- as.numeric( factor(dimension) ) - 1
    if (length(discrete) != m)
        stop("length(discrete) is not equal to ncol(data).")
    disc <- as.numeric( as.logical(discrete) )
    
    # data manipulation
    ret <- list()
    ret$rng <- matrix(NA, 2, m)
    for ( j in seq_len(m) )
    {
        if (disc[j])  # discrete
        {
            data[, j] <- as.numeric( factor(data[, j]) ) - 1
            ret$rng[, j] <- range(data[, j], na.rm = T)
        }
        else  # continuous
        {
            ret$rng[, j] <- range(data[, j], na.rm = T)
            data[, j] <- (data[, j] - ret$rng[1, j]) / (ret$rng[2, j] - ret$rng[1, j])
        }
        data[, j][is.na(data[, j])] <- -999
    }
    
    # default control configuration
    ctrl <- list(
        shortpar = NULL, 
        pos = NULL, 
        lmbd = rep( 1, ifelse(d == 1, m, m + 1) ), 
        n_basis = 11, 
        n_quad = 21, 
        maxit_em = 500, 
        maxit_mstep = 20, 
        maxit_start = 20, 
        tol_em = 1e-4, 
        tol_mstep = 1e-6, 
        tol_start = 1e-6,
        n_thrd = 1)
    names.ctrl <- names(ctrl)
    names.user <- names(control)
    ctrl[names.user] <- control
    if ( length(invalid <- names.user[!names.user %in% names.ctrl]) ) 
        warning("unknown options in control: ", paste(invalid, collapse = ", ") )
    
    # check starting values and positivity constraints
    if( is.null(ctrl$shortpar) )
        ctrl$shortpar <- vector("list", m)
    if( is.null(ctrl$pos) )  # positive constraints
        ctrl$pos <- vector("list", m)
    for (j in seq_len(m))
    {
        if (disc[j])  # discrete items
            nshortpar <- ctrl$n_basis * ret$rng[2, j]
        else  # continuous items
            nshortpar <-  ctrl$n_basis * (ctrl$n_basis - 1)
        ctrl$shortpar[[j]] <- rep(0.1, nshortpar)
        ctrl$pos[[j]] <- rep(0, nshortpar)
    }
    if (d == 2)  # group
        ctrl$shortpar[[m + 1]] <- rep(0.1, ctrl$n_basis * ctrl$n_basis)
    
    # estimation
    ret <- c(ret, spfa_main2(
        # arguments that users cannot set via control
        dat = data,
        na = -999,
        update_group = ifelse(d == 1, F, T),
        item_type = disc,
        dim = dim,
        # arguments that users can set via control
        shortpar = ctrl$shortpar,
        pos = ctrl$pos, 
        n_basis = ctrl$n_basis, 
        lmbd = ctrl$lmbd, 
        n_quad = ctrl$n_quad, 
        maxit_em = ctrl$maxit_em, 
        maxit_mstep = ctrl$maxit_mstep, 
        maxit_start = ctrl$maxit_start, 
        tol_em = ctrl$tol_em, 
        tol_mstep = ctrl$tol_mstep, 
        tol_start = ctrl$tol_start,
        n_thrd = ctrl$n_thrd) )
    
    # compute marginal log-likelihood
    ret$loglik <- marg_loglik2(
        # arguments that users cannot set via control
        dat = data,
        na = -999,
        update_group = ifelse(d == 1, F, T),
        item_type = disc,
        dim = dim,
        shortpar = ret$shortpar,
        # arguments that users can set via control
        n_basis = ctrl$n_basis,
        n_quad = ctrl$n_quad,
        n_thrd = ctrl$n_thrd)
    
    # return
    ret
}
