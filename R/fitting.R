
#' Semi-parametric Factor Analysis (spfa) Model 
#' 
#' 
#' 
#' \code{spfa} fits a unidimensional or two-dimension factor analysis model using penalized maximum likelihood estimation. A unidimensional spfa model can handle  
#  discrete response data (i.e., item responses including binary responses and polytomous responses)
#' or continuous response data (e.g., response time). A two-dimensional spfa model can only handle simple structure model with two latent factors loads to continuous and discrete response data, respectively.
#' @param data a matrix that consists of item responses with missing data coded as \code{NA}.
#' @param dimension a vector of integers containing indicators of the latent factor. The default is rep(0, ncol(data)) indicating all item loads on the same latent factor. 
#' @param discrete a vector of \code{TRUE} or \code{FALSE} indicating whether the item is discrete with \code{TRUE} indicating the item is a discrete variable. The length of the vector should be the same as the number of columns of the input data. 
#' @param control a list containing technical parameters for estimation. For example:
#'   \itemize{
#'     \item \code{shortpar} a list containing the starting values of spfa model parameters for each item.
#'     \item  \code{pos} a list containing the positivity constraints
#'      \item  \code{lmbd} a vector of penalty parameter (lambda). Default values a vector of 1s. 
#'      \item \code{n_basis} number of basis. Default is 11. 
#'      \item  \code{n_quad} number of quadrature points. Default is 21
#'      \item  \code{maxit_em} the maximum number of iterations for the EM cycles. Default is 500.
#'      \item \code{maxit_mstep} the maximum number of iterations for the mstep optimizer.
#'       \item \code{maxit_start} 
#'       \item \code{tol_em} tolerance for the EM convergence. Default is 1e-4.
#'       \item \code{tol_mstep} tolerance for the mstep optimizer. Default is 1e-6.
#'       \item \code{n_thrd} number of cores used for the penalized EM algorithm to run. Default is 1. 
#'    
#'   
#'   }        
#'
#' @return
#' @export
#'
#' @examples 
#' # load simulated data 
#' 
#' 
#' #RT <-spfa::simdata[,1:8]
#' #R <-spfa::simdata[,9:16]
#' 
#' 
#' # Fit a unidimensional spfa model with discrete data 
#' 
#' 
#' \dontrun{spfa(data = R, dimension = rep(1, ncol(R)), discrete = rep(T, ncol(R)))}
#' 
#' 
#' # Fit a unidimensional spfa model with continuous data
#' 
#' 
#' \dontrun{spfa(data = RT, dimension = rep(0, ncol(RT)), discrete = rep(F, ncol(RT)))}
#' 
#' 
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
    item_type = disc,
    shortpar = NULL, 
    pos = NULL, 
    lmbd = rep( 1, ifelse(d == 1, m, m + 1) ), 
    n_basis = 11, 
    n_quad = 21, 
    dim = dim,
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
    # arguments that users cannot set
    dat = data,
    na = -999,
    update_group = ifelse(d == 1, F, T),
    # arguments that users can set
    item_type = ctrl$item_type,
    shortpar = ctrl$shortpar,
    pos = ctrl$pos, 
    n_basis = ctrl$n_basis, 
    lmbd = ctrl$lmbd, 
    n_quad = ctrl$n_quad, 
    dim = ctrl$dim,
    maxit_em = ctrl$maxit_em, 
    maxit_mstep = ctrl$maxit_mstep, 
    maxit_start = ctrl$maxit_start, 
    tol_em = ctrl$tol_em, 
    tol_mstep = ctrl$tol_mstep, 
    tol_start = ctrl$tol_start,
    n_thrd = ctrl$n_thrd) )
   
  ## compute marginal log-likelihood
  #ret$loglik <- marg_loglik2(
  #  dat = data,
  #  na = -999,
  #  update_group = ifelse(d == 1, F, T),
  #  item_type = ctrl$item_type,
  #  shortpar = ret$shortpar,
  #  dim = ctrl$dim,
  #  n_basis = ctrl$n_basis,
  #  n_quad = ctrl$n_quad,
  #  n_thrd = ctrl$n_thrd)

  # return
  ret
}
