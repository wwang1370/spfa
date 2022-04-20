# fitting functions
#
# -Weimeng: Please write the help document; after that you can delete my
# comments on the arguments 

spfa <- function(
  data,                             # data set
  dimension = rep( 0, ncol(data) ), # dimension indicators (integer)
  discrete = rep(F, ncol(data) ),   # discrete indicator (logical)
  control = list()                  # control list
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
      nshortpar <- control$n_basis * ret$rng[2, j]
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
   
  # compute marginal log-likelihood
  ret$loglik <- marg_loglik2(
    dat = data,
    na = -999,
    update_group = ifelse(d == 1, F, T),
    item_type = ctrl$item_type,
    shortpar = ret$shortpar,
    dim = ctrl$dim,
    n_basis = ctrl$n_basis,
    n_quad = ctrl$n_quad,
    n_thrd = ctrl$n_thrd)

  # return
  ret
}
