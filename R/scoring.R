#' Computing EAP scores
#'
#' @param data data to be scored
#' @param fit the function return from fitting a \code{spfa} model.
#' @param dimension a vector of integers containing indicators of the latent factor. The default is rep(0, ncol(data)) indicating all item loads on the same latent factor. 
#' @param discrete a vector of \code{TRUE} or \code{FALSE} indicating whether the item is discrete with \code{TRUE} indicating the item is a discrete variable. The length of the vector should be the same as the number of columns of the input data. 
#' @param normal a logical value \code{TRUE} or \code{FALSE} indicating which density is used to rescale y.
#' @param control a list of technical control variables. See \code{\link{spfa}}.
#'
#' @return EAP scores for the fitted spfa model and reliability
#' @export 
#' @examples
#' RT <- spfa::simdata[,1:8]
#' myeaps <- fscores(data = RT, fit = spfa:::spfa_example, 
#' dimension = rep( 0, ncol(RT)), discrete = rep(FALSE, ncol(RT) ))
fscores <- function(
  data,  # data to be scored                            
  fit,   # function return from spfa 
  dimension = rep( 0, ncol(data) ),  # dimension indicator
  discrete = rep(FALSE, ncol(data) ),    # discrete indicator
  normal = TRUE,   # normal scale?
  control = list()               
)
{
  # check validity of arguments
  m <- ncol(data)
  if (ncol(fit$rng) != m)
    stop("numbers of items in data and fit do not match.")
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
  for ( j in seq_len(m) )
  {
    if (disc[j])  # discrete
    {
      data[, j] <- as.numeric( factor(data[, j]) ) - 1
      if ( max(data[, j], na.rm = T) > fit$rng[2, j] )
        stop( paste0("responses to discrete item ", j, "out of bound") )
    }
    else  # continuous
    {
      if ( min(data[, j], na.rm = T) < fit$rng[1, j] || 
           max(data[, j], na.rm = T) > fit$rng[2, j] )
        warning( paste0("responses to continuous item ", j, "out of bound") )
      data[, j] <- (data[, j] - fit$rng[1, j]) / (fit$rng[2, j] - fit$rng[1, j])
    }
    data[, j][is.na(data[, j])] <- -999
  }

  # default control configuration
  ctrl <- list(
    lmbd = rep( 1, ifelse(d == 1, m, m + 1) ), 
    n_basis = 11, 
    n_quad = 21, 
    n_thrd = 1)
  names.ctrl <- names(ctrl)
  names.user <- names(control)
  ctrl[names.user] <- control
  if ( length(invalid <- names.user[!names.user %in% names.ctrl]) ) 
    warning("unknown options in control: ", paste(invalid, collapse = ", ") )
  
  # scoring
  ret <- list()
  capture.output( sco <- spfa_score2(
    # arguments that users cannot set via control
    dat = data,
    na = -999,
    dim = dim,
    update_group = ifelse(d == 1, FALSE, TRUE),
    item_type = disc,
    # arguments that users can set via control
    shortpar = fit$shortpar,
    n_basis = ctrl$n_basis,
    n_quad = ctrl$n_quad,
    n_thrd = ctrl$n_thrd,
    mode = as.numeric(normal)
  ) )
  cov.indx <- outer(1:d, 1:d, paste0)
  colnames(sco) <- c( paste0('x', 1:d), 
    paste0('v', cov.indx[lower.tri(cov.indx, diag = TRUE)]) )
  ret$eap <- sco[, 1:d, drop = FALSE]      # EAP scores
  ret$pcov <- sco[, -(1:d), drop = FALSE]  # posterior covariance

  # reliability
  ret$rel <- numeric(d)
  for ( k in seq_len(d) )
  {
    eap.indx <- paste0('x', k)
    pvar.indx <- paste0('v', k, k)
    var.eap <- var(ret$eap[, eap.indx])
    ex.pvar <- mean(ret$pcov[, pvar.indx])
    ret$rel[k] <- var.eap / (var.eap + ex.pvar)
  }

  ret
}
