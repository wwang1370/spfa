# functions in development...

# function to transform data into U[0,1]

#' marginal probability integral transform
#' 
#' Transform a random variable into a scale of 0 and 1 using the marginal 
#' probability integral transformation
#' 
#' @param data Response time data or log response time data on a scale different from 0 to 1.
#' @param logT a logical value indicating whether response time data is on a log scale. The default is \code{TRUE}
#' @param na.strings a character vector of strings which are to be interpreted as NA values. Has to be outside the range of 0 to 1 to be recognized as a missing value to be used in the \code{spfa} package. The default is -999.
#' @return return response time data scaled from 0 to 1 with -999 as the missing value. 
#' @export
#' @examples

DenTrans <- function(data, logT = T, na.strings = -999){
        
        logrt <- if(logT == T){ data }else{ log(data)}
        
        nobsn <- nrow(logrt) 
        nitem <- ncol(logrt)
        
        mdns <- vector("list", nitem)
        tlogrt <- logrt
        for (j in 1:nitem)
        {
                mdns[[j]] <- density(logrt[, j], na.rm = T)
                tlogrt[, j] <- spatstat::CDF(mdns[[j]])(logrt[, j])
        }
        tlogrt[is.na(tlogrt)] <- -999
        
        return(tlogrt)
}



#' min max method 
#' 
#' To scale a manifest variable into the range from 0 to 1. 
#' 
#' @param data Data frame. Not necessarily on the [0, 1] scale.
#' @param na A numeric value which is to be interpreted as NA values. The default is -999.
#' @return
#' @export
#'
#' @examples

minmax <- function(data, na = -999)
{
  data[data == na] <- NA
  tdata <- apply( data, 2L, function(y) 
    ( y - min(y, na.rm = T) ) / diff( range(y, na.rm = T) ) )
  attr(tdata, "range") <- apply(data, 2L, range, na.rm = T)
  tdata[ is.na(tdata) ] <- na
  tdata
}

#' multi-fold cross-validation
#' 
#' Wrapper function to perform multi-fold cross-validation.
#' 
#' @param data Data frame. Not necessarily on the [0, 1] scale.
#' @param n.fold Number of folds The default is 10.
#' @param na A numeric value which is to be interpreted as NA values. The default is -999.
#' @param lambda Vector of penalty weights. The default is 1.
#' @param control List of control variables. 
#' @return
#' @export
#' @examples

cross.val <- function( 
  data, na = -999, n.fold = 10L, lambda = 1, control = list()
)
{
  # default control list
  n <- nrow(data)
  m <- ncol(data)
  maxit.start <- 20L
  tdata <- minmax(data, na = na)  # transformed data
  lambda <- sort(lambda, decreasing = T)  # decreasing lambda
  n.lambda <- length(lambda)
  ctrl <- list(
    dat = NULL,
    item_type = rep(0L, m),
    shortpar = NULL, 
    pos = NULL, 
    n_basis = 11L, 
    lmbd = lambda[1L], 
    n_quad = 21L, 
    dim = rep(1L, m),
    maxit_em = 5000L, 
    maxit_mstep = 20L, 
    maxit_start = maxit.start, 
    tol_em = 1e-4, 
    tol_mstep = 1e-6, 
    tol_start = 1e-6,
    na = na,
    n_thrd = 1L)
  start.val <- ctrl$pos <- vector("list", m)
  for ( j in seq_len(m) )
  {
    n.shortpar <- ctrl$n_basis * (ctrl$n_basis - 1L)
    start.val[[j]] <- rep(0.1, n.shortpar)
    ctrl$pos[[j]] <- rep(0, n.shortpar)
    if (j == 1L) 
      ctrl$pos[[j]][-seq_len(ctrl$n_basis - 1)] <- 1
  }

  # user-supplied control variables
  names.ctrl <- names(ctrl)
  names.user <- names(control)
  ctrl[names.user] <- control
  if ( length(invalid <- names.user[!names.user %in% names.ctrl]) ) 
    warning("unknown options in control: ", paste(invalid, collapse = ", ") )
  if ( !is.null(control$maxit_start) ) 
    maxit.start <- control$maxit_start  # save this value for later use

  # multi-fold cross-validation
  fold <- split(sample(1L:n), 1L:n.fold)  # indices for validation data in each fold
  risk <- matrix(0, n.lambda, n.fold)  # KL risk
  for ( k in seq_len(n.fold) )  # loop over folds
  {
    fit <- NULL
    ctrl$dat <- tdata[-fold[[k]], ]  # calibration data
    nv <- sum( apply(tdata[fold[[k]], ], 1L, function(y) !all( is.na(y) ) ) )
    for ( l in seq_along(lambda) )
    {
      cat("***** Fold ", k, ", ", "lambda = ", lambda[l], " *****\n", sep = '')
      # controls
      ctrl$lmbd <- lambda[l]
      if (l > 1) 
      {
        ctrl$shortpar <- fit$shortpar
        ctrl$maxit_start <- 0L
      }
      else 
      {
        ctrl$shortpar <- start.val
        ctrl$maxit_start <- maxit.start
      }
      # fit
      fit <- do.call(spfa_main, args = ctrl)
      # risk computation with validation data
      risk[l, k] <- - marg_loglik1(
        tdata[fold[[k]], ], 
        na = ctrl$na,
        item_type = ctrl$item_type,
        shortpar = fit$shortpar,
        n_basis = ctrl$n_basis,
        n_quad = ctrl$n_quad,
        n_thrd = ctrl$n_thrd) / nv + 
        sum( log( apply(attr(tdata, "range"), 2L, diff) ) )
    }
  }

  # return
  risk
}
