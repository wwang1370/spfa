# Fit one-dimensional SPFA for discrete data (item responses)
#' Fit one-dimensional SPFA for discrete data (item responses)
#'
#' @param dat A matrix of item response data which contains \code{ncol(dat)} of discrete item responses.
#' @param na.string A numeric value which is to be interpreted as NA values. The default is -999.
#' @param control A list of control variables.
#'
#' @return
#' @export
#'
#' @examples
spfa_discrete <- function(dat = NULL, 
                           na.string = -999, 
                           control = list(n_basis = 11, lmbd = 10)){
    nitem <- ncol(dat)
    control$item_type <- rep(1, nitem)
    for (j in seq_len(nitem))
    {
        ncat <- length( unique(dat[, j]) ) - 1
        nshortpar <- control$n_basis * (ncat - 1)
        control$shortpar[[j]] <- rep(0.1, nshortpar)
        control$pos[[j]] <- rep(0, nshortpar)
        if (j == 1) control$pos[[j]][-1] <- 1
    }
    
   spfa_continous(RT = dat, 
                   na.string = na.string,
                   control = control)
    
    
}
