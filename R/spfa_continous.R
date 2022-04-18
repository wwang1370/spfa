

#' Fit one-dimensional SPFA for continuous data (response time)
#'
#' @param RT Response time data in a matrix format. 
#' @param na.string A numeric value which is to be interpreted as NA values. The default is -999.
#' @param control A list of control variables.
#'
#' @return 
#' @export
#'
#' @examples
spfa_continous <- function(RT = NULL, 
                           na.string = -999, 
                           control = list()){
    if (is.null(RT)){
        abort("Attempted to fit a unidimensional SPFA model without data input.")
    }
    if (is.matrix(RT)==F){ RT <- as.matrix(RT)}
    nitem <- ncol(RT)
    ctrl <- list(
        item_type = rep(0L, nitem),
        shortpar = NULL, 
        pos = NULL, 
        lmbd = 1, 
        n_basis = 11L, 
        n_quad = 21L, 
        dim = rep(1L, nitem),
        maxit_em = 10000, 
        maxit_mstep = 20, 
        maxit_start = 0, 
        tol_em = 1e-3, 
        tol_mstep = 1e-5, 
        tol_start = 1e-5,
        n_thrd = 4L)

    
   
    
    # user-supplied control variables
    names.ctrl <- names(ctrl)
    names.user <- names(control)
    ctrl[names.user] <- control
    if ( length(invalid <- names.user[!names.user %in% names.ctrl]) ) {warning("unknown options in control: ", paste(invalid, collapse = ", ") )}
    if(is.null(ctrl$shortpar)){
    ctrl$shortpar <- ctrl$pos <- vector("list", nitem)
    for (j in seq_len(nitem))
    {
        nshortpar <-  ctrl$n_basis * ( ctrl$n_basis - 1)
        ctrl$shortpar[[j]] <- rep(0.1, nshortpar)
        ctrl$pos[[j]] <- rep(0, nshortpar)
        if (j == 1) ctrl$pos[[j]][-(1:( ctrl$n_basis - 1))] <- 1
    }}
    
  
    item.names <- if(is.null(colnames(RT))){paste0("RT", 1: nitem)}else{colnames(RT)}

     
resultRT <- spfa_main(
    dat = RT,
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
    na = na.string,
    n_thrd = ctrl$n_thrd)

   names(resultRT$shortpar) <- names(resultRT$par) <- item.names
   
return(resultRT)


}