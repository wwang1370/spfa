# Fit two-dimensional independent-cluster SPFA
spfa_joint <- function(dat = NULL, RT = NULL, 
                          resultRT = NULL,
                          resultR = NULL,
                           na.string = -999, 
                           control = list()){
    if(ncol(dat)!=ncol(RT)){
        abort("Please check the item response data and item response time data input. The number of items does not match in the two datasets.")
    }
    
    nitem = ncol(dat)
    
    ctrl <- list(
        item_type = c( rep(0, nitem), rep(1, nitem) ),
        shortpar = NULL, 
        pos = NULL, 
        lmbdR = 10, 
        lmbdRT = 1,
        lmbdJoint = 1,
        n_basis = 11L, 
        n_quad = 21L, 
        dim = c( rep(0, nitem), rep(1, nitem) ),
        maxit_em = 10000, 
        maxit_mstep = 20, 
        maxit_start = 20, 
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
    ctrl$shortpar <- ctrl$pos <- vector("list", 2*nitem)
    for (j in 1:nitem)  # RT
    {
        nshortpar <- ctrl$n_basis * (ctrl$n_basis - 1)
        ctrl$shortpar[[j]] <- resultRT$shortpar[[j]]   # start from separate fitting result
        ctrl$pos[[j]] <- rep(0, nshortpar)
        if (j == 1) ctrl$pos[[j]][-(1:(ctrl$n_basis - 1))] <- 1
    }

    for (j in nitem + 1:nitem)  # responses
    {
        ncat <- length( unique(dat[, j - nitem]) ) - 1
        nshortpar <- ctrl$n_basis * (ncat - 1)
        ctrl$shortpar[[j]] <- resultR$shortpar[[j - nitem]] # start from separate fitting result
        ctrl$pos[[j]] <- rep(0, nshortpar)
        if (j == nitem + 1) ctrl$pos[[j]][-1] <- 1
    }
    ctrl$shortpar[[2 * nitem + 1]] <- 
        rep(1, ctrl$n_basis * ctrl$n_basis)
    
    }
    result <- spfa_main2(
        dat = cbind(RT, dat),
        na = na.string,
        item_type = ctrl$item_type,
        shortpar = ctrl$shortpar, 
        pos = ctrl$pos, 
        n_basis = ctrl$n_basis, 
        lmbd = c( rep(ctrl$lmbdRT,nitem), rep(ctrl$lmbdR, nitem), ctrl$lmbdJoint), 
        n_quad = ctrl$n_quad, 
        dim = ctrl$dim,
        update_group = T,
        maxit_em = ctrl$maxit_em, 
        maxit_mstep = ctrl$maxit_mstep, 
        maxit_start = ctrl$maxit_start, 
        tol_em = ctrl$tol_em, 
        tol_mstep = ctrl$tol_mstep, 
        tol_start = ctrl$tol_start,
        n_thrd = ctrl$n_thrd) 
    
    return(result)
    
}
