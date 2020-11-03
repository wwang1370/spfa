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
#' @param data Response time data or log response time data on a scale different from 0 to 1.
#' @param logT a logical value indicating whether response time data is on a log scale. The default is \code{TRUE}
#' @param na.strings a character vector of strings which are to be interpreted as NA values. Has to be outside the range of 0 to 1 to be recognized as a missing value to be used in the \code{spfa} package. The default is -999.
#' @return
#' @export
#'
#' @examples

minmax <- function(data, logT = T, na.strings = -999){
        
        logrt <- if(logT == T){ data }else{ log(data)}
        
       tlogrt <-  apply(logrt,2, FUN = function(x) if(max(x, na.rm = T)==min(x, na.rm = T)){(x- min(x, na.rm = T))}else{
               (x- min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))})
      
        tlogrt[is.na(tlogrt)] <- -999
        return(tlogrt)
}
