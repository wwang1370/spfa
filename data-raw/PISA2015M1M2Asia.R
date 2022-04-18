## code to prepare `PISA2015M1M2Asia` dataset goes here

usethis::use_data(PISA2015M1M2Asia, overwrite = TRUE)
# clean PSIA2015M1M2Asian.csv

# Example R code
#
# Last modified: 05/31/2022

#----------------------------------------
# Basic fitting functions
#----------------------------------------
# 1. data preprocessing (only need to implement the rescaling part)
fulldata <- read.csv("data-raw/PISA2015M1M2Asia.csv")
sel.RT <- grepl("^CM.*T$", names(fulldata) )  # response time
sel.R <- grepl("^CM.*S$", names(fulldata) )  # responses
RT <- as.matrix(fulldata[, sel.RT]) 
R <- as.matrix(fulldata[, sel.R])
# create testlets
#to.rm <- c(3, 4, 5, 6, 11, 12, 16, 17)
# responses
R[, 3] <- R[, 3] + 2 * R[, 4]
R[, 5] <- R[, 5] + 2 * R[, 6]
R[, 11] <- R[, 11] + 2 * R[, 12]
R[, 16] <- R[, 16] + 2 * R[, 17]
R <- R[, -c(4, 6, 12, 17)]
# RT
RT[, 3] <- RT[, 3] + RT[, 4]
RT[, 5] <- RT[, 5] + RT[, 6]
RT[, 11] <- RT[, 11] + RT[, 12]
RT[, 16] <- RT[, 16] + RT[, 17]
RT <- RT[, -c(4, 6, 12, 17)]
#RT <- log10(RT)
nobsn <- nrow(RT)
nitem <- ncol(RT)

# delete extreme RT
for ( j in seq_len(nitem) )
{
    cutj <- median(RT[, j]) + 4.4478 * mad(RT[, j], constant = 1)
    outj <- RT[, j] > cutj
    #cutj <- quantile( RT[, j], prob = c(0.025, 0.975) )
    #outj <- RT[, j] < cutj[1] | RT[, j] > cutj[2]
    RT[outj, j] <- -999  # missing code
    R[outj, j] <- -999  # missing code
    #RT[outj, j] <- NA  # missing code
    #R[outj, j] <- NA  # missing code
    rgj <- range(RT[!outj, j])
    RT[!outj, j] <- (RT[!outj, j] - rgj[1]) / diff(rgj)  # scale to [0, 1]
}



