#**************************************************************************
#* LFMC - A remote sensing approach to predict LFMC at large scale
#* Quan time-series filter for outlier detection
#* 
#* Inputs:
#*  x = numeric vector representing a time series
#*  k = window length 2*k+1 in indices
#*  
#* Author: Angel Cunill Camprubi <acunill86@gmail.com>
#* Last update: 11 Aug 2021
#**************************************************************************

quan <- function(x, k) {
    n <- length(x)
    y <- x
    S0 <- sd(x)
    val <- rep(0, n)
    for (i in (k + 1):(n - k)) {
        x0 <- median(x[(i - k):(i + k)])
        iden <- abs(x[i] - x0) / S0
        val[i] <- iden
    }
    val
}
