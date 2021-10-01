#**************************************************************************
#* LFMC - A remote sensing approach to predict LFMC at large scale
#* Hampel identifier function
#* Modified from pracma::hampel
#* 
#* Inputs:
#*  x = numeric vector representing a time series
#*  k = window length 2*k+1 in indices
#*  
#* Author: Angel Cunill Camprubi <acunill86@gmail.com>
#* Last update: 11 Aug 2021
#**************************************************************************

hampel <- function (x, k) 
{
    n <- length(x)
    y <- x
    val <- rep(0, n)
    L <- 1.4826
    for (i in (k + 1):(n - k)) {
        x0 <- median(x[(i - k):(i + k)])
        S0 <- L * median(abs(x[(i - k):(i + k)] - x0))
        iden <- abs(x[i] - x0) / S0
        val[i] <- iden
    }
    val
}
