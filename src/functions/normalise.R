# *************************************************************************
# LFMC-RF - A remote sensing approach to predict LFMC at large scale using
# Random Forests
# Normalise data function
# Adapted from https://github.com/lynn-miller/LFMC_from_MODIS
#  - x : vector data to normalise
#  - in_range : range of source data; calculated from the source data if not specified
#  - out_range : range of normalised data; range is 0-1 if not specified
# Author: Angel Cunill Camprubi <acunill86@gmail.com>
# Last update: 20 Oct 2021
# *************************************************************************

normalise <- function(x, in_range = NULL, out_range = c(0, 1)) {
    if (is.null(in_range)) {
        in_range <- c(min(x), max(x))
    }
    tmp <- (x - in_range[1]) / (in_range[2] - in_range[1])
    if (identical(out_range, c(0, 1))) {
        tmp
    } else {
        tmp * (out_range[2] - out_range[1]) + out_range[1]
    }
}

# END ---
