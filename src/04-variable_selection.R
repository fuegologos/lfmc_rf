# *************************************************************************
# LFMC-RF - A remote sensing approach to predict LFMC at large scale using
# Random Forests
# Variable selection using Forward Feature Selection (FFS) as described in
# Meyer et al. (2018)
# Author: Angel Cunill Camprubi <acunill86@gmail.com>
# Last update: 03 Mar 2022
# *************************************************************************

# SETTINGS ----------------------------------------------------------------

# Packages
source("src/functions/ffs.R")

# Others
options(stringsAsFactors = FALSE)
RNGkind("L'Ecuyer-CMRG")
    
# Working directory
dir_out <- "./data/vs"
dir.create(dir_out)

# INITIAL VARIABLE COLLECTION APPROACHES ----------------------------------

# Predictors
si <- c('NDVI', 'EVI', 'SAVI', 'VARI', 'VIgreen', 'Gratio', 'NDII6', 'NDII7',
        'NDWI', 'GVMI', 'MSI', 'NDTI', 'STI')
bands <- c('NR1', 'NR2', 'NR3', 'NR4', 'NR5', 'NR6', 'NR7')

vars_list <- list(
    "appr_1" = c(si, bands, "LST", "DOY_COS", "DOY_SIN"),
    "appr_2" = c(si, bands, "LST8day", "DOY_COS", "DOY_SIN"),
    "appr_3" = c(si, "LST", "DOY_COS", "DOY_SIN"),
    "appr_4" = c(si, "LST8day", "DOY_COS", "DOY_SIN")
)

# COMPARING APPROACHES ----------------------------------------------------
# Comparison with the same number of cal/test samples
# ~10 hours

# Load database and filtering
db <- readRDS("./data/db03_lfmc.rds")
vars <- unique(do.call(c, vars_list))
db <- db[complete.cases(db[, vars]), ]

# Variable selection over the list of approaches
for (appr in 1:4) {
    cat(sprintf("\nApproach %i\n", appr))
    
    # Model variables
    vars <- vars_list[[appr]]
    mdl_vars <- c("LFMC",  vars)
    
    # Temporary directory
    tmpdir <- file.path(dir_out, paste0("appr", appr))
    dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
    
    # Variable selection
    varsel <- ffs_iter(
        predictors = vars,
        response = "LFMC",
        df = db,
        wdir = tmpdir,
        ncores = parallel::detectCores() - 1
    )
    saveRDS(varsel, file.path(dir_out, paste0("appr", appr, ".rds")))
}

# VARIABLE SELECTION ------------------------------------------------------
# Repeat the process with the whole available data for the best approach

# Variable selection over the best approach
appr <- 2

# Model variables
vars <- vars_list[[appr]]
mdl_vars <- c("LFMC",  vars)

# Load database and filtering
db <- readRDS("./data/db03_lfmc.rds")
db <- db[db$LFMC >= 20 & db$LFMC <= 250, ]
db <- db[complete.cases(db[, mdl_vars]), ]

# Temporary directory
final_dir <- file.path(dir_out, "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

# Variable selection
varsel <- ffs_iter(
    predictors = vars,
    response = "LFMC",
    df = db,
    wdir = final_dir,
    ncores = 10
)
saveRDS(varsel, file.path(dir_out, "varsel_final.rds"))

# END ---
