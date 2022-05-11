# *************************************************************************
# LFMC-RF - A remote sensing approach to predict LFMC at large scale using
# Random Forests
# Calibration and validation of the LFMC-RF model
# NDVIcv filter applied only to the training data.
# Author: Angel Cunill Camprubi <acunill86@gmail.com>
# Last update: 13 Mar 2022
# *************************************************************************

# SETTINGS ----------------------------------------------------------------

# Packages
library(ranger)
library(Rcpp)
library(CAST)
library(MLmetrics)

# Cross-validation folds
source("src/functions/create_llocv_folds.R")

# Cross-validation settings
niter <- 25              #number of iterations
k <- 5                   #cross-validation folds
filt_opt <- "nofilt"     #(nofilt, filt)

# Others
options(stringsAsFactors = FALSE)
ncores <- parallel::detectCores() - 1
RNGkind("L'Ecuyer-CMRG")

# Create working directory
wdir <- "data/calext"
dir.create(wdir, showWarnings = FALSE)

# MODEL INPUTS ------------------------------------------------------------

# Predictors
pred <- readRDS("data/vs/varsel_final.rds")$predictors
pred <- pred[[length(pred)]]
mdl_vars <- c("LFMC", pred)

# Filters
sp_filt <- seq(0.4, 0.2, -0.05)
no_filt <- 1000

# Hyperparameters
num_trees <- c(500)
min_node_size <- c(5:10)
sample_fraction <- seq(0.2, 0.7, 0.05)
mtry <- 2:3

# ITERATIVE MODEL TRAINING ------------------------------------------------
# 30 mins

# Create temporary directory
tmpdir <- file.path(wdir, filt_opt)
dir.create(tmpdir)

# Load database
db <- readRDS("data/db03_lfmc.rds")
db <- db[complete.cases(db[, mdl_vars]), ]

# Split data into development and validation
dev_df <- db[db$YEAR < 2015, ]
val_df <- db[db$YEAR >= 2015, ]

# Filter settings
if (filt_opt == "filt") {
    dev_df <- dev_df[complete.cases(dev_df$NDVIcv), ]
    spfilt <- sp_filt
} else {
    spfilt <- no_filt
}

# Hyperparameter grid
hyp <- expand.grid(
    "spfilt" = spfilt,
    "num_trees" = num_trees,
    "mtry" = mtry,
    "min_node_size" = min_node_size,
    "sample_fraction" = sample_fraction
)

# Seeds for iterative cross-validation folds
set.seed(42)
seeds <- sample(1:9999, niter)

# Model trainning
for (iter_i in 1:niter) {
    
    cat(sprintf("\nIteration %i of %i: %s\n",
                iter_i,
                niter,
                as.character(Sys.time())))
    
    # Define outer cross-validation folds
    indices <- CreateSpacetimeFolds(dev_df,
                                    spacevar = "SITE",
                                    k = k,
                                    seed = seeds[iter_i])
    
    # Train and test random forest
    hyp_error <- c()
    for (h in seq_len(nrow(hyp))) {
        
        # Settings
        hyp_i <- hyp[h, ]
        obs <- c()
        pred <- c()
        
        # Cross-validation
        for (fold in seq_along(indices$index)) {
            
            # Split data into train/test
            train_df <- dev_df[indices$index[[fold]], ]
            test_df <- dev_df[indices$indexOut[[fold]], mdl_vars]
            
            # Apply filter to training data
            if (filt_opt == "filt") {
                #train_df <- train_df[complete.cases(train_df$NDVIcv), ]
                train_df <- train_df[train_df$NDVIcv < hyp_i$spfilt, ]    
            }
            
            # Model tune
            mdl <- ranger(
                formula = LFMC ~ .,
                data = train_df[, mdl_vars],
                num.trees = hyp_i$num_trees,
                mtry = hyp_i$mtry,
                min.node.size = hyp_i$min_node_size,
                sample.fraction = hyp_i$sample_fraction,
                replace = FALSE,
                oob.error = FALSE,
                splitrule = "variance",
                verbose = FALSE,
                seed = 42,
                num.threads = ncores
            )
            
            # Predictions
            obs <- c(obs, test_df$LFMC)
            pred <- c(pred, predict(mdl, test_df, num.threads = ncores)$predictions)
        }
        
        # Calculate CV error
        hyp_error <- c(hyp_error,
                       RMSE(y_pred = pred, y_true = obs))
    }
    
    # Save current iteration results
    saveRDS(hyp_error, file.path(tmpdir, paste0("iter_", iter_i, ".rds")))
}

# Save error profile
iter_files <- paste0("iter_", 1:niter, ".rds")
tmp <- lapply(iter_files, function(x) {
    readRDS(file.path(tmpdir, x))
})
iter_err <- as.data.frame(do.call(cbind, tmp))
colnames(iter_err) <- paste0("iter", 1:niter)
iter_err$rmse <- apply(iter_err, 1, mean)
hyp_prof <- cbind(hyp, iter_err)
saveRDS(hyp_prof, file.path(tmpdir, "hyp_eprof.rds"))
remove(tmp)

# MODEL VALIDATION: ADJUSTED MODEL ----------------------------------------

# Best hyperparameters
hyp_prof <- readRDS(file.path(tmpdir, "hyp_eprof.rds"))
best_hyp <- hyp_prof[which.min(hyp_prof$rmse), c(1:5)]

# Define cross-validation folds
set.seed(42)
dev_df$FOLD <- create_llocv_folds(locations = c("LON", "LAT"),
                                  id = "SITE",
                                  data = dev_df,
                                  kfold = k)

# Setting outputs
site_id <- c()
cv_fold <- c()
cv_obs <- c()
cv_pred <- c()

# Best model cross-validation
for (val_fold in sort(unique(dev_df$FOLD))) {
    
    # Split data into train/test
    train_df <- dev_df[dev_df$FOLD != val_fold, ]
    test_df <- dev_df[dev_df$FOLD == val_fold, ]
    
    # Apply filter to training data
    if(filt_opt == "filt") {
        train_df <- train_df[train_df$NDVIcv < best_hyp$spfilt, ]   
    }
    
    # Calibrate model
    best_mdl <- ranger(
        formula = LFMC ~ .,
        data = train_df[, mdl_vars],
        num.trees = best_hyp$num_trees,
        mtry = best_hyp$mtry,
        min.node.size = best_hyp$min_node_size,
        sample.fraction = best_hyp$sample_fraction,
        replace = FALSE,
        oob.error = FALSE,
        splitrule = "variance",
        verbose = FALSE,
        seed = 42,
        num.threads = ncores
    )
    
    # Predictions
    site_id <- c(site_id, test_df$ID)
    cv_fold <- c(cv_fold, test_df$FOLD)
    cv_obs <- c(cv_obs, test_df$LFMC)
    cv_pred <- c(cv_pred,
                 predict(best_mdl, test_df[, mdl_vars], num.threads = ncores)$predictions)
}

# Save results
dev_out <- data.frame(
    "ID" = site_id,
    "fold" = cv_fold,
    "obs" = cv_obs,
    "pred" = cv_pred
)

# MODEL VALIDATION: EXTRAPOLATION -----------------------------------------

# Apply filter to development data
if(filt_opt == "filt") {
    dev_df <- dev_df[dev_df$NDVIcv < best_hyp$spfilt, ]
}

# Calibrate best model
best_mdl <- ranger(
    formula = LFMC ~ .,
    data = dev_df[, mdl_vars],
    num.trees = best_hyp$num_trees,
    mtry = best_hyp$mtry,
    min.node.size = best_hyp$min_node_size,
    sample.fraction = best_hyp$sample_fraction,
    replace = FALSE,
    oob.error = FALSE,
    splitrule = "variance",
    verbose = FALSE,
    seed = 42,
    num.threads = ncores
)

# Save results
output <- list(
    "dev" = dev_out,
    "val" = data.frame(
        "ID" = val_df$ID,
        "obs" = val_df$LFMC,
        "pred" = predict(best_mdl, val_df[, mdl_vars], num.threads = ncores)$predictions
    ),
    "hyp" = best_hyp,
    "sum" = data.frame(
        "data" = c("dev", "val"),
        "n" = c(nrow(dev_df), nrow(val_df)),
        "min" = c(min(dev_df$LFMC), min(val_df$LFMC)),
        "max" = c(max(dev_df$LFMC), max(val_df$LFMC)),
        "sd" = c(sd(dev_df$LFMC), sd(val_df$LFMC)),
        "mean" = c(mean(dev_df$LFMC), mean(val_df$LFMC)),
        "median" = c(median(dev_df$LFMC), median(val_df$LFMC))
    )
)
saveRDS(output, file.path(wdir, paste0(filt_opt, ".rds")))

# END ---
