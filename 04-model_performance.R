#**************************************************************************
#* LFMC - A remote sensing approach to predict LFMC at large scale
#* Modelling approach:
#*    Model performance assessment and feature selection.
#*    Include Hybrid Forest Mask as input.
#* Author: Angel Cunill Camprubi <acunill86@gmail.com>
#* Last update: 27 Aug 2021
#**************************************************************************

# SETTINGS ----------------------------------------------------------------

# Packages
library(ranger)
library(Rcpp)
library(VSURF)
library(CAST)
library(doParallel)
library(foreach)

# Others
options(stringsAsFactors = FALSE)
RNGkind("L'Ecuyer-CMRG")

# Output directory
wdir <- "./data/tmp/mp"
if (!dir.exists(wdir)) dir.create(wdir)

# DATA PREPROCESSING ------------------------------------------------------

# Load database
db <- readRDS("./data/tmp/db03_lfmc.rds")

# Features collection
si <- c('NDVI', 'NDII6', 'NDII7', 'GVMI', 'NDWI', 'EVI', 'SAVI', 'VARI', 'VIgreen', 'NDTI', 'MSI', 'Gratio')
bands <- c('NR1', 'NR2', 'NR3', 'NR4', 'NR5', 'NR6', 'NR7')
vars <- c('HYBFOR', si, bands)
mdl_vars <- c("LFMC_AVG", vars)

# Define cross-validation folds
k <- 5
source("src/functions/create_llocv_folds.R")
set.seed(42)
db$FOLD <- create_llocv_folds(
    locations = c("LON", "LAT"),
    id = "SITE",
    data = db,
    kfold = k
)

# VARIABLE SELECTION ------------------------------------------------------
# Computational time ~24 hours (1000 ntree)

# Setting outputs
var_sel_cv <- list()

# 5-fold variable selection
t0 <- Sys.time()
for (val_fold in sort(unique(db$FOLD))) {
    cat(sprintf("\nValidation fold %i\n", val_fold))
    cat(sprintf("%s\n", as.character(Sys.time())))
    
    # Get data subset for variable selection
    train_df <- db[db$FOLD != val_fold, mdl_vars]
    
    # Remove noisy variables
    set.seed(42)
    var_sel <- VSURF(
        x = train_df[, vars],
        y = train_df[, "LFMC_AVG"],
        ntree = 500,
        parallel = TRUE
    )
    
    # Outputs
    var_sel_cv <- c(var_sel_cv, list(var_sel))
}
Sys.time() - t0

# Save results
#saveRDS(var_sel_cv, "./data/mp_varsel.rds")
remove(val_fold, var_sel, var_sel_cv)

# MODEL ASSESSMENT: ALLP-NO FILTERS ---------------------------------------
# Computational time ~36.5 (5 folds/50 hyp)

#* Hyperparameters ----
num_trees <- c(250, 500) #c(250, 500, 800, 1000)
mtry <- seq(1, ceiling(length(vars) * 0.3), 1)
min_node_size <- c(5:70)
sample_fraction <- c(0.632, 0.8, 1) #seq(1, .65, -.05)

# Seeds for reproducibility
set.seed(42)
seeds <- sample(1:1000, k)

# Setting outputs
cv_fold <- c()
cv_obs <- c()
cv_pred <- c()
out_params <- list()
var_importance <- list()

#* Model calibration and validation ----
t0 <- Sys.time()
for (val_fold in sort(unique(db$FOLD))) {
    cat(sprintf("\nValidation fold %i\n", val_fold))
    cat(sprintf("%s\n", as.character(Sys.time())))

    # Split data into train/test
    train_df <- db[db$FOLD != val_fold, ]
    test_df <- db[db$FOLD == val_fold, ]
    
    # Create indices for nested-cv folders
    indices <- CreateSpacetimeFolds(train_df, spacevar = "SITE", k = 5, seed = 42)
    
    # Random subset of hyperparameter combinations
    hyp <- expand.grid(
        "num_trees" = num_trees,
        "mtry" = mtry,
        "min_node_size" = min_node_size,
        "sample_fraction" = sample_fraction
    )
    set.seed(seeds[val_fold])
    hyp <- hyp[sample(nrow(hyp), 100), ]
    hyp_order <- order(
        hyp$mtry, hyp$num_trees, hyp$min_node_size,
        method = "radix",
        decreasing = c(FALSE, FALSE, TRUE)
    )
    hyp <- hyp[hyp_order, ]
    
    # Tune random forest (parallel)
    ncores <- parallel::detectCores() - 1
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    tmp <- foreach(h = seq_len(nrow(hyp))) %dopar% {
        library(ranger)
        
        hyp_i <- hyp[h, ]
        
        obs <- c()
        pred <- c()
        
        for (f in seq_along(indices$index)) {
            
            # Split training data into tune/assess
            tune_df <- train_df[indices$index[[f]], mdl_vars]
            asses_df <- train_df[indices$indexOut[[f]], mdl_vars]
            
            # Model tune
            mdl <- ranger(
                formula = LFMC_AVG ~ .,
                data = tune_df,
                num.trees = hyp_i[, "num_trees"],
                mtry = hyp_i[, "mtry"],
                min.node.size = hyp_i[, "min_node_size"],
                sample.fraction = hyp_i[, "sample_fraction"],
                replace = ifelse(hyp_i[, "sample_fraction"] == 1, TRUE, FALSE),
                oob.error = FALSE,
                splitrule = "variance",
                verbose = FALSE,
                seed = 42,
                num.threads = 1
            )
            
            # Inner predictions
            obs <- c(obs, asses_df$LFMC_AVG)
            pred <- c(pred, predict(mdl, asses_df, num.threads = 1)$predictions)
        }
        
        # Calculate 5-fold CV error rate
        hyp_i$rmse <- sqrt(mean((obs - pred) ^ 2, na.rm = TRUE))
        hyp_i
    }
    stopCluster(cl)
    hyp <- do.call(rbind, tmp)
    
    # Calibrate best model with the training data
    params <- hyp[which.min(hyp$rmse), ]
    t_mdl <- ranger(
        formula = LFMC_AVG ~ .,
        data = train_df[, mdl_vars],
        num.trees = params[, "num_trees"],
        mtry = params[, "mtry"],
        min.node.size = params[, "min_node_size"],
        sample.fraction = params[, "sample_fraction"],
        replace = ifelse(params[, "sample_fraction"] == 1, TRUE, FALSE),
        importance = "impurity",
        oob.error = FALSE,
        splitrule = "variance",
        verbose = FALSE,
        seed = 42,
        num.threads = ncores
    )
    
    # Outer predictions
    cv_fold <- c(cv_fold, test_df$FOLD)
    cv_obs <- c(cv_obs, test_df$LFMC_AVG)
    cv_pred <- c(cv_pred, predict(t_mdl, test_df[, mdl_vars], num.threads = ncores)$predictions)
    
    # Outputs
    out_params <- c(out_params, list(unlist(params)))
    var_importance <- c(var_importance, list(importance(t_mdl)))
}
Sys.time() - t0

#* Save results ----
nestedCV_output <- list(
    "validation" = data.frame("fold" = cv_fold, "obs" = cv_obs, "pred" = cv_pred),
    "parameters" = out_params,
    "variable_importance" = var_importance
)
saveRDS(nestedCV_output, "./data/mp_1_allp_nof.rds")

# MODEL ASSESSMENT SELECTED P-NO FILTERS ----------------------------------
# Computational time ~21.2 (5 folds)

# Clean memory
obj_list <- c("wdir", "db", "vars", "k", "si", "bands")
remove(list = ls()[!ls() %in% obj_list])

#* Hyperparameters ----
num_trees <- c(250, 500)
min_node_size <- c(5:70)
sample_fraction <- c(0.632, 0.8, 1)

#* Selected features ----
var_sel_cv <- readRDS("./data/mp_varsel.rds")

# Seeds for reproducibility
set.seed(42)
seeds <- sample(1:1000, k)

# Setting outputs
cv_fold <- c()
cv_obs <- c()
cv_pred <- c()
out_params <- list()
var_importance <- list()

#* Model calibration and validation ----
t0 <- Sys.time()
for (val_fold in sort(unique(db$FOLD))) {
    cat(sprintf("\nValidation fold %i\n", val_fold))
    cat(sprintf("%s\n", as.character(Sys.time())))
    
    # Subset of variables
    vars2 <- vars[var_sel_cv[[val_fold]]$varselect.pred]
    mdl_vars <- c("LFMC_AVG", vars2)
    
    # Split data into train/test
    train_df <- db[db$FOLD != val_fold, ]
    test_df <- db[db$FOLD == val_fold, ]
    
    # Create indices for nested-cv folders
    indices <- CreateSpacetimeFolds(train_df, spacevar = "SITE", k = 5, seed = 42)
    
    # Random subset of hyperparameter combinations
    if(ceiling(length(vars2) * 0.2) < 4) {
        mtry <- seq(1, 4, 1)
    } else {
        mtry <- seq(1, ceiling(length(vars2) * 0.2), 1)
    }
    hyp <- expand.grid(
        "num_trees" = num_trees,
        "mtry" = mtry,
        "min_node_size" = min_node_size,
        "sample_fraction" = sample_fraction
    )
    set.seed(seeds[val_fold])
    hyp <- hyp[sample(nrow(hyp), 100), ]
    hyp_order <- order(
        hyp$mtry, hyp$num_trees, hyp$min_node_size,
        method = "radix",
        decreasing = c(FALSE, FALSE, TRUE)
    )
    hyp <- hyp[hyp_order, ]
    
    # Tune random forest (parallel)
    ncores <- parallel::detectCores() - 1
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    tmp <- foreach(h = seq_len(nrow(hyp))) %dopar% {
        library(ranger)
        
        hyp_i <- hyp[h, ]
        
        obs <- c()
        pred <- c()
        
        for (f in seq_along(indices$index)) {
            
            # Split training data into tune/assess
            tune_df <- train_df[indices$index[[f]], mdl_vars]
            asses_df <- train_df[indices$indexOut[[f]], mdl_vars]
            
            # Model tune
            mdl <- ranger(
                formula = LFMC_AVG ~ .,
                data = tune_df,
                num.trees = hyp_i[, "num_trees"],
                mtry = hyp_i[, "mtry"],
                min.node.size = hyp_i[, "min_node_size"],
                sample.fraction = hyp_i[, "sample_fraction"],
                replace = ifelse(hyp_i[, "sample_fraction"] == 1, TRUE, FALSE),
                oob.error = FALSE,
                splitrule = "variance",
                verbose = FALSE,
                seed = 42,
                num.threads = 1
            )
            
            # Inner predictions
            obs <- c(obs, asses_df$LFMC_AVG)
            pred <- c(pred, predict(mdl, asses_df, num.threads = 1)$predictions)
        }
        
        # Calculate 5-fold CV error rate
        hyp_i$rmse <- sqrt(mean((obs - pred) ^ 2, na.rm = TRUE))
        hyp_i
    }
    stopCluster(cl)
    hyp <- do.call(rbind, tmp)
    
    # Calibrate best model with the training data
    params <- hyp[which.min(hyp$rmse), ]
    t_mdl <- ranger(
        formula = LFMC_AVG ~ .,
        data = train_df[, mdl_vars],
        num.trees = params[, "num_trees"],
        mtry = params[, "mtry"],
        min.node.size = params[, "min_node_size"],
        sample.fraction = params[, "sample_fraction"],
        replace = ifelse(params[, "sample_fraction"] == 1, TRUE, FALSE),
        importance = "impurity",
        oob.error = FALSE,
        splitrule = "variance",
        verbose = FALSE,
        seed = 42,
        num.threads = ncores
    )
    
    # Outer predictions
    cv_fold <- c(cv_fold, test_df$FOLD)
    cv_obs <- c(cv_obs, test_df$LFMC_AVG)
    cv_pred <- c(cv_pred, predict(t_mdl, test_df[, mdl_vars], num.threads = ncores)$predictions)
    
    # Outputs
    out_params <- c(out_params, list(unlist(params)))
    var_importance <- c(var_importance, list(importance(t_mdl)))
}
Sys.time() - t0

#* Save results ----
nestedCV_output <- list(
    "validation" = data.frame("fold" = cv_fold, "obs" = cv_obs, "pred" = cv_pred),
    "parameters" = out_params,
    "variable_importance" = var_importance
)
saveRDS(nestedCV_output, "./data/mp_2_selp_nof.rds")

# FILTERING SETTINGS ------------------------------------------------------

# Remove missing data
db <- db[complete.cases(db$NDVIcv), ]
mdl_vars <- c("LFMC_AVG", vars)

# Re-define cross-validation folds
source("src/functions/create_llocv_folds.R")
set.seed(42)
db$FOLD <- create_llocv_folds(
    locations = c("LON", "LAT"),
    id = "SITE",
    data = db,
    kfold = k
)

# MODEL ASSESSMENT: ALLP-FILTERS ------------------------------------------
# Computational time ~16.4 (5 folds)

# Clean memory
obj_list <- c("wdir", "db", "mdl_vars", "vars", "k", "si", "bands")
remove(list = ls()[!ls() %in% obj_list])

#* Hyperparameters ----
sp_filter <- seq(0.5, 0.1, -0.05)
num_trees <- c(250, 500)
min_node_size <- c(5:70)
sample_fraction <- c(0.632, 0.8, 1)

# Seeds for reproducibility
set.seed(42)
seeds <- sample(1:1000, k)

# Setting outputs
cv_fold <- c()
cv_obs <- c()
cv_pred <- c()
out_params <- list()
var_importance <- list()

#* Model calibration and validation ----
t0 <- Sys.time()
for (val_fold in sort(unique(db$FOLD))) {
    cat(sprintf("\nValidation fold %i\n", val_fold))
    cat(sprintf("%s\n", as.character(Sys.time())))
    
    # Split data into train/test
    train_df <- db[db$FOLD != val_fold, ]
    test_df <- db[db$FOLD == val_fold, ]
    
    # Random subset of hyperparameter combinations
    if(ceiling(length(vars) * 0.2) < 4) {
        mtry <- seq(1, 4, 1)
    } else {
        mtry <- seq(1, ceiling(length(vars) * 0.2), 1)
    }
    hyp <- expand.grid(
        "spf" = sp_filter,
        "num_trees" = num_trees,
        "mtry" = mtry,
        "min_node_size" = min_node_size,
        "sample_fraction" = sample_fraction
    )
    set.seed(seeds[val_fold])
    hyp <- hyp[sample(nrow(hyp), 150), ]
    hyp_order <- order(
        hyp$mtry, hyp$num_trees, hyp$min_node_size,
        method = "radix",
        decreasing = c(FALSE, FALSE, TRUE)
    )
    hyp <- hyp[hyp_order, ]
    
    # Tune random forest (parallel)
    ncores <- parallel::detectCores() - 1
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    tmp <- foreach(h = seq_len(nrow(hyp))) %dopar% {
        library(ranger)
        library(CAST)
        
        hyp_i <- hyp[h, ]
        
        obs <- c()
        pred <- c()
        
        # Apply filter to training data
        f_train_df <- train_df[train_df$NDVIcv < hyp_i$spf, ]
        
        # Create indices for nested-cv folders
        indices <- CreateSpacetimeFolds(f_train_df, spacevar = "SITE", k = 5, seed = 42)
        
        # Train model
        for (f in seq_along(indices$index)) {
            
            # Split training data into tune/assess
            tune_df <- f_train_df[indices$index[[f]], mdl_vars]
            asses_df <- f_train_df[indices$indexOut[[f]], mdl_vars]
            
            # Model tune
            mdl <- ranger(
                formula = LFMC_AVG ~ .,
                data = tune_df,
                num.trees = hyp_i[, "num_trees"],
                mtry = hyp_i[, "mtry"],
                min.node.size = hyp_i[, "min_node_size"],
                sample.fraction = hyp_i[, "sample_fraction"],
                replace = ifelse(hyp_i[, "sample_fraction"] == 1, TRUE, FALSE),
                oob.error = FALSE,
                splitrule = "variance",
                verbose = FALSE,
                seed = 42,
                num.threads = 1
            )
            
            # Inner predictions
            obs <- c(obs, asses_df$LFMC_AVG)
            pred <- c(pred, predict(mdl, asses_df, num.threads = 1)$predictions)
        }
        
        # Calculate 5-fold CV error rate
        hyp_i$rmse <- sqrt(mean((obs - pred) ^ 2, na.rm = TRUE))
        hyp_i
    }
    stopCluster(cl)
    hyp <- do.call(rbind, tmp)
    
    # Calibrate best model with the training data
    params <- hyp[which.min(hyp$rmse), ]
    t_mdl <- ranger(
        formula = LFMC_AVG ~ .,
        data = train_df[train_df$NDVIcv < params$spf, mdl_vars],
        num.trees = params[, "num_trees"],
        mtry = params[, "mtry"],
        min.node.size = params[, "min_node_size"],
        sample.fraction = params[, "sample_fraction"],
        replace = ifelse(params[, "sample_fraction"] == 1, TRUE, FALSE),
        importance = "impurity",
        oob.error = FALSE,
        splitrule = "variance",
        verbose = FALSE,
        seed = 42,
        num.threads = ncores
    )
    
    # Outer predictions
    cv_fold <- c(cv_fold, test_df$FOLD)
    cv_obs <- c(cv_obs, test_df$LFMC_AVG)
    cv_pred <- c(cv_pred, predict(t_mdl, test_df[, mdl_vars], num.threads = ncores)$predictions)
    
    # Outputs
    out_params <- c(out_params, list(unlist(params)))
    var_importance <- c(var_importance, list(importance(t_mdl)))
}
Sys.time() - t0

#* Save results ----
nestedCV_output <- list(
    "validation" = data.frame("fold" = cv_fold, "obs" = cv_obs, "pred" = cv_pred),
    "parameters" = out_params,
    "variable_importance" = var_importance
)
saveRDS(nestedCV_output, paste0("./data/mp_3_allp_filt.rds"))

# MODEL ASSESSMENT: SELECTED P-FILTERS ------------------------------------
# Computational time ~16.4 (5 folds)

# Clean memory
obj_list <- c("wdir", "db", "vars", "k", "si", "bands")
remove(list = ls()[!ls() %in% obj_list])

#* Hyperparameters ----
sp_filter <- seq(0.5, 0.1, -0.05)
num_trees <- c(250, 500)
min_node_size <- c(5:70)
sample_fraction <- c(0.632, 0.8, 1)

#* Selected features ----
var_sel_cv <- readRDS("./data/mp_varsel.rds")

# Seeds for reproducibility
set.seed(42)
seeds <- sample(1:1000, k)

# Setting outputs
cv_fold <- c()
cv_obs <- c()
cv_pred <- c()
out_params <- list()
var_importance <- list()

#* Model calibration and validation ----
t0 <- Sys.time()
for (val_fold in sort(unique(db$FOLD))) {
    cat(sprintf("\nValidation fold %i\n", val_fold))
    cat(sprintf("%s\n", as.character(Sys.time())))
    
    # Subset of variables
    vars2 <- vars[var_sel_cv[[val_fold]]$varselect.pred]
    mdl_vars <- c("LFMC_AVG", vars2)
    
    # Split data into train/test
    train_df <- db[db$FOLD != val_fold, ]
    test_df <- db[db$FOLD == val_fold, ]
    
    # Random subset of hyperparameter combinations
    if(ceiling(length(vars2) * 0.2) < 4) {
        mtry <- seq(1, 4, 1)
    } else {
        mtry <- seq(1, ceiling(length(vars2) * 0.2), 1)
    }
    hyp <- expand.grid(
        "spf" = sp_filter,
        "num_trees" = num_trees,
        "mtry" = mtry,
        "min_node_size" = min_node_size,
        "sample_fraction" = sample_fraction
    )
    set.seed(seeds[val_fold])
    hyp <- hyp[sample(nrow(hyp), 150), ]
    hyp_order <- order(
        hyp$mtry, hyp$num_trees, hyp$min_node_size,
        method = "radix",
        decreasing = c(FALSE, FALSE, TRUE)
    )
    hyp <- hyp[hyp_order, ]
    
    # Tune random forest (parallel)
    ncores <- parallel::detectCores() - 1
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    tmp <- foreach(h = seq_len(nrow(hyp))) %dopar% {
        library(ranger)
        library(CAST)
        
        hyp_i <- hyp[h, ]
        
        obs <- c()
        pred <- c()
        
        # Apply filter to training data
        f_train_df <- train_df[train_df$NDVIcv < hyp_i$spf, ]
        
        # Create indices for nested-cv folders
        indices <- CreateSpacetimeFolds(f_train_df, spacevar = "SITE", k = 5, seed = 42)
        
        # Train model
        for (f in seq_along(indices$index)) {
            
            # Split training data into tune/assess
            tune_df <- f_train_df[indices$index[[f]], mdl_vars]
            asses_df <- f_train_df[indices$indexOut[[f]], mdl_vars]
            
            # Model tune
            mdl <- ranger(
                formula = LFMC_AVG ~ .,
                data = tune_df,
                num.trees = hyp_i[, "num_trees"],
                mtry = hyp_i[, "mtry"],
                min.node.size = hyp_i[, "min_node_size"],
                sample.fraction = hyp_i[, "sample_fraction"],
                replace = ifelse(hyp_i[, "sample_fraction"] == 1, TRUE, FALSE),
                oob.error = FALSE,
                splitrule = "variance",
                verbose = FALSE,
                seed = 42,
                num.threads = 1
            )
            
            # Inner predictions
            obs <- c(obs, asses_df$LFMC_AVG)
            pred <- c(pred, predict(mdl, asses_df, num.threads = 1)$predictions)
        }
        
        # Calculate 5-fold CV error rate
        hyp_i$rmse <- sqrt(mean((obs - pred) ^ 2, na.rm = TRUE))
        hyp_i
    }
    stopCluster(cl)
    hyp <- do.call(rbind, tmp)
    
    # Calibrate best model with the training data
    params <- hyp[which.min(hyp$rmse), ]
    t_mdl <- ranger(
        formula = LFMC_AVG ~ .,
        data = train_df[train_df$NDVIcv < params$spf, mdl_vars],
        num.trees = params[, "num_trees"],
        mtry = params[, "mtry"],
        min.node.size = params[, "min_node_size"],
        sample.fraction = params[, "sample_fraction"],
        replace = ifelse(params[, "sample_fraction"] == 1, TRUE, FALSE),
        importance = "impurity",
        oob.error = FALSE,
        splitrule = "variance",
        verbose = FALSE,
        seed = 42,
        num.threads = ncores
    )
    
    # Outer predictions
    cv_fold <- c(cv_fold, test_df$FOLD)
    cv_obs <- c(cv_obs, test_df$LFMC_AVG)
    cv_pred <- c(cv_pred, predict(t_mdl, test_df[, mdl_vars], num.threads = ncores)$predictions)
    
    # Outputs
    out_params <- c(out_params, list(unlist(params)))
    var_importance <- c(var_importance, list(importance(t_mdl)))
}
Sys.time() - t0

#* Save results ----
nestedCV_output <- list(
    "validation" = data.frame("fold" = cv_fold, "obs" = cv_obs, "pred" = cv_pred),
    "parameters" = out_params,
    "variable_importance" = var_importance
)
saveRDS(nestedCV_output, paste0("./data/mp_4_selp_filt.rds"))

# END ---
