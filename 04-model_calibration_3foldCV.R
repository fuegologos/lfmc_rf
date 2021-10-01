#**************************************************************************
#* LFMC - A remote sensing approach to predict LFMC at large scale
#* Modelling approach:
#*    Model performance assessment using 3-fold LLOCV.
#* Author: Angel Cunill Camprubi <acunill86@gmail.com>
#* Last update: 06 Sept 2021
#**************************************************************************

# SETTINGS ----------------------------------------------------------------

# Packages
library(ranger)
library(Rcpp)
library(CAST)
library(VSURF)
library(doParallel)
library(foreach)

# Others
options(stringsAsFactors = FALSE)
ncores <- parallel::detectCores() - 1
RNGkind("L'Ecuyer-CMRG")

# Output directory
wdir <- "./data/tmp/mp"
if (!dir.exists(wdir)) dir.create(wdir)

# DATA PREPROCESSING ------------------------------------------------------

# Load database
db <- readRDS("./data/tmp/db03_lfmc.rds")

# Initial features
si <- c('NDVI', 'NDII6', 'NDII7', 'GVMI', 'NDWI', 'EVI', 'SAVI', 'VARI', 'VIgreen', 'NDTI', 'MSI', 'Gratio')
bands <- c('NR1', 'NR2', 'NR3', 'NR4', 'NR5', 'NR6', 'NR7')
vars <- c('HYBFOR', si, bands)

# VARIABLE SELECTION ------------------------------------------------------
# Computational time ~7.2 hours (1000 ntree)

# Variable selection
t0 <- Sys.time()
set.seed(42)
var_sel <- VSURF(
    x = db[, vars],
    y = db[, "LFMC_AVG"],
    ntree = 1000,
    parallel = TRUE
)
Sys.time() - t0

# Save results
#saveRDS(var_sel, "./data/3fcv_varsel.rds")

# Feature options
vs <- readRDS("./data/3fcv_varsel.rds")
vars_df <- list(
    "name" = c("allp", "selp"),
    "variables" = list(vars, vars[vs$varselect.pred])
)
remove(var_sel, vs)

# HYPERPARAMETERS ---------------------------------------------------------

num_trees <- 250
mtry <- 1:2
min_node_size <- 15:75
sample_fraction <- c(0.632, 1)

# 3-FOLD CROSS VALIDATION: NO FILTERS -------------------------------------
# Computational time: ~3 hours

#* Hyperparameters ----
hyp <- expand.grid(
    "num_trees" = num_trees,
    "mtry" = mtry,
    "min_node_size" = min_node_size,
    "sample_fraction" = sample_fraction
)

# Train and test
t0 <- Sys.time()
for (vopt in seq_along(vars_df$name)) {
    cat(sprintf("\nNo filters - %s", vars_df$name[vopt]))
    
    #* Select model vars ----
    vars2 <- vars_df$variables[[vopt]]
    mdl_vars <- c("LFMC_AVG", vars2)
    
    #* Train model ----
    # Seeds for reproducibility
    set.seed(42)
    seeds <- sample(1:1000, 25)
    
    # Iterative cross-validation process
    for (i in 1:25) {
        cat(sprintf("\nIteration %i of 25: %s\n", i, as.character(Sys.time())))
        
        # Create cross-validation folds
        indices <- CreateSpacetimeFolds(db, spacevar = "SITE", k = 3, seed = seeds[i])
        
        # Train and test random forest
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        tmp <- foreach(h = seq_len(nrow(hyp))) %dopar% {
            library(ranger)
            
            hyp_i <- hyp[h, ]
            
            obs <- c()
            pred <- c()
            
            # Cross-validation
            for (fold in seq_along(indices$index)) {
                
                # Split data into train/test
                train_df <- db[indices$index[[fold]], mdl_vars]
                test_df <- db[indices$indexOut[[fold]], mdl_vars]
                
                # Model tune
                mdl <- ranger(
                    formula = LFMC_AVG ~ .,
                    data = train_df,
                    num.trees = hyp_i[, "num_trees"],
                    mtry = hyp_i[, "mtry"],
                    min.node.size = hyp_i[, "min_node_size"],
                    sample.fraction = hyp_i[, "sample_fraction"],
                    replace = FALSE,
                    oob.error = FALSE,
                    splitrule = "variance",
                    verbose = FALSE,
                    seed = 42,
                    num.threads = 1
                )
                
                # Predictions
                obs <- c(obs, test_df$LFMC_AVG)
                pred <- c(pred, predict(mdl, test_df, num.threads = 1)$predictions)
            }
            
            # Calculate CV error
            sqrt(mean((obs - pred) ^ 2, na.rm = TRUE))
        }
        stopCluster(cl)
        
        # Save partial outputs
        cv_error <- do.call(rbind, tmp)
        saveRDS(cv_error, file.path(wdir, paste0("cve_", i, ".rds")))
    }
    
    #* Join partial outputs ----
    files_list <- paste0("cve_", 1:25, ".rds")
    tmp <- lapply(files_list, function(x) {
        readRDS(file.path(wdir, x))
    })
    errors <- as.data.frame(do.call(cbind, tmp))
    colnames(errors) <- paste0("Iter", 1:25)
    df <- cbind(hyp, errors)
    
    #* Compute mean RMSE and save error profile ----
    df$rmse <- apply(df[, (1:25) + 4], 1, mean)
    saveRDS(df, paste0("./data/3fcv_eprof_", vars_df$name[vopt], "_nof.rds"))
    
    #* Calibrate final model ----
    # Define cross-validation folds
    source("src/functions/create_llocv_folds.R")
    set.seed(42)
    db$FOLD <- create_llocv_folds(
        locations = c("LON", "LAT"),
        id = "SITE",
        data = db,
        kfold = 3
    )
    
    # Parameters
    params <- df[which.min(df$rmse), c(1:4)]
    
    # Setting outputs
    site_id <- c()
    cv_fold <- c()
    cv_obs <- c()
    cv_pred <- c()
    var_importance <- list()
    
    # Best model calibration
    for (val_fold in sort(unique(db$FOLD))) {
        cat(sprintf("\nValidation fold %i\n", val_fold))
        cat(sprintf("%s\n", as.character(Sys.time())))
        
        # Split data into train/test
        train_df <- db[db$FOLD != val_fold, mdl_vars]
        test_df <- db[db$FOLD == val_fold, ]
        
        # Calibrate model
        t_mdl <- ranger(
            formula = LFMC_AVG ~ .,
            data = train_df,
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
        
        # Predictions and variable importance
        site_id <- c(site_id, test_df$ID)
        cv_fold <- c(cv_fold, test_df$FOLD)
        cv_obs <- c(cv_obs, test_df$LFMC_AVG)
        cv_pred <- c(cv_pred, predict(t_mdl, test_df[, mdl_vars], num.threads = ncores)$predictions)
        var_importance <- c(var_importance, list(importance(t_mdl)))
    }
    
    #* Save results ----
    cv_output <- list(
        "validation" = data.frame("ID" = site_id,
                                  "fold" = cv_fold,
                                  "obs" = cv_obs,
                                  "RFcv" = cv_pred),
        "variable_importance" = var_importance
    )
    file_name <- paste("./data/3fcv", vopt, vars_df$name[vopt], "nof.rds", sep = "_")
    saveRDS(cv_output, file_name)
}
Sys.time() - t0

# 3-FOLD CROSS VALIDATION: FILTERS ----------------------------------------
# Computational time: 8.832517 hours

# Clean memory
obj_list <- c("wdir", "db", "vars_df", "ncores",
              "num_trees", "mtry", "min_node_size", "sample_fraction")
remove(list = ls()[!ls() %in% obj_list])

#* Filters settings ----
db <- db[complete.cases(db$NDVIcv), ]
sp_filter <- seq(0.5, 0.3, -0.05)
ts_filter <- max(db$Q3) + 1 #seq(1, 0.2, -0.1)

#* Hyperparameters ----
hyp <- expand.grid(
    "sp" = sp_filter,
    "ts" = ts_filter,
    "num_trees" = num_trees,
    "mtry" = mtry,
    "min_node_size" = min_node_size,
    "sample_fraction" = sample_fraction
)

# Train and test with all p and selected p'
t0 <- Sys.time()
for (vopt in seq_along(vars_df$name)) {
    cat(sprintf("\nFilters - %s", vars_df$name[vopt]))
    
    #* Select model vars ----
    vars2 <- vars_df$variables[[vopt]]
    mdl_vars <- c("LFMC_AVG", vars2)
    
    #* Train model ----
    # Seeds for reproducibility
    set.seed(42)
    seeds <- sample(1:1000, 25)
    
    # Iterative cross-validation process
    for (i in 1:25) {
        cat(sprintf("\nIteration %i of 25: %s\n", i, as.character(Sys.time())))
        
        # Create cross-validation folds
        indices <- CreateSpacetimeFolds(db, spacevar = "SITE", k = 3, seed = seeds[i])
        
        # Train and test random forest
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        tmp <- foreach(h = seq_len(nrow(hyp))) %dopar% {
            library(ranger)
            
            hyp_i <- hyp[h, ]
            
            obs <- c()
            pred <- c()
            
            # Cross-validation
            for (fold in seq_along(indices$index)) {
                
                # Split data into train/test
                train_df <- db[indices$index[[fold]], ]
                test_df <- db[indices$indexOut[[fold]], mdl_vars]
                
                # Apply filters combination to training data
                train_df <- train_df[train_df$NDVIcv < hyp_i[, "sp"], ]
                train_df <- train_df[train_df$Q3 < hyp_i[, "ts"], mdl_vars]
                
                # Model tune
                mdl <- ranger(
                    formula = LFMC_AVG ~ .,
                    data = train_df,
                    num.trees = hyp_i[, "num_trees"],
                    mtry = hyp_i[, "mtry"],
                    min.node.size = hyp_i[, "min_node_size"],
                    sample.fraction = hyp_i[, "sample_fraction"],
                    replace = FALSE,
                    oob.error = FALSE,
                    splitrule = "variance",
                    verbose = FALSE,
                    seed = 42,
                    num.threads = 1
                )
                
                # Predictions
                obs <- c(obs, test_df$LFMC_AVG)
                pred <- c(pred, predict(mdl, test_df, num.threads = 1)$predictions)
            }
            
            # Calculate CV error
            sqrt(mean((obs - pred) ^ 2, na.rm = TRUE))
        }
        stopCluster(cl)
        
        # Save partial outputs
        cv_error <- do.call(rbind, tmp)
        saveRDS(cv_error, file.path(wdir, paste0("cve_", i, ".rds")))
    }
    
    #* Join partial outputs ----
    files_list <- paste0("cve_", 1:25, ".rds")
    tmp <- lapply(files_list, function(x) {
        readRDS(file.path(wdir, x))
    })
    errors <- as.data.frame(do.call(cbind, tmp))
    colnames(errors) <- paste0("Iter", 1:25)
    df <- cbind(hyp, errors)
    
    #* Compute mean RMSE and save error profile ----
    df$rmse <- apply(df[, (1:25) + 6], 1, mean)
    saveRDS(df, paste0("./data/3fcv_eprof_", vars_df$name[vopt], "_filt.rds"))
    
    #* Calibrate final model ----
    # Define cross-validation folds
    source("src/functions/create_llocv_folds.R")
    set.seed(42)
    db$FOLD <- create_llocv_folds(
        locations = c("LON", "LAT"),
        id = "SITE",
        data = db,
        kfold = 3
    )
    
    # Parameters
    params <- df[which.min(df$rmse), c(1:6)]
    
    # Setting outputs
    site_id <- c()
    cv_fold <- c()
    cv_obs <- c()
    cv_pred <- c()
    var_importance <- list()
    
    # Best model calibration
    for (val_fold in sort(unique(db$FOLD))) {
        cat(sprintf("\nValidation fold %i\n", val_fold))
        cat(sprintf("%s\n", as.character(Sys.time())))
        
        # Split data into train/test
        train_df <- db[db$FOLD != val_fold, ]
        test_df <- db[db$FOLD == val_fold, ]
        
        # Apply filters combination to training data
        train_df <- train_df[train_df$NDVIcv < params[, "sp"], ]
        train_df <- train_df[train_df$Q3 < params[, "ts"], mdl_vars]
        
        # Calibrate model
        t_mdl <- ranger(
            formula = LFMC_AVG ~ .,
            data = train_df,
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
        
        # Predictions and variable importance
        site_id <- c(site_id, test_df$ID)
        cv_fold <- c(cv_fold, test_df$FOLD)
        cv_obs <- c(cv_obs, test_df$LFMC_AVG)
        cv_pred <- c(cv_pred, predict(t_mdl, test_df[, mdl_vars], num.threads = ncores)$predictions)
        var_importance <- c(var_importance, list(importance(t_mdl)))
    }
    
    #* Save results ----
    cv_output <- list(
        "validation" = data.frame("ID" = site_id,
                                  "fold" = cv_fold,
                                  "obs" = cv_obs,
                                  "RFcv" = cv_pred),
        "variable_importance" = var_importance
    )
    file_name <- paste("./data/3fcv",
                       vopt + length(vars_df$name),
                       vars_df$name[vopt],
                       "filt.rds",
                       sep = "_")
    saveRDS(cv_output, file_name)
}
Sys.time() - t0

#file.remove(file.path(wdir, paste0("cve_", 1:25, ".rds")))

# END ---
