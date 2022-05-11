# *************************************************************************
# LFMC-RF - A remote sensing approach to predict LFMC at large scale using
# Random Forests
# Model performance assessment independently from a specific hyperparameters
# NDVIcv filter applied to the whole dataset
# Author: Angel Cunill Camprubi <acunill86@gmail.com>
# Last update: 06 Mar 2022
# *************************************************************************

# SETTINGS ----------------------------------------------------------------

# Packages
library(ranger)
library(Rcpp)
library(DescTools)
library(MLmetrics)
library(CAST)

# Functions
source("src/functions/misc.R")

# Cross-validation folds
source("src/functions/create_llocv_folds.R")

# Others
options(stringsAsFactors = FALSE)
RNGkind("L'Ecuyer-CMRG")
ncores <- parallel::detectCores() - 1

# Cross-validation settings
niter <- 100               #number of iterations
k <- 5                     #cross-validation folds
pred_opt <- "selp"         #(allp, selp)
fname <- paste0(pred_opt, "")

# Create temporary directory
tmpdir <- paste0("./data/mp/", fname)
dir.create(tmpdir, recursive = TRUE)

# MODEL INPUTS ------------------------------------------------------------

# Predictors
si <- c('NDVI', 'EVI', 'SAVI', 'VARI', 'VIgreen', 'Gratio', 'NDII6', 'NDII7',
        'NDWI', 'GVMI', 'MSI', 'NDTI', 'STI')
bands <- c('NR1', 'NR2', 'NR3', 'NR4', 'NR5', 'NR6', 'NR7')
aux <- c("LST8day", "DOY_COS", "DOY_SIN")
predictors <- list(
    "allp" = c(si, bands, aux),
    "selp" = lapply(1, function(x) {
        p <- readRDS("data/vs/varsel_final.rds")$predictors
        p[[length(p)]]
    })[[1]]
)

# Filter
sp_filt <- c(1000, seq(0.2, 0.6, 0.05))
#sp_filt <- c(1000, seq(0.2, 0.45, 0.05)) #selp
#sp_filt <- c(1000, seq(0.25, 0.35, 0.05)) #allp

# Hyperparameters
num_trees <- c(500) #c(250, 500, 800, 1000)
min_node_size <- c(5:30)
sample_fraction <- seq(0.2, 0.95, 0.05)

# MODEL ASSESSMENT --------------------------------------------------------
# 35 hours (selp)

# Subset of variables
vars <- predictors[[pred_opt]]
mdl_vars <- c("LFMC", vars)

# Load database
db0 <- readRDS("./data/db03_lfmc.rds")
db0 <- db0[complete.cases(db0[, mdl_vars]), ]

# Hyperparameter combinations
max_mtry <- ceiling(length(vars) * 0.4)
max_mtry <- ifelse(max_mtry < 4, 4, max_mtry)
hyp0 <- expand.grid(
    "num_trees" = num_trees,
    "mtry" = seq(2, max_mtry, 1),
    "min_node_size" = min_node_size,
    "sample_fraction" = sample_fraction
)

# Seeds for outer cross-validation folds
set.seed(42)
seeds <- sample(1:9999, niter)

# Evaluate models across filters
for (spf_i in seq_along(sp_filt)) {
    
    cat(sprintf("\nFilter %i of %i: %s\n",
                spf_i,
                length(sp_filt),
                as.character(Sys.time())))
    
    # Apply NDVIcv filter
    if (sp_filt[spf_i] != 1000) {
        db <- db0[complete.cases(db0$NDVIcv), ]
        db <- db[db$NDVIcv < sp_filt[spf_i], ]
    } else {
        db <- db0
    }
    
    # Iterative cross-validation process
    for (iter_i in seq_len(niter)) {
        
        cat(sprintf("\nIteration %i of %i: %s\n",
                    iter_i,
                    niter,
                    as.character(Sys.time())))
        
        # Setting outputs
        cv_obs <- c()
        cv_pred <- c()
        out_hyp <- list()
        
        # Define outer cross-validation folds
        out_idx <- CreateSpacetimeFolds(db,
                                        spacevar = "SITE",
                                        k = k,
                                        seed = seeds[iter_i])
        
        # Set seeds for hyperparameter selection
        set.seed(seeds[iter_i])
        seeds2 <- sample(1:1000, k)
        
        # Nested cross-validation
        for (outfold in seq_along(out_idx$index)) {
            
            # Split data into train/test
            train_df <- db[out_idx$index[[outfold]], ]
            test_df <- db[out_idx$indexOut[[outfold]], ]
            
            # Random subset of hyperparameters
            set.seed(seeds2[outfold])
            hyp <- hyp0[sample(nrow(hyp0), 50), ]
            hyp_order <- order(hyp$mtry, hyp$num_trees, hyp$min_node_size)
            hyp <- hyp[hyp_order, ]
            
            # Inner calibration
            hyp_error <- c()
            for (h in seq_len(nrow(hyp))) {
                
                # Calibration settings
                hyp_i <- hyp[h, ]
                obs <- c()
                pred <- c()
                
                # Create indices for nested-cv folds
                inn_idx <- CreateSpacetimeFolds(train_df,
                                                spacevar = "SITE",
                                                k = k,
                                                seed = 42)
                
                # Train model
                for (innfold in seq_along(inn_idx$index)) {
                    
                    # Split training data into tune/assess
                    tune_df <- train_df[inn_idx$index[[innfold]], mdl_vars]
                    asses_df <- train_df[inn_idx$indexOut[[innfold]], mdl_vars]
                    
                    # Model tune
                    mdl <- ranger(
                        formula = LFMC ~ .,
                        data = tune_df,
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
                    
                    # Inner predictions
                    obs <- c(obs, asses_df$LFMC)
                    pred <- c(pred, predict(mdl,
                                            asses_df,
                                            num.threads = ncores)$predictions)
                }
                
                # Calculate calibration error
                hyp_error <- c(hyp_error,
                               RMSE(y_true = obs, y_pred = pred))
            }
            
            # Calibrate best model with the training data
            best_hyp <- hyp[which.min(hyp_error), ]
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
            
            # Outer predictions
            cv_obs <- c(cv_obs, test_df$LFMC)
            cv_pred <- c(cv_pred, predict(best_mdl,
                                          test_df[, mdl_vars],
                                          num.threads = ncores)$predictions)
            
            # Output hyperparameters
            out_hyp <- c(out_hyp, list(unlist(best_hyp)))
        }
        
        # Calculate model performance metrics
        mdl_perf <- data.frame(
            "MBE" = MBE(obs = cv_obs, pred = cv_pred),
            "MAE" = MAE(y_true = cv_obs, y_pred = cv_pred),
            "RMSE" = RMSE(y_true = cv_obs, y_pred = cv_pred),
            "ubRMSE" = ubRMSE(obs = cv_obs, pred = cv_pred),
            "CCC" = CCC(x = cv_pred, y = cv_obs)$rho.c[[1]],
            "Cb" = CCC(x = cv_pred, y = cv_obs)$C.b[[1]],
            "VEcv" = R2_Score(y_true = cv_obs, y_pred = cv_pred),
            "n" = length(cv_obs),
            "sdObs" = sd(cv_obs),
            "minObs" = min(cv_obs),
            "maxObs" = max(cv_obs)
        )
        
        # Save partial results from current filter application
        filt_dir <- file.path(tmpdir, paste0("filt", spf_i))
        dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)
        iter_filename <- paste0("mp_iter", iter_i, ".rds")
        iter_output <- list(
            "options" = c(pred_opt, ifelse(sp_filt[spf_i] == 1000, "nofilt", "filt"), niter),
            "performance" = mdl_perf,
            "hyperparams" = do.call(rbind, out_hyp),
            "spfilt" = sp_filt[spf_i]
        )
        saveRDS(iter_output, file.path(filt_dir, iter_filename))
    }
    
    # Load final outputs from current filter application
    part_files <- paste0("mp_iter", 1:niter, ".rds")
    part_res <- lapply(part_files, function(x) readRDS(file.path(filt_dir, x)))
    
    # Model performance from current filter application
    res_mp <- do.call(rbind, lapply(part_res, function(x) x$performance))
    saveRDS(res_mp, file.path(tmpdir, paste0("mp_", pred_opt, spf_i, ".rds")))
    
    # Selected hyperparameters from current filter application
    res_hyp <- do.call(rbind, lapply(part_res, function(x) x$hyperparams))
    res_hyp <- cbind("spfilt" = sp_filt[spf_i], res_hyp)
    saveRDS(res_hyp, file.path(tmpdir, paste0("hyp_", pred_opt, spf_i, ".rds")))
}

# RESULTS -----------------------------------------------------------------

# Calculate summary metrics from NDVIcv filter results
out <- lapply(list(mean, se, min, max), function(x) {
    tmp <- lapply(seq_along(sp_filt), function(spf_i) {
        mpfile <- file.path(tmpdir, paste0("mp_", pred_opt, spf_i, ".rds"))
        mp <- readRDS(mpfile)
        apply(mp, 2, x)
    }) 
    mp <- do.call(rbind, tmp)
    as.data.frame(cbind("spfilt" = sp_filt, mp))
})
names(out) <- c("mean", "se", "min", "max")
saveRDS(out, paste0("data/mp/mp_", fname, ".rds"))

# Hyperparameters of the best result
best <- 2
hyp <- readRDS(file.path(tmpdir, paste0("hyp_", pred_opt, "_hom", best, ".rds")))
summary(hyp)
{
    par(mfrow = c(2, 2))
    for (i in 2:5) {
        step <- ifelse(i %in% c(1, 5), 0.025,
                       ifelse(i %in% c(3, 4), 0.5, 0))
        if (step == 0) {
            brk <- unique(hyp[, i])
        } else {
            brk <- c(min(hyp[, i]) - step, sort(unique(hyp[, i])) + step)
        }
        hist(
            hyp[, i],
            breaks = brk,
            main = colnames(hyp)[i]
        )
    }
    par(op)
}

# END ---
