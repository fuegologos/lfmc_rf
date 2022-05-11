# *************************************************************************
# LFMC-RF - A remote sensing approach to predict LFMC at large scale using
# Random Forests
# Forward Feature Selection (FFS) function
# Adapted from Meyer et al. (2018)
# Author: Angel Cunill Camprubi <acunill86@gmail.com>
# Last update: 30 Oct 2021
# *************************************************************************

ffs_iter <- function(predictors,
                     response,
                     df,
                     spacevar = "SITE",
                     wdir,
                     iseed = 42,
                     ncores = 4,
                     niter = 25,
                     cvfolds = 5,
                     serule = FALSE) {
    # Packages
    require(ranger)
    require(CAST)
    require(doParallel)
    require(foreach)
    
    # Standard error function
    se <- function(x) {
        sd(x, na.rm = TRUE) / sqrt(length(na.exclude(x)))
    }
    
    # Initial settings
    names(df)[names(df) == response] <- "RESP"
    
    # Forward selection
    for (run in seq_len(length(predictors) - 1)) {
        
        cat(sprintf("\nRun %i\n", run))
        
        # Predictors combination
        if (run == 1) {
            pred_comb <- data.frame(t(combn(predictors, 2)))
        } else {
            pred_comb <- data.frame(
                "X1" = predictors[!predictors %in% sel_predictors]
            )
        }
        
        # Current run settings
        set.seed(iseed)
        seeds <- sample(1:1000, niter)
        n <- nrow(pred_comb)
        ncores <- ifelse(n < ncores, n, ncores)
        
        # Cross-validation iteration
        cv_error <- list()
        for (i in 1:niter) {
            cat(sprintf("\nIteration %i of %i: %s\n", i, niter, as.character(Sys.time())))
            
            indices <- CreateSpacetimeFolds(df,
                                            spacevar = spacevar,
                                            k = cvfolds,
                                            seed = seeds[i])
            
            cl <- makeCluster(ncores)
            registerDoParallel(cl)
            tmp <- foreach(comb_i = seq_len(n),
                           .packages = c("ranger")) %dopar% {
                
                if (run == 1) {
                    mdl_vars <- c("RESP", pred_comb$X1[comb_i], pred_comb$X2[comb_i])
                } else {
                    mdl_vars <- c("RESP", sel_predictors, pred_comb$X1[comb_i])
                }
                
                obs <- c()
                pred <- c()
                
                for (fold in seq_along(indices$index)) {
                    train_df <- df[indices$index[[fold]], mdl_vars]
                    test_df <- df[indices$indexOut[[fold]], mdl_vars]
                    mdl <- ranger(
                        formula = RESP ~ .,
                        data = train_df,
                        num.trees = 250,
                        mtry = 2,
                        oob.error = FALSE,
                        splitrule = "variance",
                        verbose = FALSE,
                        seed = 42,
                        num.threads = 1
                    )
                    obs <- c(obs, test_df$RESP)
                    pred <- c(pred, predict(mdl, test_df, num.threads = 1)$predictions)
                }
                sqrt(mean((obs - pred) ^ 2, na.rm = TRUE))
            }
            stopCluster(cl)
            cv_error <- c(cv_error, list(do.call(rbind, tmp)))
        }
        
        # Best model of the current run
        res <- do.call(cbind, cv_error)
        perf <- data.frame(
            "id" = seq_len(nrow(pred_comb)),
            "rmse" = apply(res, 1, mean, na.rm = TRUE),
            "sd" = apply(res, 1, sd, na.rm = TRUE),
            "se" = apply(res, 1, se)
        )
        best <- perf[perf$rmse == min(perf$rmse), ]
        best <- best[which.min(best$sd), ]
        
        # Accuracy threshold
        acc_thr <- ifelse(serule, best$rmse + best$se, best$rmse)
        
        # Update current best model
        if (run == 1) {
            sel_predictors <- c(pred_comb$X1[best$id], pred_comb$X2[best$id])
        } else {
            # Continue?
            if (acc_best <= acc_thr) {
                cat(sprintf("\nFinal model with %i predictors\n",
                            length(sel_predictors)))
                break
            }
            sel_predictors <- c(sel_predictors, pred_comb$X1[best$id])
        }
        acc_best <- best$rmse
        cat(sprintf("\nBest model: %s (RMSE %.2f)\n",
                    paste(sel_predictors, collapse = " "),
                    acc_best))
        
        # Save partial results
        saveRDS(list("predictors" = sel_predictors,
                     "performance" = best),
                file.path(wdir, paste0("run", run, ".rds")))
        
    }
    
    # Joint partial results
    partfiles <- list.files(wdir, pattern = "^run", full.names = TRUE)
    partres <- lapply(partfiles, readRDS)
    final_predictors <- lapply(partres, function(x) x$predictors)
    final_acc <- lapply(partres, function(x) x$performance)
    
    # Output
    list(
        "predictors" = final_predictors,
        "perfomance" = do.call(rbind, final_acc)
    )
}

# END ---
