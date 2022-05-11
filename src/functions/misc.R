# *************************************************************************
# LFMC-RF - A remote sensing approach to predict LFMC at large scale using
# Random Forests
# Miscelaneaous functions
# Author: Angel Cunill Camprubi <acunill86@gmail.com>
# Last update: 03 Mar 2022
# *************************************************************************

# Standard error of the mean
se <- function(x) {
    sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
}

# Mean bias error
MBE <- function(obs, pred) {
    mean(pred - obs, na.rm = TRUE)
}

# Unbiased mean bias error
ubRMSE <- function(obs, pred) {
    m_obs <- mean(obs)
    m_pred <- mean(pred)
    sqrt(mean((obs - pred - (m_obs - m_pred)) ^ 2))
}

# Summary of accuracy metrics for the analysis of results
summarize <- function(df, obs = "obs", pred = "pred", interval = NULL) {
    df1 <- df[, c(obs, pred)]
    if(!is.null(interval)) {
        df1 <- df1[df1$obs > interval[1] & df1$obs <= interval[2], ]
    }
    names(df1) <- c("obs", "pred")
    list(
        with(df1, data.frame(
            "MBE" = MBE(obs = obs, pred = pred),
            "MAE" = MAE(y_true = obs, y_pred = pred),
            "RMSE" = RMSE(y_true = obs, y_pred = pred),
            "ubRMSE" = ubRMSE(obs = obs, pred = pred),
            "CCC" = DescTools::CCC(x = obs, y = pred)$rho.c[[1]],
            "Cb" = DescTools::CCC(x = obs, y = pred)$C.b[[1]],
            "VEcv" = R2_Score(y_true = obs, y_pred = pred),
            "R" = sqrt(R2_Score(y_true = obs, y_pred = pred)),
            "r" = cor(x = obs, y = pred),
            "nObs" = length(obs),
            "meanObs" = mean(obs),
            "medObs" = mean(obs),
            "sdObs" = sd(obs),
            "minObs" = min(obs),
            "maxObs" = max(obs),
            "rangeObs" = diff(range(obs))
        )),
        summary(lm(obs ~ pred, data = df1))
    )
}
