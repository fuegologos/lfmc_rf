# *************************************************************************
# LFMC-RF - A remote sensing approach to predict LFMC at large scale using
# Random Forests
# Analysis of the results
# Author: Angel Cunill Camprubi <acunill86@gmail.com>
# Last update: 13 Mar 2022
# *************************************************************************

# SETTINGS ----------------------------------------------------------------

# Packages
library(stringr)
library(lubridate)
library(grid)
library(ggplot2)
library(ggpointdensity) 
library(ggridges)
library(patchwork)
library(ggh4x)                 # Conditional area time-series
library(cowplot)               # Extract ggplot legend/axis_canvas
library(RColorBrewer)
library(MASS)                  # True histogram function
library(PerformanceAnalytics)  # Correlogram and moments
library(corrplot)              # Correlogram
library(DescTools)
library(MLmetrics)
library(sf)
library(raster)

# Functions
source("src/functions/misc.R")

# Others
options(stringsAsFactors = FALSE)

# GRAPHICAL PARAMETERS ----------------------------------------------------

op <- par(no.readonly = TRUE)
convertUnit(unit(.28, "points"), "mm", valueOnly = TRUE)
convertUnit(unit(60, "mm"), "inches", valueOnly = TRUE)
convertUnit(unit(190, "mm"), "inches", valueOnly = TRUE)

# Line width (mm_real * 1/0.75 = mm_R)
.cf <- 1/.75
lw_s <- .18 * .cf
lw_m <- .25 * .cf
lw_l <- .35 * .cf
lw_xl <- .6 * .cf

# Colour palette
YlGnBu <- colorRampPalette(brewer.pal(3, "YlGnBu"))
GnBu <- colorRampPalette(brewer.pal(3, "GnBu"))
red_line <- "#de2d26"

# Custom theme (lfmc)
custom_theme <-
    theme_bw(
        base_size = 8,
        base_line_size = lw_m,
        base_rect_size = 0
    ) +
    theme(
        axis.ticks.length = unit(1, "mm"),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(margin = margin(.6, 0, 0, 0, "mm")),
        axis.text.y = element_text(margin = margin(0, .6, 0, 0, "mm")),
        axis.title.x.top = element_text(margin = margin(0, 0, 1, 0, "mm")),
        axis.title.x.bottom = element_text(margin = margin(1, 0, 0, 0, "mm")),
        axis.title.y.left = element_text(margin = margin(0, .6, 0, 0, "mm")),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, size = lw_m * 1.5),
        panel.ontop = FALSE,
        plot.background = element_blank()
    )

# LOAD DATABASE -----------------------------------------------------------

# Load database
db <- readRDS("./data/db03_lfmc.rds")

# Features collection
si <- c('NDVI', 'EVI', 'SAVI', 'VARI', 'VIgreen', 'Gratio', 'NDII6', 'NDII7',
        'NDWI', 'GVMI', 'MSI', 'NDTI', 'STI')
bands <- c('NR1', 'NR2', 'NR3', 'NR4', 'NR5', 'NR6', 'NR7')
aux <- c("LST8day", "DOY_COS", "DOY_SIN")
mdl_vars <- c("LFMC", si, bands, aux)

# Missing data due to covariates
n0 <- nrow(db)
n1 <- sum(complete.cases(db[, c(si, bands)]))
1 - (n1 / n0)
n2 <- sum(complete.cases(db[, c(si, bands, aux)]))
1 - (n2 / n1)
n3 <- sum(complete.cases(db[, c(si, bands, "DOY_COS", "DOY_SIN", "LST")]))
1 - (n3 / n1)
length(unique(db[complete.cases(db[, c(si, bands)]), "SITE"]))

# Remove missing data from covariates
db <- db[complete.cases(db[, mdl_vars]), ]
nrow(db)
length(unique(db$SITE))

# LFMC DATA DESCRIPTION ---------------------------------------------------

# Summaries of LFMC
summary(db)
quantile(db$LFMC, probs = c(.05, .95))
sd(db$LFMC)
skewness(db$LFMC) # positive skew
kurtosis(db$LFMC) # leptokurtic (3 = mesokurtic)

# Plot settings
descr_theme <- theme(
    axis.title.x.bottom = element_text(margin = margin(-2, 0, 0, 0, "mm")),
    panel.grid = element_blank(),
    plot.margin = margin(0, 1, 0, 1, "mm")
)

# LFMC univariate distribution
lfmc_hist <- ggplot(data = db, aes(x = LFMC)) +
    geom_histogram(aes(y = ..density..),
                   binwidth = 5,
                   fill = GnBu(4)[3],
                   color = gray(.95),
                   size = .04 * .cf) +
    scale_x_continuous(name = "LFMC [%]",
                       expand = c(.02, 0)) +
    scale_y_continuous(name = "Density",
                       expand = c(.02, 0)) +
    annotate("text",
             x = Inf,
             y = Inf,
             label = "(a)",
             size = 10 / .pt, #mm
             vjust = 1.5,
             hjust = 1.3) +
    custom_theme +
    descr_theme

# Number of records by country
ncountry <- as.data.frame(sort(table("country" = db$COUNTRY), decreasing = TRUE))
ncountry$perc <- ncountry$Freq / sum(ncountry$Freq)
ncountry$cumsum <- cumsum(ncountry$perc)
ncountry$count_ord <- seq_len(nrow(ncountry))

# LFMC by country
db$count_ord <- 
    as.factor(ncountry$count_ord[match(db$COUNTRY, levels(ncountry$country))])
xlabs <- paste0(ncountry$country, "\n(", ncountry$Freq, ")")
lfmc_violin <- ggplot(data = db,
                      aes(x = count_ord, y = LFMC, fill = count_ord),
                      color = "grey20") +
    geom_violin(draw_quantiles = .5,
                size = lw_s) +
    labs(x = NULL, y = "LFMC [%]") +
    annotate("text",
             x = Inf,
             y = Inf,
             label = "(b)",
             size = 10 / .pt,
             vjust = 1.5,
             hjust = 1.3) +
    scale_y_continuous(expand = c(.04, .04), limits = c(20, 250)) +
    scale_x_discrete(labels = xlabs) +
    scale_fill_brewer(type = "seq", palette = "GnBu") +
    custom_theme +
    descr_theme +
    theme(
        legend.position = "none",
    )

# Output fig
fig <- lfmc_hist | lfmc_violin
ggsave(plot = fig,
       filename = "figS1.pdf",
       device = "pdf",
       path = "./output",
       width = 146,
       height = 68,
       units = "mm")

# CORRELATION PLOTS -------------------------------------------------------

# Correlogram with scatterplots
# png("./output/figs/00-correlogram.png", width = 28, height = 20, units = "cm", res = 150)l
# chart.Correlation(db[, mdl_vars], histogram = FALSE)
# dev.off()

# Correlation matrix plot (supplementary)
vars <- c("LFMC", si, bands, c("LST", "LST8day", "DOY_COS", "DOY_SIN"))
lfmc_corr <- round(cor(db[complete.cases(db[, vars]), vars]), 2)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
pdf("./output/figS3.pdf", width = 6, height = 6)
corrplot(
    corr = lfmc_corr,
    method = "circle",
    type = "upper",
    col = col(200),
    diag = FALSE,
    addgrid.col = "grey30",
    tl.cex = .7,
    tl.col = "black",
    cl.cex = .7,
    cl.align.text = "l"
    # addCoef.col = "gray20",
    # number.cex = .4
)
dev.off()

# VARIABLE SELECTION ------------------------------------------------------

# Comparative of daily and composite LST products
# vs <- readRDS("data/vs/varsel1.rds")
# vs <- readRDS("data/vs/varsel2.rds")
# vs <- readRDS("data/vs/varsel3.rds")
# vs <- readRDS("data/vs/varsel4.rds")

# Selected variables
vs <- readRDS("data/vs/varsel_final.rds")
#vs <- readRDS("data/vs/varsel_final-alldata.rds")

# Stepwise results
vs_perf <- vs$perfomance
vs_perf$dif <- 0
for (i in 2:nrow(vs_perf)) {
    vs_perf$dif[i] <- vs_perf$rmse[i-1] - vs_perf$rmse[i]
}
vs_perf$cum <- cumsum(vs_perf$dif)

# Plot of selected variables
pred <- vs$predictors[[length(vs$predictors)]]
pred[pred == "LST8day"] <- "LST"
xlabs <- c(str_c(pred[1:2], collapse = "\n+"), paste0("+", pred[3:length(pred)]))

vs <- vs$perfomance
vs$id <- seq_len(nrow(vs))
upper <- vs$rmse + vs$se
lower <- vs$rmse - vs$se

fig <- ggplot(data = vs) +
    geom_line(aes(x = id, y = rmse),
              size = lw_m) +
    geom_errorbar(aes(x = id, ymin = lower, ymax = upper),
                  width = .15,
                  size = lw_s) +
    geom_point(aes(x = id, y = rmse),
               shape = 16,
               size = .8) +
    labs(x = NULL, y = "RMSE [%]") +
    scale_x_continuous(
        labels = xlabs,
        breaks = 1:length(xlabs),
        expand = c(.02, .02)
    ) +
    # scale_y_continuous(
    #     limits = c(19.5, 22.1),
    #     #limits = c(21.5, 25.1),
    #     breaks = seq(19, 25, 1),
    #     expand = c(0, 0)
    # ) +
    custom_theme +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1
        ),
        plot.margin = margin(.5, .5, .5, .5, "mm")
    )
fig

# Save plot
ggsave(plot = fig,
       filename = "fig3.pdf",
       device = "pdf",
       path = "./output",
       width = 90,
       height = 60,
       units = "mm")

# MODEL PERFORMANCE -------------------------------------------------------

# Load model performance results
df1 <- readRDS("data/mp/mp_allp.rds")
df2 <- readRDS("data/mp/mp_selp.rds")
df3 <- readRDS("data/mp/mp_selp-train.rds")

# Summary of accuracy metrics from all fitted models (MP iterations)
round(apply(rbind(df1$min, df2$min, df3$min), 2, min), 2)
round(apply(rbind(df1$max, df2$max, df3$max), 2, max), 2)

# Bias and large errors analysis
bias <- do.call(rbind, lapply(list(df1$mean, df2$mean), function(df) {
    d1 <- df$RMSE - df$ubRMSE
    d2 <- df$RMSE - df$MAE
    data.frame(
        "minB" = min(d1),
        "maxB" = max(d1),
        "minL" = min(d2),
        "maxL" = max(d2)
    )
}))
min(bias$minB); max(bias$maxB)
min(bias$minL); max(bias$maxL)

# NDVIcv FILTER ANALYSIS --------------------------------------------------

# Load model performance results
df1 <- readRDS("data/mp/mp_allp.rds")
df2 <- readRDS("data/mp/mp_selp.rds")

# Explore NDVIcv filter application
db2 <- db[complete.cases(db$NDVIcv), ]
(1 - nrow(db2) / nrow(db)) * 100
(1 - sum(db2$NDVIcv < 0.3) / nrow(db2)) * 100
(1 - sum(db2$NDVIcv < 0.35) / nrow(db2)) * 100
idsNDVIcv <- db2[db2$NDVIcv >= 0.3, "ID"]
db2 <- db2[db2$NDVIcv < 0.3, ]
length(unique(db2$SITE))
summary(db2$LFMC)
remove(db2)

# Differences between Selp-nofilter and Selp-filter
df <- df2$mean
(df[1, "sdObs"] - df[4, "sdObs"]) / df[4, "sdObs"]
(df[1, "RMSE"] - df[4, "RMSE"]) / df[4, "RMSE"]
df_ci <- qt(0.975, 99) * df2$se
df_cilow <- df - df_ci
df_ciupp <- df + df_ci
round(df[c(1, 4), ], 3)
round(df_cilow[c(1, 4), ], 3)
round(df_ciupp[c(1, 4), ], 3)
round(df_ci, 3)

# Plot results
plot_filt <- function(data1, data2, variable, ylabel, xlabel = NULL, plot_tag) {
    df1 <- lapply(data1, function(x) {
        s <- x[, c("spfilt", variable)]
        names(s) <- c("filt", "xvar")
        s[s$filt != 1000, ]
    })
    m1 <- df1$mean
    se1 <- df1$se
    upper1 <- m1$xvar + se1$xvar
    lower1 <- m1$xvar - se1$xvar
    df2 <- lapply(data2, function(x) {
        s <- x[, c("spfilt", variable)]
        names(s) <- c("filt", "xvar")
        s[s$filt != 1000, ]
    })
    m2 <- df2$mean
    se2 <- df2$se
    upper2 <- m2$xvar + se2$xvar
    lower2 <- m2$xvar - se2$xvar
    ggplot() +
        geom_line(
            data = m1,
            aes(x = filt, y = xvar),
            size = lw_m,
            colour = GnBu(4)[4]
        ) +
        geom_errorbar(
            data = m1,
            aes(x = filt, ymin = lower1, ymax = upper1),
            width = .01,
            size = lw_s,
            colour = GnBu(4)[4]
        ) +
        geom_point(
            data = m1,
            aes(x = filt, y = xvar),
            shape = 16,
            size = .8,
            colour = GnBu(4)[4]
        ) +
        geom_line(
            data = m2,
            aes(x = filt, y = xvar),
            size = lw_m,
            colour = red_line
        ) +
        geom_errorbar(
            data = m2,
            aes(x = filt, ymin = lower2, ymax = upper2),
            width = .01,
            size = lw_s,
            colour = red_line
        ) +
        geom_point(
            data = m2,
            aes(x = filt, y = xvar),
            shape = 16,
            size = .8,
            colour = red_line
        ) +
        labs(x = xlabel, y = ylabel) +
        scale_x_continuous(
            limits = c(0.18, 0.47),
            breaks = seq(0.1, 0.5, 0.05),
            expand = c(0, 0)
        ) +
        scale_y_continuous(expand = c(.01, .01)) +
        annotate(
            "text",
            x = Inf,
            y = Inf,
            label = plot_tag,
            size = 10 / .pt,
            vjust = 1.5,
            hjust = 1.3
        ) +
        custom_theme +
        theme(
            axis.text.x = element_text(
                angle = 90,
                vjust = 0.5,
                hjust = 1
            ),
            panel.grid = element_blank(),
            plot.margin = margin(0, 1, 1, 0, "mm")
        )
}
p1 <- plot_filt(df1, df2, "MBE", "MBE [%]", plot_tag = "(a)")
p2 <- plot_filt(df1, df2, "RMSE", "RMSE [%]", plot_tag = "(b)")
p3 <- plot_filt(df1, df2, "CCC", "CCC", "NDVIcv", plot_tag = "(c)")
p4 <- plot_filt(df1, df2, "VEcv", "VEcv", "NDVIcv", plot_tag = "(d)")

# Save plot
fig <- p1 + p2 + p3 + p4 + plot_layout(nrow = 2)
ggsave(plot = fig,
       filename = "figS4.pdf",
       device = "pdf",
       path = "./output",
       width = 150,
       height = 110,
       units = "mm")

# MODEL CALIBRATION AND EXTRAPOLATION -------------------------------------

# Load data
calext <- readRDS("data/calext/nofilt.rds")

# Accuracy metrics (calibration)
summarize(calext$dev)[[1]]

# Accuracy metrics (extrapolation)
summarize(calext$val)[[1]]

# Accuracy metrics (critical LFMC)
summarize(calext$dev, interval = c(30, 120))
summarize(calext$val, interval = c(30, 120))

# Hyperparameters
calext$hyp

# Data summaries
calext$sum
d1 <- unique(db[db$ID %in% calext$dev$ID, "SITE"])
length(d1)
d2 <- unique(db[db$ID %in% calext$val$ID, "SITE"])
length(d2)
sum(d2 %in% d1)

# PREDICTIONS OUTSIDE NDVIcv THRESHOLD ------------------------------------

# Predictions dataframe
calext <- readRDS("data/calext/nofilt.rds")
df <- calext$dev[, c("ID", "obs", "pred")]
s <- summarize(df)[[1]]

# Predictions rejected for the theoric NDVIcv filter
idsNDVIcv <- na.exclude(db[db$NDVIcv > 0.30, "ID"])
predNDVIcv <- df[df$ID %in% idsNDVIcv, ]
predNDVIcv$diff <- predNDVIcv$obs - predNDVIcv$pred
predNDVIcv$cols <- ifelse(abs(predNDVIcv$diff) > s$MAE, red_line, GnBu(4)[4])
{
    pdf("./output/figS5.pdf", width = 4, height = 4)
    d <- summarize(predNDVIcv)
    lbl <- sprintf("n: %i\nMBE: %.02f%%\nRMSE: %.02f%%\nMAE: %.02f%%\nCCC: %.02f\nCb: %.02f",
                   nrow(predNDVIcv), d[[1]]$MBE, d[[1]]$RMSE, d[[1]]$MAE, d[[1]]$CCC, d[[1]]$Cb)
    par(mar = c(4, 4, 1, 1), mgp = c(1.9, .8, 0))
    plot(predNDVIcv$pred,
         predNDVIcv$obs,
         col = predNDVIcv$cols,
         xlim = c(0, 200),
         ylim = c(0, 200),
         xaxs = "i", yaxs = "i",
         pch = 16,
         xlab = "Predicted LFMC [%]",
         ylab = "Measured LFMC [%]")
    text(x = 2, y = 198, labels = lbl, adj = c(0, 1), cex = 0.8)
    legend(
        x = "bottomright",
        legend = c("dif > MAE", "dif <= MAE"),
        pch = 16,
        col = c(red_line, GnBu(4)[4]),
        cex = 0.8
    )
    abline(0, 1)
    par(op)
    dev.off()
}

# OBSERVED VS PREDICTED ---------------------------------------------------

# Load data
calext <- readRDS("data/calext/nofilt.rds")
cal <- calext$dev
ext <- calext$val
rtm <- db[complete.cases(db$RTM), c("ID", "RTM", "SITE")]
rtm <- rtm[rtm$RTM > 0, ]
rtm_cal <- merge(cal, rtm, by = "ID")
rtm_ext <- merge(ext, rtm, by = "ID")

# Observed vs predicted plot function
plot_op <- function(df, lfmc_obs, lfmc_RF, plotid, xlabel = FALSE, ylabel = FALSE) {
    df <- df[, c(lfmc_obs, lfmc_RF)]
    form <- formula(paste(lfmc_obs, lfmc_RF, sep = " ~ "))
    cv_lm <- lm(form, data = df)
    cv_lm_pred <- data.frame(c(min(df[, lfmc_RF]), max(df[, lfmc_RF])))
    names(cv_lm_pred) <- lfmc_RF
    cv_lm_pred$y <- predict(cv_lm, cv_lm_pred)
    cv_fun <- sprintf("y=%.2f+%.2fx", cv_lm$coefficients[1], cv_lm$coefficients[2])
    names(df)[names(df) %in% c(lfmc_obs, lfmc_RF)] <- c("obs", "pred")
    p <- ggplot(data = df, aes(x = pred, y = obs)) +
        geom_pointdensity(size = .3, adjust = 1) +
        #geom_hex(bins = 50) +
        geom_abline(slope = 1,
                    intercept = 0,
                    colour = "grey40",
                    linetype = "dashed",
                    size = lw_m) +
        geom_segment(x = cv_lm_pred[1, lfmc_RF],
                     y = cv_lm_pred[1, "y"],
                     xend = cv_lm_pred[2, lfmc_RF],
                     yend = cv_lm_pred[2, "y"],
                     colour = red_line,
                     linetype = "longdash",
                     size = lw_m) +
        annotate("text",
                 x = 262,
                 y = 10,
                 label = cv_fun,
                 size = 8 / .pt,
                 vjust = 0,
                 hjust = 1) +
        annotate("text",
                 x = 8,
                 y = 262, #292
                 label = plotid,
                 size = 10 / .pt,
                 vjust = 1,
                 hjust = 0) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 270)) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 270)) +
        scale_colour_viridis_c() +
        guides(colour = guide_colourbar(
            title = NULL,
            barwidth = unit(2, "mm"),
            barheight = unit(16, "mm"),
            frame.linewidth = lw_m,
            ticks.linewidth = lw_m,
            frame.colour = "gray30",
            ticks.colour = "gray90",
            draw.ulim = FALSE,
            draw.llim = FALSE,
            default.unit = "mm"
        )) +
        custom_theme +
        theme(
            panel.grid = element_blank(),
            legend.position = c(1, 0),
            legend.justification = c(-.1, 0),
            legend.margin = margin(0, 0, 0, 0, "mm"),
            legend.text = element_text(margin = margin(0, 0, 0, -.5, "mm")),
            legend.background = element_blank(),
            plot.margin = margin(0, 6.5, 2, 0, "mm")
        )
    if(xlabel) {
        p <- p + xlab("Predicted LFMC [%]")
    } else {
        p <- p + xlab(NULL)
    }
    if(ylabel) {
        p <- p + ylab("Measured LFMC [%]")
    } else {
        p <- p + ylab(NULL)
    }
    p
}

# OP output plots
p1 <- plot_op(cal, "obs", "pred", "(a) CAL", ylabel = TRUE)
p2 <- plot_op(rtm_cal, "obs", "pred", "(b) LFMCRF")
p3 <- plot_op(rtm_cal, "obs", "RTM", "(c) RTM")
p4 <- plot_op(ext, "obs", "pred", "(d) EXT", xlabel = TRUE, ylabel = TRUE)
p5 <- plot_op(rtm_ext, "obs", "pred", "(e) LFMCRF", xlabel = TRUE)
p6 <- plot_op(rtm_ext, "obs", "RTM", "(f) RTM", xlabel = TRUE)
fig <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(widths = unit(c(48, 48, 48), "mm"),
                                                 heights = unit(c(48, 48), "mm"),
                                                 ncol = 3)
ggsave(plot = fig,
       filename = "fig4.pdf",
       device = "pdf",
       path = "./output",
       width = 190, #160
       height = 150, #160
       units = "mm")

# RF vs RTM summaries
lapply(list(rtm_cal, rtm_ext), function(dat) {
    do.call(rbind, lapply(c("pred", "RTM"), function(x) {
        cbind(
            "method" = ifelse(x == "pred", "RF", x),
            summarize(dat, obs = "obs", pred = x)[[1]]
        )
    }))    
})
length(unique(rtm_cal$SITE))
length(unique(rtm_ext$SITE))
summarize(rtm, obs = "obs", pred = "pred")[[2]]

# Inflexion point under/over-prediction of RTM
lapply(list(rtm_cal, rtm_ext), function(dat) {
    rtm_lm <- lm(obs ~ RTM, data = dat)
    coef(rtm_lm)[[1]] / (1 - coef(rtm_lm)[[2]])    
})

# ERROR BY INTERVALS OF LFMC ----------------------------------------------

# Load data
calext <- readRDS("data/calext/nofilt.rds")
cal <- calext$dev
ext <- calext$val

# Calculate interval errors for CAL
cal$err <- cal$pred - cal$obs
cal$int <- cut(
    cal$obs,
    breaks = c(-Inf, 30, 120, Inf),
    labels = c("<30", "30-120", ">120"),
    include.lowest = FALSE,
    right = TRUE,
    ordered_result = TRUE
)
cal_count <- data.frame(table(cal$int))
cal_count$perc <- (cal_count$Freq / sum(cal_count$Freq)) * 100
cal_count

# Calculate interval errors for EXT
ext$err <- ext$pred - ext$obs
ext$int <- cut(
    ext$obs,
    breaks = c(-Inf, 30, 120, Inf),
    labels = c("<30", "30-120", ">120"),
    include.lowest = FALSE,
    right = TRUE,
    ordered_result = TRUE
)
ext_count <- data.frame(table(ext$int))
ext_count$perc <- (ext_count$Freq / sum(ext_count$Freq)) * 100
ext_count

# Counts
counts <- paste(cal_count$Freq, ext_count$Freq, sep = " / ")

# Scatterplot function with marginal densities
plot_scatden <- function(df0, plot_tag, plot_legend = TRUE, yaxis = TRUE) {
    scatter_p <- ggplot(df0,
                        aes(x = obs, y = err, color = int, fill = int)) +     
        geom_hline(yintercept = 0, size = lw_m, linetype = "dashed") +
        geom_point(
            size = 1, shape = 21, stroke = .2 
        ) +
        annotate(
            "text",
            x = Inf,
            y = Inf,
            label = plot_tag,
            size = 10 / .pt,
            vjust = 1.4,
            hjust = 1.3
        ) +
        scale_color_manual(
            values = c("#56B4E9", "#E69F00", "#009E73"),
            name = NULL,
            labels = paste0(levels(df0$int), "%"),
            drop = FALSE
        ) +
        scale_fill_manual(
            values = c("#56B4E980", "#E69F0080", "#009E7380"),
            name = NULL,
            labels = paste0(levels(df0$int), "%"),
            drop = FALSE
        ) +
        scale_x_continuous(
            limits = c(15, 250), expand = c(0, 0),
            name = "LFMC [%]"
        ) +
        scale_y_continuous(
            limits = c(-150, 150), expand = c(0, 0)
        ) +
        custom_theme +
        theme(
            panel.grid = element_blank(),
            plot.margin = margin(0, 0, 0, 0, "mm")
        )
    if (plot_legend) {
        scatter_p <- scatter_p + theme(legend.position = c(.01, .01),
                                       legend.justification = c(0, 0),
                                       legend.text = element_text(size = 8),
                                       legend.box.spacing = unit(0, "mm"),
                                       legend.key.size = unit(3.5, "mm"))
    } else {
        scatter_p <- scatter_p + theme(legend.position = "none")
    }
    if (yaxis) {
        scatter_p <- scatter_p + ylab("Predicted - Observed [%]")
    } else {
        scatter_p <- scatter_p + ylab(NULL)
    }
    # Compute densities
    densX <- ggplot2:::compute_density(df0$obs, NULL)
    densY <- ggplot2:::compute_density(df0$err, NULL)
    # Plots
    xdens <- axis_canvas(scatter_p, axis = "x") +
        geom_density_line(
            data = densX,
            aes(x = x, y = density),
            colour = "gray40",
            fill = "#d6d6d680",
            stat = "identity",
            size = .2
        ) +
        scale_x_continuous(limits = c(15, 250), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, max(densX$density) * 1.01), expand = c(0, 0))
    ydens <- axis_canvas(scatter_p, axis = "y") +
        geom_density_line(
            data = densY,
            aes(x = x, y = density),
            colour = "gray40",
            fill = "#d6d6d680",
            stat = "identity",
            size = .2
        ) +
        scale_x_continuous(limits = c(-150, 150), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, max(densY$density) * 1.01), expand = c(0, 0)) +
        coord_flip()
    p1 <- insert_xaxis_grob(
        plot = scatter_p,
        grob = xdens,
        height = grid::unit(10, "mm"),
        position = "top"
    )
    p2 <- insert_yaxis_grob(
        plot = p1,
        grob = ydens,
        width = grid::unit(10, "mm"),
        position = "right"
    )
    ggdraw(p2)
}

# Plots of residuals
p1 <- plot_scatden(cal, "(a)")
p2 <- plot_scatden(ext, "(b)", plot_legend = FALSE, yaxis = FALSE)
fig <- p1 | p2
ggsave(plot = fig,
       filename = "figS6.pdf",
       device = "pdf",
       path = "./output",
       width = 180,
       height = 80,
       units = "mm")

# Error distribution at the critical LFMC
quantile(cal[, "err"], probs = c(.05, .5, .95))
quantile(ext[, "err"], probs = c(.05, .5, .95))
quantile(cal[cal$int == "30-120", "err"], probs = c(.05, .5, .95))
quantile(ext[ext$int == "30-120", "err"], probs = c(.05, .5, .95))
with(cal[cal$obs < 30, ], MBE(obs = obs, pred = pred))
with(cal[cal$obs > 120, ], MBE(obs = obs, pred = pred))
with(ext[ext$obs > 120, ], MBE(obs = obs, pred = pred))

# ERRORS BY FUEL TYPE AND MONTH -------------------------------------------

# Load data
calext <- readRDS("data/calext/nofilt.rds")

# Assign LC to predictions data
cal <- merge(calext$dev,
             db[, c("ID", "MONTH", "SITE", "IGBP_NAME", "FUEL")],
             by = "ID")
ext <- merge(calext$val,
             db[, c("ID", "MONTH", "SITE", "IGBP_NAME", "FUEL", "LFMC")],
             by = "ID")

# Overall performance by LC
lcp_fun <- function(df) {
    tmp <- lapply(unique(df$FUEL), function(lc) {
            df_i <- df[df$FUEL == lc, ]
            data.frame(
                "LC" = lc,
                "MBE" = MBE(obs = df_i$obs, pred = df_i$pred),
                "MAE" = MAE(y_true = df_i$obs, y_pred = df_i$pred),
                "RMSE" = RMSE(y_true = df_i$obs, y_pred = df_i$pred),
                "ubRMSE" = ubRMSE(obs = df_i$obs, pred = df_i$pred),
                "CCC" = DescTools::CCC(x = df_i$obs, y = df_i$pred)$rho.c[[1]],
                "VEcv" = R2_Score(y_true = df_i$obs, y_pred = df_i$pred),
                "Cb" = DescTools::CCC(x = df_i$obs, y = df_i$pred)$C.b[[1]],
                "r" = cor(df_i$obs, df_i$pred),
                "N" = nrow(df_i),
                "Sites" = length(unique(df_i$SITE)),
                "SDobs" = sd(df_i$obs),
                "minobs" = min(df_i$obs),
                "maxobs" = max(df_i$obs)
            )
    })
    out <- do.call(rbind, tmp)
    out[order(out$LC), ]
}

lcp_cal <- lcp_fun(cal)
lcp_cal$Nperc <- lcp_cal$N / sum(lcp_cal$N)
lcp_cal$Ncum <- cumsum(lcp_cal$Nperc)

lcp_ext <- lcp_fun(ext)
lcp_ext$Nperc <- lcp_ext$N / sum(lcp_ext$N)
lcp_ext$Ncum <- cumsum(lcp_ext$Nperc)

# Why differences between CAL and EXT?
cal$M <- month(cal$MONTH)
d1 <- table(cal$FUEL, cal$M)
apply(d1[, 6:9], 1, sum) / apply(d1, 1, sum)
ext$M <- month(ext$MONTH)
d2 <- table(ext$FUEL, ext$M)
apply(d2[, 7:8], 1, sum) / apply(d2, 1, sum)

quantile(cal$obs, prob = c(.05, .95))
quantile(ext$obs, prob = c(.05, .95))
db$CALEXT <- ifelse(db$YEAR < 2015, "cal", "ext")
boxplot(db$LFMC ~ db$CALEXT)

# Calculate RMSE by fuel type and month
lcm_fun <- function(df) {
    tmp <- lapply(c(unique(df$FUEL), "0_All"), function(lc) {
        if (lc == "0_All") {
            df2 <- df
        } else {
            df2 <- df[df$FUEL == lc, ]
        }
        out <- lapply(1:13, function(mth) {
            if (mth == 13) {
                df_i <- df2
            } else {
                df_i <- df2[df2$MONTH == mth, ]   
            }
            data.frame(
                "LC" = lc,
                "MONTH" = mth,
                "RMSE" = with(df_i, RMSE(y_true = obs, y_pred = pred)),
                "MBE" = with(df_i, MBE(obs = obs, pred = pred)),
                "meanObs" = mean(df_i$obs),
                "sdObs" = sd(df_i$obs),
                "N" = nrow(df_i)
            )
        })
        do.call(rbind, out)
    })
    out2 <- do.call(rbind, tmp)
    out2[order(out2$LC), ]
}
lcm_cal <- lcm_fun(cal)
lcm_ext <- lcm_fun(ext)

# Plot settings
lbls <- str_replace(sort(unique(lcm_cal$LC)), "[:digit:]_", "")
lbls2 <- str_sub(lbls, 1, 3)
lbls3 <- sprintf("%s (%s)", lbls, lbls2)[-1]
lbls_m <- month.abb
lbls_m[seq(2, 12, 2)] <- ""
lbls_m2 <- c(lbls_m, "All")
cols <- c("#800000", viridis::viridis(5)[2], viridis::viridis(5)[5:4])

# RMSE heatmap plot function
plot_hm <- function(df0, tag, variable, ylabs, xlabs) {
    df <- df0[, c("MONTH", "LC", variable)]
    names(df)[3] <- "var1"
    p <- ggplot(data = df,
                aes(x = as.factor(MONTH), y = LC, fill = var1)) +
        geom_tile() +
        geom_hline(yintercept = 1.5,
                   colour = "white",
                   size = lw_m) +
        scale_y_discrete(expand = c(0, 0), labels = ylabs) +
        scale_x_discrete(expand = c(0, 0), labels = xlabs) +
        labs(tag = tag, x = NULL, y = NULL) +
        custom_theme +
        theme(
            plot.tag = element_text(size = 10,
                                    vjust = -3,
                                    margin = margin(0, 0, 0, 2)),
            plot.tag.position = "topright"
        )
    }

# RMSE heatmaps
p1 <- plot_hm(lcm_cal[lcm_cal$MONTH != 13,], "(c)", "RMSE", lbls2, lbls_m) +
    scale_fill_viridis_c(
        limits = c(0, 40)
    ) +
    theme(
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "mm")
    )

p3 <- plot_hm(lcm_ext[lcm_ext$MONTH != 13,], "(d)", "RMSE", lbls2, lbls_m) +
    scale_fill_viridis_c(
        limits = c(0, 40),
        guide = guide_colorbar(
            title = "RMSE [%]",
            barwidth = unit(2, "mm"),
            barheight = unit(25, "mm"),
            frame.linewidth = unit(.4, "mm"),
            ticks.linewidth = unit(.4, "mm"),
            frame.colour = "black",
            ticks.colour = "black",
            ticks = TRUE,
            draw.ulim = FALSE,
            draw.llim = FALSE,
            default.unit = "mm"
        )
    ) +
    theme(
        legend.position = c(1, 0),
        legend.justification = c(0, 0),
        legend.margin = margin(0, 0, 0, 1, "mm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 9,
                                   margin = margin(0, 0, 0, 0, "mm")),
        legend.title = element_text(size = 9),
        plot.margin = margin(0, 10, 0, 0, "mm")
    )

# Number of testing samples plot function
plot_nsam <- function(df0, tag) {
    ggplot(data = df0, aes(x = as.factor(MONTH))) +
        geom_bar(aes(fill = FUEL),
                 width = .7) +
        scale_fill_manual(
            labels = lbls3,
            values = cols
        ) +
        scale_y_continuous(expand = c(.01, 10)) +
        scale_x_discrete(expand = c(.05, .05),
                         labels = lbls_m) +
        labs(fill = "Fuel Type", x = NULL) +
        labs(tag = tag) +
        custom_theme +
        theme(
            legend.justification = c(-.05, 1.05),
            legend.margin = margin(.5, .5, .5, .5, "mm"),
            legend.text = element_text(size = 7),
            legend.title = element_text(size = 7),
            legend.key.size = unit(3, "mm"),
            panel.grid = element_blank(),
            plot.tag = element_text(size = 10,
                                    vjust = -3,
                                    margin = margin(0, 0, 0, 0)),
            plot.tag.position = "topright"
            
        )
}
p2 <- plot_nsam(cal, "(a)") + ylab("LFMC samples") + theme(legend.position = c(0, 1))
p4 <- plot_nsam(ext, "(b)") + ylab(NULL) + theme(legend.position = "none")

# Save plot
fig <- p2 + p4 + p1 + p3 + plot_layout(ncol = 2)
ggsave(plot = fig,
       filename = "fig5.pdf",
       device = "pdf",
       path = "./output",
       width = 160,
       height = 90,
       units = "mm")

# Land Cover classification anomalies
table(db[, c("LAND_COVER", "FUEL")])
table(db[, c("LAND_COVER", "IGBP_NAME")])
table(db[db$LFMC < 30, "FUEL"])
table(db[db$LFMC > 120, "FUEL"])

# LFMC sample variability (DATA DESCRIPTION)
p1 <- plot_hm(lcm_cal, "(a)", "meanObs", lbls2, lbls_m2) +
    geom_vline(xintercept = 12.5,
               colour = "white",
               size = lw_m) +
    scale_fill_gradient2(
        limits = c(55, 130),
        low = "#d73027",
        mid = "#ffffbf",
        high = "#4575b4",
        midpoint = round(mean(db$LFMC), 0)
    ) +
    theme(
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "mm")
    )

p2 <- plot_hm(lcm_ext, "(b)", "meanObs", lbls2, lbls_m2) +
    geom_vline(xintercept = 12.5,
               colour = "white",
               size = lw_m) +
    scale_fill_gradient2(
        limits = c(55, 130),
        low = "#d73027",
        mid = "#ffffbf",
        high = "#4575b4",
        midpoint = round(mean(db$LFMC), 0),
        guide = guide_colorbar(
            title = "LFMC [%]",
            barwidth = unit(2, "mm"),
            barheight = unit(25, "mm"),
            frame.linewidth = unit(.4, "mm"),
            ticks.linewidth = unit(.4, "mm"),
            frame.colour = "black",
            ticks.colour = "black",
            ticks = TRUE,
            draw.ulim = FALSE,
            draw.llim = FALSE,
            default.unit = "mm"
        )
    ) +
    theme(
        legend.position = c(1, 0),
        legend.justification = c(0, 0),
        legend.margin = margin(0, 0, 0, 1, "mm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 9,
                                   margin = margin(0, 0, 0, 0, "mm")),
        legend.title = element_text(size = 9),
        plot.margin = margin(0, 10, 0, 0, "mm")
    )

p3 <- plot_hm(lcm_cal, "(c)", "sdObs", lbls2, lbls_m2) +
    geom_vline(xintercept = 12.5,
               colour = "white",
               size = lw_m) +
    scale_fill_gradient(
        limits = c(2, 45),
        low = "#fff7ec",
        high = "#7f0000",
    ) +
    theme(
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "mm")
    )

p4 <- plot_hm(lcm_ext, "(d)", "sdObs", lbls2, lbls_m2) +
    geom_vline(xintercept = 12.5,
               colour = "white",
               size = lw_m) +
    scale_fill_gradient(
        limits = c(2, 45),
        low = "#fff7ec",
        high = "#7f0000",
        guide = guide_colorbar(
            title = "SD [%]",
            barwidth = unit(2, "mm"),
            barheight = unit(25, "mm"),
            frame.linewidth = unit(.4, "mm"),
            ticks.linewidth = unit(.4, "mm"),
            frame.colour = "black",
            ticks.colour = "black",
            ticks = TRUE,
            draw.ulim = FALSE,
            draw.llim = FALSE,
            default.unit = "mm"
        )
    ) +
    theme(
        legend.position = c(1, 0),
        legend.justification = c(0, 0),
        legend.margin = margin(0, 0, 0, 1, "mm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 9,
                                   margin = margin(0, 0, 0, 0, "mm")),
        legend.title = element_text(size = 9),
        plot.margin = margin(0, 10, 0, 0, "mm")
    )

# Save plot
fig <- p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
ggsave(plot = fig,
       filename = "figS2.pdf",
       device = "pdf",
       path = "./output",
       width = 170,
       height = 90,
       units = "mm")

# Num. samples LFMC >120%
vc <- unique(db$FUEL)
d1 <- do.call(c, lapply(vc, function(x) {sum(db[db$FUEL == x, "LFMC"] > 120)}))
d2 <- do.call(c, lapply(vc, function(x) sum(db$FUEL == x)))
round((d1/d2) * 100, 1)

# PREDICTIONS IN TIME -----------------------------------------------------
# (Not included in the article)

# Load data
ext <- readRDS("data/calext/nofilt.rds")$val
ext <- merge(ext, db[, c("ID", "DATE", "IGBP_NAME", "FUEL", "SITE")])
ext$MONTH <- round_date(ext$DATE, "month")
ext$LC <- str_replace(ext$FUEL, "[:digit:]_", "")

# Calculate RMSE at site scale
site_err <- lapply(unique(ext$SITE), function(site) {
    df <- ext[ext$SITE == site, ]
    data.frame(
        "SITE" = site,
        "IGBP_NAME" = modeest::mfv1(df$IGBP_NAME),
        "LC" = modeest::mfv1(df$LC),
        "N" = nrow(df),
        "RMSE" = RMSE(y_pred = df$pred, y_true = df$obs),
        "bias" = MBE(obs = df$obs, pred = df$pred)
    )
})
site_err <- do.call(rbind, site_err)
site_err <- site_err[order(site_err$RMSE), ]

# Are errors related to number of samples or land cover?
cor.test(site_err$N, site_err$RMSE)
summary(aov(RMSE ~ LC, site_err))

# Calculate RMSE quantiles
q_err <- quantile(site_err$RMSE, prob = c(.05, .5, .95))

# Get sites with RMSE nearest to error quantiles and N >30 samples
site_err$q05 <- site_err$RMSE - q_err[1]
site_err$q5 <- site_err$RMSE - q_err[2]
site_err$q95 <- site_err$RMSE - q_err[3]
site_err <- site_err[site_err$N > 30, ]
s1 <- site_err[which.min(abs(site_err$q05)), ]
s2 <- site_err[which.min(abs(site_err$q5)), ]
s3 <- site_err[which.min(abs(site_err$q95)), ]
rbind(s1, s2, s3)

# Plot function
plot_ts <- function(site,
                    tag,
                    xlim = NULL,
                    legend_pos = "none",
                    plt_yrs = FALSE) {
    df <- ext[ext$SITE == site,]
    if (is.null(xlim)) {
        xlim <- c(min(df$MONTH), max(df$DATE))
    }
    # Plot settings
    tseq <- seq(xlim[1], xlim[2], "months")
    xlabels <- tseq[month(tseq) %in% c(1, 4, 7, 10)]
    xgrid <- tseq[month(tseq) == 1]
    xlabels_sec <- tseq[month(tseq) == 7]
    if (plt_yrs) {
        sec_x <- element_text(margin = margin(0, 0, .2, 0, "mm"))
    } else {
        sec_x <- element_blank()
    }
    rmse <- RMSE(y_true = df$obs, y_pred = df$pred)
    bias <- MBE(obs = df$obs, pred = df$pred)
    mtext <- sprintf("RMSE: %.1f%%; MBE: %.1f%%", rmse, bias)
    # Plot
    ggplot(data = df, aes(x = DATE)) +
        geom_line(aes(y = obs, colour = "Measured LFMC"), size = lw_l) +
        geom_line(aes(y = pred, colour = "Predicted LFMC"), size = lw_l) +
        stat_difference(aes(ymin = pred, ymax = obs), alpha = 0.3) +
        geom_point(aes(y = obs, colour = "Measured LFMC"), shape = 20) +
        geom_point(aes(y = pred, colour = "Predicted LFMC"), shape = 20) +
        geom_vline(
            xintercept = xgrid[-1],
            linetype = "dashed",
            size = lw_m,
            colour = "gray40"
        ) +
        scale_colour_manual(name = NULL,
                            values = c("#43A2CA", "#CA6B43")) +
        scale_fill_manual(
            name = NULL,
            labels = c("Underpredict", "Overpredict"),
            values = c("#43A2CA", "#CA6B43")
        ) +
        scale_x_date(
            name = NULL,
            breaks = xlabels,
            labels = month.abb[month(xlabels)],
            limits = xlim,
            sec.axis = dup_axis(breaks = xlabels_sec,
                                labels = year(xlabels_sec)),
            expand = c(.01, .01)
        ) +
        scale_y_continuous(
            name = "LFMC [%]",
            limits = c(30, 180),
            breaks = seq(40, 160, 40),
            expand = c(0, 0)
        ) +
        annotate(
            "text",
            x = xlim[1],
            y = Inf,
            label = tag,
            size = 10 / .pt,
            vjust = 1.4,
            hjust = .2
        ) +
        annotate(
            "label",
            x = xlim[1] + 45,
            y = Inf,
            fill = "white",
            label.padding = unit(1, "mm"),
            label.size = NA,
            label = mtext,
            size = 9 / .pt,
            vjust = 1.1,
            hjust = 0
        ) +
        custom_theme +
        theme(
            panel.grid = element_blank(),
            axis.ticks.x.top = element_blank(),
            axis.text.x.top = sec_x,
            legend.title = element_blank(),
            legend.position = legend_pos,
            legend.justification = c(1.01, 1.05),
            legend.text = element_text(size = 9),
            legend.key.size = unit(4, "mm"),
            legend.direction = "horizontal",
            legend.box = "horizontal",
            plot.margin = margin(0, .5, 0, 0, "mm")
        )
}

# Plot time-series
p1 <- plot_ts(
    site = s1$SITE,
    tag = "(a)",
    xlim = c(min(ext$MONTH), max(ext$DATE)),
    legend_pos = c(1, 1),
    plt_yrs = TRUE
)
p2 <- plot_ts(
    site = s2$SITE,
    tag = "(b)",
    xlim = c(min(ext$MONTH), max(ext$DATE))
)
p3 <- plot_ts(
    site = s3$SITE,
    tag = "(c)",
    xlim = c(min(ext$MONTH), max(ext$DATE)),
)
fig <- p1 + p2 + p3 + plot_layout(ncol = 1)
# ggsave(plot = fig,
#        filename = "fig6.pdf",
#        device = "pdf",
#        path = "./output",
#        width = 190,
#        height = 120,
#        units = "mm")

# Additional analysis
df <- ext[ext$SITE == s3$SITE,]
sum(df$pred - df$obs > 0) / nrow(df)  #overpredict
sum(df$pred - df$obs < 0) / nrow(df)  #underpredict
aggregate(site_err$bias, list(site_err$LC), mean)

# EXPLORE DIFFERENCES BETWEEN SITES FROM PREVIOUS ANALYSIS STEP -----------

# Subset sites
sites <- c(s1$SITE, s2$SITE, s3$SITE)
df <- db[db$SITE %in% sites, ]
df$SITE <- factor(df$SITE, levels = sites, ordered = TRUE)

# Species collected and landcover
lapply(sites, function(site) {
    d <- df[df$SITE == site, c("SITE", "SPECIES_COLLECTED", "IGBP_NAME")]
    d <- d[!duplicated(d$SPECIES_COLLECTED), ]
})

# Exploration of some variables
levels(df$SITE) <- c("p5", "p50", "p95")
p1 <- ggplot(df, aes(x = SITE, y = LFMC)) +
    geom_boxplot(fill = "grey90",
                 size = lw_s,
                 outlier.size = 1,
                 outlier.shape = 20) +
    annotate(
        "text",
        x = -Inf,
        y = Inf,
        label = "(a)",
        size = 10 / .pt,
        vjust = 1.4,
        hjust = -.5
    ) +
    labs(x = NULL, y = "LFMC [%]") +
    custom_theme +
    descr_theme
p2 <- ggplot(df, aes(x = SITE, y = NDVIcv)) +
    geom_boxplot(fill = "grey90",
                 size = lw_s,
                 outlier.size = 1,
                 outlier.shape = 20) +
    annotate(
        "text",
        x = -Inf,
        y = Inf,
        label = "(b)",
        size = 10 / .pt,
        vjust = 1.4,
        hjust = -.5
    ) +
    labs(x = NULL, y = "NDVIcv") +
    custom_theme +
    descr_theme
fig <- p1 | p2
ggsave(plot = fig,
       filename = "figS7.pdf",
       device = "pdf",
       path = "./output",
       width = 146,
       height = 65,
       units = "mm")

by(df$LFMC, list(df$SITE), median)
by(df$LFMC, list(df$SITE), range)
by(df$NDVIcv, list(df$SITE), median, na.rm = TRUE)
by(df$NDVIcv, list(df$SITE), range, na.rm = TRUE)

# MARGINAL EFFECTS --------------------------------------------------------

# Calibrate RF model
hyp <- readRDS("data/calext/nofilt.rds")$hyp
vs <- readRDS("data/vs/varsel_final.rds")$predictors
vs <- vs[[length(vs)]]
db$LST8day <- db$LST8day - 273.15
db2 <- db[db$YEAR < 2015, c("LFMC", vs)]
set.seed(42)
mdl <- randomForest::randomForest(
    formula = LFMC ~ .,
    data = db2,
    ntree = hyp$num_trees,
    mtry = hyp$mtry,
    nodesize = hyp$min_node_size,
    replace = FALSE
)

# Partial dependence plots
pdf(file = "output/fig6.pdf", width = 6, height = 6)
par(mfrow = c(3, 3), mar = c(3, 3, .5, .5), mgp = c(2, .7, 0))
for (i in 1:7) {
    randomForest::partialPlot(
        mdl,
        pred.data = db[, c("LFMC", vs)],
        x.var = vs[i],
        col = GnBu(4)[4],
        lwd = 2,
        main = NULL,
        xlab = ifelse(vs[i] == "LST8day", "LST [°C]", vs[i]),
        ylab = "LFMC [%]",
        ylim = c(70, 115)
    )
}
par(op)
dev.off()

# LST vs DOY
{
    par(mfrow = c(2, 4), mar = c(4, 4, 1, .5), pch = 16)
    plot(db$DOY_COS, db$LST8day)
    plot(db$DOY_COS, db$LFMC)
    plot(db$MONTH, db$DOY_COS)
    plot(db$MONTH, db$LST8day)
    plot(db$DOY_SIN, db$LFMC)
    plot(db$MONTH, db$DOY_SIN)
    #plot(db$DOY_SIN, db$LST8day)
    plot(db$LST8day, db$LFMC)
    par(op)   
}
quantile(db$LST8day, prob = c(0, .1))

# NR5
cal <- merge(cal, db[, c("ID", vs)], by = "ID")
cal$dif <- cal$pred - cal$obs
summarize(cal[cal$NR5 > .3, ])
with(cal[cal$NR5 > .3, ], hist(dif))
with(cal[cal$NR5 > .3, ], hist(MONTH))
as.data.frame(with(cal[cal$NR5 > .3, ], table(LC_NAME)))

# SAMPLE PLOTS ------------------------------------------------------------

# Sites
sites <- st_as_sf(
    db[, c("SITE", "LON", "LAT")],
    coords = c("LON", "LAT"),
    crs = "epsg:4326"
)
sites <- aggregate(sites, list(sites$SITE), function(x) length(x))
sites$int <- cut(sites$SITE,
                 breaks = c(0, 50, 100, 150, 200, 300),
                 labels = c("1-50", "50-100", "100-150", "150-200", "200-300"))

# Biomes
biomes <- st_read("data/shp/biomes.shp")
sites$biome <- extract(as_Spatial(biomes), as_Spatial(sites))$BIOME_NAME
d <- aggregate(sites$SITE, list(sites$biome), sum)
d$perc <- d$x / sum(d$x)

# Preparing data for plotting
biomes[!biomes$BIOME_NUM %in% c(4:5, 12), "BIOME_NAME"] <- "Others"
biomes[!biomes$BIOME_NUM %in% c(4:5, 12), "BIOME_NUM"] <- 99
biomes[!biomes$BIOME_NUM %in% c(4:5, 12), "COLOR_BIO"] <- "#ebebeb"
biomes2 <- aggregate(biomes, list(biomes$BIOME_NAME), max)
biomes2 <- biomes2[order(biomes2$BIOME_NUM), ]

# Location plot
fig <- ggplot(data = sites) +
    geom_sf(data = biomes2,
            aes(fill = as.factor(BIOME_NUM)),
            colour = "gray60",
            size = lw_s) +
    geom_sf(aes(size = as.numeric(int)),
            colour = alpha("black", .7),
            shape = 16) +
    coord_sf(
        xlim = c(-19, 43),
        ylim = c(26.5, 46.5),
        expand = FALSE
    ) +
    xlab(NULL) + ylab(NULL) +
    scale_size(name = "Number of Samples",
               labels = levels(sites$int),
               range = c(.5, 3),
               guide = guide_legend(title.position = "top", order = 1)) +
    scale_fill_manual(
        name = "Biome",
        values = c("#AAD1AC", "#daeed1", "#FEDC6D", "#ebebeb"),
        labels = biomes2$BIOME_NAME,
        guide = guide_legend(ncol = 2, title.position = "top", order = 2)
    ) +
    custom_theme +
    theme(
        axis.text.y = element_text(angle = 90, hjust = .5),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7, face = "bold"),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(1, 1, 1, 1, "mm"),
        legend.spacing = unit(4, "mm"),
        plot.margin = margin(1, 1, 0, 1, "mm")
    )
fig

# Save plot
ggsave(plot = fig,
       filename = "fig2.pdf",
       device = "pdf",
       path = "./output",
       width = 192,
       units = "mm")

# END ---




