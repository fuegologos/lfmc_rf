#**************************************************************************
#* LFMC - A remote sensing approach to predict LFMC at large scale
#* Exploratory data analysis
#* Author: Angel Cunill Camprubi <acunill86@gmail.com>
#* Last update: 20 Sep 2021
#**************************************************************************

# SETTINGS ----------------------------------------------------------------

# Packages
library(stringr)
library(lubridate)
library(grid)
library(ggplot2)
library(ggpointdensity)
library(patchwork)
library(RColorBrewer)
library(MASS)                  # True histogram function
library(PerformanceAnalytics)  # Correlogram and moments
library(corrplot)              # Correlogram
library(DescTools)
library(MLmetrics)

# Functions
MBE <- function(obs, pred) {
    mean(obs - pred, na.rm = TRUE)
}
ubRMSE <- function(obs, pred) {
    m_obs <- mean(obs)
    m_pred <- mean(pred)
    sqrt(mean((obs - pred - (m_obs - m_pred)) ^ 2))
}

# Others
options(stringsAsFactors = FALSE)

# GRAPHICAL PARAMETERS ----------------------------------------------------

op <- par(no.readonly = TRUE)
convertUnit(unit(.28, "points"), "mm", valueOnly = TRUE)
convertUnit(unit(60, "mm"), "inches", valueOnly = TRUE)
convertUnit(unit(190, "mm"), "inches", valueOnly = TRUE)

#* Line width (mm_real * 1/0.75 = mm_R) ----
.cf <- 1/.75
lw_s <- .18 * .cf
lw_m <- .25 * .cf
lw_l <- .35 * .cf
lw_xl <- .6 * .cf

#* Colour palette ----
YlGnBu <- colorRampPalette(brewer.pal(3, "YlGnBu"))
GnBu <- colorRampPalette(brewer.pal(3, "GnBu"))
red_line <- "#de2d26"

#* Custom theme ----

# Custom theme (lfmc)
custom_theme <-
    theme_bw(
        base_size = 8,
        base_line_size = lw_m,
        base_rect_size = 0
    ) +
    theme(
        axis.ticks.length = unit(1, "mm"),
        axis.text = element_text(size = 7),
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
db <- readRDS("./data/tmp/db03_lfmc.rds")

# Features collection
si <- c('NDVI', 'NDII6', 'NDII7', 'GVMI', 'NDWI', 'EVI', 'SAVI', 'VARI', 'VIgreen', 'NDTI', 'MSI', 'Gratio')
bands <- c('NR1', 'NR2', 'NR3', 'NR4', 'NR5', 'NR6', 'NR7')
vars <- c('HYBFOR', si, bands)
mdl_vars <- c("LFMC", vars)
names(db)[names(db) == "LFMC_AVG"] <- "LFMC"

# EXPLORATORY DATA ANALYSIS -----------------------------------------------

# Summaries
summary(db$LFMC)
quantile(db$LFMC, probs = c(.05, .95))
sd(db$LFMC)
skewness(db$LFMC) # positive skew
kurtosis(db$LFMC) # leptokurtic (3 = mesokurtic)

# Number of records
nrow(db)
n_sp <- nrow(db[db$COUNTRY == "Spain", ])
n_fr <- nrow(db[db$COUNTRY == "France", ])
(n_sp + n_fr) / nrow(db)

# Missing data due to NDVIcv filter
n1 <- sum(complete.cases(db$NDVIcv))
nrow(db) - n1
(nrow(db) - n1) / nrow(db)

# Time period
summary(db$DATE)
{
    par(mfrow = c(2, 2))
    for(i in unique(db$COUNTRY)) {
        db_c <- db[db$COUNTRY == i, ]
        db_c$ym <- ym(paste0(db_c$YEAR, "-", db_c$MONTH))
        tmp <- aggregate(db_c$LFMC, list(db_c$ym), length)
        plot(tmp$Group.1, tmp$x, main = i)
        arrows(x0 = tmp$Group.1, x1 = tmp$Group.1, y0 = 0, y1 = tmp$x, angle = 0)
    }
    par(op)
}
tapply(db$LFMC, list(db$MONTH, db$YEAR, db$COUNTRY), length)

# Sites
length(unique(db$SITE))

# Records and sites >2010
sum(db$YEAR > 2010)/ nrow(db) # 70/30 %
length(unique(db[db$YEAR > 2010, "SITE"]))
length(unique(db[db$YEAR <= 2010, "SITE"]))
sum(unique(db[db$YEAR <= 2010, "SITE"]) %in% unique(db[db$YEAR > 2010, "SITE"]))

# Plot settings
descr_theme <- theme(
    axis.title.x.bottom = element_text(margin = margin(-2, 0, 0, 0, "mm")),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "gray90", size = lw_s),
    panel.grid.minor.y = element_line(colour = "gray90", size = lw_s),
    plot.margin = margin(0, 1, 0, 1, "mm")
)

#* LFMC univariate distribution ----
truehist(db$LFMC, h = 1)
lfmc_hist <- ggplot(data = db, aes(x = LFMC)) +
    geom_histogram(aes(y = ..density..),
                   binwidth = 5,
                   fill = GnBu(4)[3],
                   color = gray(.95),
                   size = .04 * .cf) +
    geom_vline(xintercept = 20,
               linetype = "longdash",
               colour = red_line,
               size = lw_m) +
    scale_x_continuous(name = "LFMC (%)",
                       breaks = seq(0, 400, 50),
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

#* LFMC by country ----
n <- as.data.frame(sort(table(db$COUNTRY)))
names(n)[1] <- "count_ord"
db$count_ord <- as.factor(n$Freq[match(db$COUNTRY, levels(n$count_ord))])
xlabs <- paste0(n$count_ord, "\n(", n$Freq, ")")
lfmc_violin <- ggplot(data = db,
                      aes(x = count_ord,
                          y = LFMC,
                          fill = count_ord),
                      color = "grey20") +
    geom_violin(draw_quantiles = .5,
                size = lw_s) +
    geom_hline(yintercept = 20,
               linetype = "longdash",
               colour = red_line,
               size = lw_m) +
    labs(x = NULL, y = "LFMC (%)") +
    annotate("text",
             x = Inf,
             y = Inf,
             label = "(b)",
             size = 10 / .pt,
             vjust = 1.5,
             hjust = 1.3) +
    scale_y_continuous(expand = c(.04, .04), limits = c(0, 400)) +
    scale_x_discrete(labels = xlabs) +
    scale_fill_brewer(type = "seq", palette = "GnBu") +
    custom_theme +
    descr_theme +
    theme(
        legend.position = "none",
    )

#* LFMC time-series ----
db2 <- db[, c("count_ord", "MONTH", "LFMC")]
ts <- as.data.frame(reshape2::melt(table(db$count_ord, db$MONTH)))
ts <- ts[ts$value == 0, ]
ts$value <- -1
names(ts) <- c("count_ord", "MONTH", "LFMC")
db2 <- rbind(db2, ts)

lfmc_ts <- ggplot(db2, aes(x = as.factor(MONTH), y = LFMC, fill = count_ord)) +
    geom_boxplot(outlier.shape = NA, lwd = lw_s) +
    geom_vline(xintercept = seq(1.5, 11.5, 1),
               colour = "gray90",
               size = lw_s) +
    scale_x_discrete(labels = month.abb) +
    coord_cartesian(ylim = c(0, 200), expand = FALSE) +
    scale_fill_brewer(type = "seq",
                      palette = "GnBu",
                      labels = n$count_ord) +
    guides(fill = guide_legend(title = NULL,
                               label.position = "right",
                               direction = "horizontal")) +
    labs(x = NULL, y = "LFMC (%)") +
    annotate("text",
             x = Inf,
             y = Inf,
             label = "(c)",
             size = 10 / .pt,
             vjust = 1.5,
             hjust = 1.3) +
    custom_theme +
    theme(
        panel.grid = element_blank(),
        legend.position = c(.01, .01),
        legend.justification = c(0, 0),
        legend.text = element_text(size = 6),
        legend.margin = margin(0, 0, 0, 0, "mm"),
        legend.key.size = unit(4, "mm"),
        plot.margin = margin(2, 0, 0, 0, "mm")
    )

# Output fig
fig <- (lfmc_hist | lfmc_violin) / lfmc_ts
ggsave(plot = fig,
       filename = "fig2.pdf",
       device = "pdf",
       path = "./output/figs",
       width = 146,
       height = 126,
       units = "mm")
remove(ts, db2)

# CORRELATION PLOTS -------------------------------------------------------

#* Correlogram with scatterplots ----
png("./output/figs/00-correlogram.png", width = 28, height = 20, units = "cm", res = 150)
chart.Correlation(db[, mdl_vars], histogram = FALSE, pch = 16)
dev.off()

#* Correlation matrix plot ----
lfmc_corr <- round(cor(db[, mdl_vars]), 2)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
pdf("./output/figs/figS1.pdf", width = 6, height = 6)
corrplot(
    corr = lfmc_corr,
    method = "color",
    type = "upper",
    col = col(200),
    diag = FALSE,
    addgrid.col = "grey",
    addCoef.col = "black",
    tl.cex = .7,
    tl.col = "black",
    cl.cex = .7,
    cl.align.text = "l",
    number.cex = .5
)
dev.off()

# VARIABLE IMPORTANCE -----------------------------------------------------

# Plot settings
vi_theme <- theme(
    panel.grid.major.x = element_line(linetype = "dashed",
                                      colour = "gray90",
                                      size = lw_s),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x.top = element_text(margin = margin(0, 0, .6, 0, "mm")),
    axis.title.x.top = element_text(margin = margin(0, 0, 1, 0, "mm")),
    plot.margin = margin(0, 1, 0, 1, "mm")
)

#* VI Model performance ----
mp_varsel <- readRDS("data/mp_varsel.rds")
mp_vi <- data.frame("id" = seq_along(vars))
tmp <- lapply(1:5, function(i) {
    df <- data.frame(
        "id" = mp_varsel[[i]]$imp.mean.dec.ind,
        "vi" = mp_varsel[[i]]$imp.mean.dec
    )
    names(df)[2] <- paste0("vi", i)
    df[order(df$id), 2, drop = FALSE]
})
mp_vi <- cbind(mp_vi, do.call(cbind, tmp))
mp_vi$avg <- apply(mp_vi[, 2:6], 1, mean)
mp_vi$min <- apply(mp_vi[, 2:6], 1, min)
mp_vi$max <- apply(mp_vi[, 2:6], 1, max)
mp_vi$name <- vars[mp_vi$id]

mp_plot <- ggplot(mp_vi, aes(x = avg, y = reorder(name, avg))) +
    geom_segment(aes(yend = name), xend = 100, colour = "grey80", size = lw_m) +
    geom_segment(aes(yend = name, x = min, xend = max),
                 colour = "grey70",
                 size = lw_l) +
    geom_point(aes(x = min), size = .6, colour = "grey60") +
    geom_point(aes(x = max), size = .6, colour = "grey60") +
    geom_point(size = .8) +
    labs(y = NULL, x = "MSE") +
    annotate("text",
             x = Inf,
             y = Inf,
             label = "(a)",
             size = 10 / .pt,
             vjust = 1.5,
             hjust = 1.3) +
    scale_x_continuous(expand = c(0, 0),
                       limits = c(100, 1100),
                       sec.axis = sec_axis(trans = ~ . ^ .5,
                                           name = "RMSE (%)")) +
    custom_theme +
    vi_theme

#* VI 3-fold CV ----
cv_varsel <- readRDS("data/3fcv_varsel.rds")
cv_vi <- data.frame(
    "id" = cv_varsel$imp.mean.dec.ind,
    "name" = vars[cv_varsel$imp.mean.dec.ind],
    "vi" = cv_varsel$imp.mean.dec,
    "visd" = cv_varsel$imp.sd.dec
)
cv_vi$minsd <- cv_vi$vi - cv_vi$visd
cv_vi$maxsd <- cv_vi$vi + cv_vi$visd

cv_plot <- ggplot(cv_vi, aes(x = vi, y = reorder(name, vi))) +
    geom_segment(aes(yend = name), xend = 100, colour = "grey80", size = lw_m) +
    geom_segment(aes(yend = name, x = minsd, xend = maxsd),
                 colour = "grey40",
                 size = lw_m,
                 arrow = arrow(ends = "both",
                               angle = 90,
                               length = unit(.5, "mm"))) +
    geom_point(size = .8) +
    labs(y = NULL, x = "MSE") +
    annotate("text",
             x = Inf,
             y = Inf,
             label = "(b)",
             size = 10 / .pt,
             vjust = 1.5,
             hjust = 1.3) +
    scale_x_continuous(expand = c(0, 0),
                       limits = c(100, 1100),
                       sec.axis = sec_axis(trans = ~ . ^ .5,
                                           name = "RMSE (%)")) +
    custom_theme +
    vi_theme

#* VI Extrapolation ----
ext_varsel <- readRDS("data/ext_varsel.rds")
ext_vi <- data.frame(
    "id" = ext_varsel$imp.mean.dec.ind,
    "name" = vars[ext_varsel$imp.mean.dec.ind],
    "vi" = ext_varsel$imp.mean.dec,
    "visd" = ext_varsel$imp.sd.dec
)
ext_vi$minsd <- ext_vi$vi - ext_vi$visd
ext_vi$maxsd <- ext_vi$vi + ext_vi$visd

ext_plot <- ggplot(ext_vi, aes(x = vi, y = reorder(name, vi))) +
    geom_segment(aes(yend = name), xend = 0, colour = "grey80", size = lw_m) +
    geom_segment(aes(yend = name, x = minsd, xend = maxsd),
                 colour = "grey40",
                 size = lw_m,
                 arrow = arrow(ends = "both",
                               angle = 90,
                               length = unit(.5, "mm"))) +
    geom_point(size = .8) +
    labs(y = NULL, x = "MSE") +
    annotate("text",
             x = Inf,
             y = Inf,
             label = "(c)",
             size = 10 / .pt,
             vjust = 1.5,
             hjust = 1.3) +
    scale_x_continuous(expand = c(0, 0),
                       limits = c(100, 1100),
                       sec.axis = sec_axis(trans = ~ . ^ .5,
                                           name = "RMSE (%)")) +
    custom_theme +
    vi_theme

# Variable importance output fig
fig <- mp_plot + cv_plot + ext_plot + plot_layout(widths = c(1, 1, 1))
ggsave(plot = fig,
       filename = "fig3.pdf",
       device = "pdf",
       path = "./output/figs",
       width = 186,
       height = 92,
       units = "mm")

# VARIABLE SELECTION ------------------------------------------------------

# Plot function
plot_vs <- function(vs_data, plot_tag, xlabels_num = FALSE, ytitle = NULL) {
    vs1 <- data.frame(
        "vs" = vs_data$varselect.thres,
        "vname" = vars[vs_data$varselect.thres],
        "vi" = vs_data$imp.mean.dec,
        "err_int" = vs_data$err.interp,
        "v_int" = ifelse(vs_data$varselect.thres %in% vs_data$varselect.interp, 2, 1)
    )
    vs2 <- data.frame(
        "vs" = vs_data$varselect.pred,
        "err_pred" = vs_data$err.pred
    )
    vs_df <- merge(vs1, vs2, by = "vs", all.x = TRUE)
    vs_df <- vs_df[order(vs_df$vi, decreasing = TRUE), ]
    vs_df$id <- seq_len(nrow(vs_df))
    
    if(xlabels_num) {
        xlabels <- vs_df$vs
    } else {
        xlabels <- vs_df$vname
    }
    ggplot(data = vs_df, aes(x = id, y = err_int)) +
        geom_line(#data = vs_df[vs_df$v_int == 2, ],
                  size = lw_m,
                  colour = "gray70") +
        geom_point(aes(shape = as.factor(v_int)), size = .8, colour = "gray70") +
        geom_line(data = vs_df[complete.cases(vs_df$err_pred), ],
                  aes(x = id, y = err_pred),
                  colour = red_line,
                  size = lw_m) +
        geom_point(data = vs_df[complete.cases(vs_df$err_pred), ],
                   aes(x = id, y = err_pred),
                   colour = red_line,
                   shape = 16,
                   size = .8) +
        labs(x = NULL, y = ytitle) +
        annotate("text",
                 x = Inf,
                 y = Inf,
                 label = plot_tag,
                 size = 10 / .pt,
                 vjust = 1.5,
                 hjust = 1.3) +
        scale_shape_manual(values = c(4, 16)) +
        scale_x_continuous(labels = xlabels,
                           breaks = 1:20,
                           expand = c(.02, .02),
                           ) +
        custom_theme +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            panel.grid = element_blank(),
            plot.margin = margin(2, 2, 0, 0, "mm")
        )    
}

#* VS Model performance ----
mp_varsel <- readRDS("data/mp_varsel.rds")
tmp <- lapply(1:5, function(i) {
    mp_varsel[[i]]$varselect.pred
})
tmp_freq <- do.call(c, tmp)
tmp_freq <- as.data.frame(table(tmp_freq))
mp_vs_freq <- merge(mp_vi, tmp_freq, by.x = "id", by.y = "tmp_freq", all.x = TRUE)
mp_vs_freq[is.na(mp_vs_freq$Freq), "Freq"] <- 0
(mp_vs_freq <- mp_vs_freq[order(mp_vs_freq$Freq, mp_vs_freq$avg, decreasing = TRUE), ])
lapply(1:5, function(i) {
    vars[mp_varsel[[i]]$varselect.pred]
})

mp1 <- plot_vs(mp_varsel[[1]], "F1", ytitle = "OOB error (MSE)")
mp2 <- plot_vs(mp_varsel[[2]], "F2")
mp3 <- plot_vs(mp_varsel[[3]], "F3")
mp4 <- plot_vs(mp_varsel[[4]], "F4", ytitle = "OOB error (MSE)")
mp5 <- plot_vs(mp_varsel[[5]], "F5")

txt1 <- paste(paste0(mp_vs_freq$name[1:10], " (", mp_vs_freq$Freq[1:10], ")"),
              collapse = "\n")
txt2 <- paste(paste0(mp_vs_freq$name[11:20], " (", mp_vs_freq$Freq[11:20], ")"),
              collapse = "\n")

mp_txt1 <- grid::textGrob(txt1, x = 0.3, hjust = 0, gp = gpar(fontsize = 8))
mp_txt2 <- grid::textGrob(txt2, x = 0.1, hjust = 0, gp = gpar(fontsize = 8))
mp_txt <- gridExtra::arrangeGrob(mp_txt1, mp_txt2, ncol = 2)
fig <- mp1 + mp2 + mp3 + mp4 + mp5 + plot_layout(ncol = 3) + mp_txt
ggsave(plot = fig,
       filename = "fig4.pdf",
       device = "pdf",
       path = "./output/figs",
       width = 190,
       height = 114,
       units = "mm")

#* VS 3FCV and EXT plot ----
cv_varsel <- readRDS("data/3fcv_varsel.rds")
vars[cv_varsel$varselect.pred]

p1 <- plot_vs(cv_varsel,
              plot_tag = "(a)",
              xlabels_num = FALSE,
              ytitle = "OOB error (MSE)")

ext_varsel <- readRDS("data/ext_varsel.rds")
vars[ext_varsel$varselect.pred]

p2 <- plot_vs(ext_varsel,
              plot_tag = "(b)",
              xlabels_num = FALSE,
              ytitle = "OOB error (MSE)")

fig <- p1 | p2
ggsave(plot = fig,
       filename = "fig5.pdf",
       device = "pdf",
       path = "./output/figs",
       width = 134,
       height = 55,
       units = "mm")

# DATAFRAME OF RESULTS ----------------------------------------------------

data_df <- data.frame(
    "method" = rep(c("MP", "3FCV", "EXT"), each = 4),
    "vars" = rep(c("allp", "p"), times = 2 * 3),
    "filt" = rep(rep(c("no", "yes"), each = 2), 3),
    "file" = c(
        # Model performance
        "mp_1_allp_nof.rds",
        "mp_2_selp_nof.rds",
        "mp_3_allp_filt.rds",
        "mp_4_selp_filt.rds",
        # 3-fold cv
        "3fcv_1_allp_nof.rds",
        "3fcv_2_selp_nof.rds",
        "3fcv_3_allp_filt.rds",
        "3fcv_4_selp_filt.rds",
        # Extrapolation
        "ext_1_allp_nof.rds",
        "ext_2_selp_nof.rds",
        "ext_3_allp_filt.rds",
        "ext_4_selp_filt.rds"
    ),
    "eprof" = c(
        rep(NA, 4),
        # 3-fold cv
        "3fcv_eprof_allp_nof.rds",
        "3fcv_eprof_selp_nof.rds",
        "3fcv_eprof_allp_filt.rds",
        "3fcv_eprof_selp_filt.rds",
        # Extrapolation
        "ext_eprof_allp_nof.rds",
        "ext_eprof_selp_nof.rds",
        "ext_eprof_allp_filt.rds",
        "ext_eprof_selp_filt.rds"
    ) 
)

# MODELS ACCURACY ---------------------------------------------------------

# Accuracy function
acc_metrics <- function(i, fd_range = NULL) {
    # Load data
    f <- file.path("data", data_df[i, "file"])
    df <- readRDS(f)
    
    # Validation values
    if (i < 9) {
        dfv <- df$validation    
    } else {
        dfv <- df
    }
    names(dfv)[ncol(dfv)] <- "pred"
    
    # Fire danger range
    if(!is.null(fd_range)) {
        dfv <- dfv[dfv$obs > min(fd_range) & dfv$obs < max(fd_range), ]
        fdr <- paste(fd_range, collapse = "-")
    }
    
    # Linear model
    op_lm <- summary(lm(obs ~ pred, data = dfv))$coefficients
    
    # Output table
    data.frame(
        "method" = data_df[i, "method"],
        "vars" = data_df[i, "vars"],
        "filt" = data_df[i, "filt"],
        "fd_range" = ifelse(is.null(fd_range), "null", fdr),
        "bias" = MBE(obs = dfv$obs, pred = dfv$pred),
        "RMSE" = RMSE(y_true = dfv$obs, y_pred = dfv$pred),
        "ubRMSE" = ubRMSE(obs = dfv$obs, pred = dfv$pred),
        "CCC" = CCC(x = dfv$obs, y = dfv$pred)$rho.c[[1]],
        "R2" = R2_Score(y_true = dfv$obs, y_pred = dfv$pred),
        "r" = sqrt(R2_Score(y_true = dfv$obs, y_pred = dfv$pred)),
        "a" = op_lm[1, 1],
        "a_pvalue" = op_lm[1, 4],
        "b" = op_lm[2, 1],
        "b_pvalue" = op_lm[2, 4]
    )
}

#* Accuracy metrics ----
acc <- lapply(seq_len(nrow(data_df)), acc_metrics)
acc_df <- do.call(rbind, acc)
write.csv(acc_df, "output/acc.csv", row.names = FALSE)

#* Accuracy within the fire danger levels ----
fd_acc <- lapply(seq_len(nrow(data_df)), acc_metrics, fd_range = c(50, 170))
fd_acc_df <- do.call(rbind, fd_acc)
write.csv(fd_acc_df, "output/acc_fd.csv", row.names = FALSE)
remove(fd_acc, fd_acc_df)

# HYPERPARAMETERS ---------------------------------------------------------

#* Parameters (MP) ----
mp_hyp <- lapply(1:4, function(i) {
    # Load data
    f <- file.path("data", data_df[i, "file"])
    df <- readRDS(f)
    hyp <- df$parameters

    # Output table
    out1 <- data.frame(
        "method" = data_df[i, "method"],
        "vars" = data_df[i, "vars"],
        "filt" = data_df[i, "filt"],
        "iter" = 1:5
    )
    out2 <- as.data.frame(do.call(rbind, hyp))
    out2$sp <- NA
    if (any(names(out2) == "spf")) {
        out2$sp <- out2$spf
        out2$spf <- NULL
    }
    cbind(out1, out2)
})
mp_hyp_df <- do.call(rbind, mp_hyp)
write.csv(mp_hyp_df, "output/hyp_mp.csv", row.names = FALSE)

#* Parameters (CV - EXT) ----
hyp <- lapply(5:12, function(i) {
    # Load data
    f <- file.path("data", data_df[i, "eprof"])
    df <- readRDS(f)
    hyp_i <- df[which.min(df$rmse), ]
    
    # Output table
    data.frame(
        "id" = i,
        "method" = data_df[i, "method"],
        "vars" = data_df[i, "vars"],
        "filt" = data_df[i, "filt"],
        "num_trees" = hyp_i$num_trees,
        "mtry" = hyp_i$mtry,
        "min_node_size" = hyp_i$min_node_size,
        "sample_fraction" = hyp_i$sample_fraction,
        "sp" = ifelse(any(names(hyp_i) == "sp"), hyp_i$sp, NA),
        "rmse" = hyp_i$rmse
    )
})
hyp_df <- do.call(rbind, hyp)
write.csv(hyp_df, "output/hyp_cvext.csv", row.names = FALSE)

# OBSERVED VS PREDICTED ---------------------------------------------------

# 3-fold CV data
data_df_cv <- data_df[data_df$method == "3FCV", ]
file_id <- which.min(acc_df[acc_df$method == "3FCV", "RMSE"])
cv_df <- readRDS(file.path("data", data_df_cv[file_id, "file"]))$validation

# Extrapolation data
data_df_ext <- data_df[data_df$method == "EXT", ]
file_id <- which.min(acc_df[acc_df$method == "EXT", "RMSE"])
ext_df <- readRDS(file.path("data", data_df_ext[file_id, "file"]))

# RTM data
rtm_df <- db[complete.cases(db$RTM), c("LFMC", "RTM")]
rtm_df <- rtm_df[rtm_df$RTM > 0, ]

# Plot function
plot_op <- function(df, lfmc_obs, lfmc_RF, plotid) {
    form <- formula(paste(lfmc_obs, lfmc_RF, sep = " ~ "))
    cv_lm <- lm(form, data = df)
    cv_lm_pred <- data.frame(c(min(df[, lfmc_RF]), max(df[, lfmc_RF])))
    names(cv_lm_pred) <- lfmc_RF
    cv_lm_pred$y <- predict(cv_lm, cv_lm_pred)
    cv_fun <- sprintf("y=%.2f+%.2fx", cv_lm$coefficients[1], cv_lm$coefficients[2])
    names(df)[names(df) %in% c(lfmc_obs, lfmc_RF)] <- c("obs", "pred")
    ggplot(data = df, aes(x = pred, y = obs)) +
        geom_pointdensity(size = .1, adjust = 1) +
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
        labs(x = "LFMCpred (%)", y = "LFMCobs (%)") +
        annotate("text",
                 x = Inf,
                 y = -Inf,
                 label = cv_fun,
                 size = 8 / .pt,
                 vjust = -.4,
                 hjust = 1.05) +
        annotate("text",
                 x = -Inf,
                 y = Inf,
                 label = plotid,
                 size = 10 / .pt,
                 vjust = 1.4,
                 hjust = -.3) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 300)) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 300)) +
        scale_colour_viridis_c() +
        guides(colour = guide_colourbar(
            title = NULL,
            barwidth = unit(2, "mm"),
            barheight = unit(12, "mm"),
            frame.linewidth = lw_m,
            ticks.linewidth = lw_m,
            draw.ulim = TRUE,
            draw.llim = TRUE,
            default.unit = "mm"
        )) +
        custom_theme +
        theme(
            panel.grid.major.x = element_line(linetype = "dashed",
                                              colour = "gray90",
                                              size = lw_s),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(linetype = "dashed",
                                              colour = "gray90",
                                              size = lw_s),
            panel.grid.minor.y = element_blank(),
            legend.position = c(1, 0),
            legend.justification = c(-.1, 0),
            legend.margin = margin(0, 0, 1, 0, "mm"),
            legend.text = element_text(margin = margin(0, 0, 0, -.5, "mm")),
            legend.background = element_blank(),
            plot.margin = margin(0, 2, 0, 0, "mm")
        )
}

# Output plot
p1 <- plot_op(cv_df, "obs", "RFcv", "(a)")
p2 <- plot_op(ext_df, "obs", "pred", "(b)")
p3 <- plot_op(rtm_df, "LFMC", "RTM", "(c)")
fig <- p1 + p2 + p3 + plot_layout(widths = unit(c(50, 50, 50), "mm"),
                                  heights = unit(50, "mm"))
ggsave(plot = fig,
       filename = "fig6.pdf",
       device = "pdf",
       path = "./output/figs",
       width = 190,
       height = 60,
       units = "mm")

# FUEL TYPE ANALYSIS ------------------------------------------------------

# Land Cover LUT
lc_lut <- read.csv("doc/landcover_fueltype_LUT.csv", sep = ";")

# Load 3FCV error and add to LFMC database
data_df_cv <- data_df[data_df$method == "3FCV", ]
file_id <- which.min(acc_df[acc_df$method == "3FCV", "RMSE"])
cv_df <- readRDS(file.path("data", data_df_cv[file_id, "file"]))$validation
db_ftype <- merge(db[, c("ID", "SITE", "LFMC", "IGBP_LC", "HYBFOR")],
                  cv_df,
                  by = "ID")
identical(db_ftype$LFMC, db_ftype$obs)

# Land cover present
lc_df <- lc_lut[lc_lut$IGBP_ID %in% unique(db_ftype$IGBP_LC), ]
lc_df$freq <- as.data.frame(table(db_ftype$IGBP_LC))$Freq
lc_df

# Anomalous classes
an_sites <- unique(db[db$IGBP_LC == 11, "SITE"])
unique(db[db$SITE %in% an_sites, "IGBP_LC"])

# Fuel type of below 20 and above 250 LFMC
lc_df[lc_df$IGBP_ID %in% unique(db_ftype$IGBP_LC[db_ftype$LFMC < 20]), "IGBP_NAME"]
lc_df[lc_df$IGBP_ID %in% unique(db_ftype$IGBP_LC[db_ftype$LFMC > 250]), "IGBP_NAME"]

# Error by fuel type (IGBP)
tmp <- lapply(c("Forest", "Shrubland", "Grassland", "DNU"), function(x) {
    lc_id <- lc_df[lc_df$FUEL_CLASS == x, "IGBP_ID"]
    df <- db_ftype[db_ftype$IGBP_LC %in% lc_id, ]
    data.frame(
        "fuel_type" = x,
        "n_samples" = nrow(df),
        "n_sites" = length(unique(df$SITE)),
        "sd_obs" = sd(df$obs),
        "bias" = MBE(obs = df$obs, pred = df$RFcv),
        "RMSE" = RMSE(y_true = df$obs, y_pred = df$RFcv),
        "ubRMSE" = ubRMSE(obs = df$obs, pred = df$RFcv),
        "CCC" = CCC(x = df$obs, y = df$RFcv)$rho.c[[1]],
        "R2" = R2_Score(y_true = df$obs, y_pred = df$RFcv),
        "r" = sqrt(R2_Score(y_true = df$obs, y_pred = df$RFcv)) 
    )
})
do.call(rbind, tmp)

# Error by HYBFOR
db_ftype$err <- abs(db_ftype$obs - db_ftype$RFcv)
db_ftype$brk <- cut(db_ftype$HYBFOR,
                    breaks = c(0, 10, 30, 50, 70, max(db_ftype$HYBFOR)),
                    include.lowest = TRUE,
                    right = FALSE)
with(db_ftype, plot(HYBFOR, err, pch = 16))
with(db_ftype, boxplot(err ~ brk,
                       xlab = "Forest Cover (%)",
                       ylab = "Absolute error",
                       main = "Error by HYBFOR"))

# Error distribution
names(cv_df)[4] <- "pred"
cv_df$dif <- cv_df$obs - cv_df$pred

{
    par(mfrow = c(1, 2))
    hist(abs(cv_df$dif), breaks = seq(0, 300, 5), xlab = "absolut.error", main = "Abs.errors distrib. 3Fcv_df (best model)")
    plot(cv_df$obs, cv_df$dif,
         xlab = "Observed LFMC", ylab = "Residual", main = "Diff. (obs-pred) vs obs")
    abline(h = 0, lty = "dashed")
    abline(h = seq(-100, 300, 50), v = seq(0, 400, 50), col = "lightgrey", lty = "dotted")
    par(mfrow = c(1, 1))
}

(q95 <- quantile(cv_df$dif, probs = 0.95))
hist(cv_df$obs[cv_df$dif < q95], breaks = seq(0, 400, 5), xlab = "LFMC (%)", main = "Obs [< q95]")
hist(cv_df$obs[cv_df$dif > q95], breaks = seq(0, 400, 5), xlab = "LFMC (%)", main = "Obs [> q95]")

# END ---
