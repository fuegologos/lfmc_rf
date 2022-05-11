# *************************************************************************
# LFMC-RF - A remote sensing approach to predict LFMC at large scale using
# Random Forests
# Create database for LFMC modelling process
#  - Target variable: LFMC_AVG
#  - Extraction protocol: Simple pixel extraction
# Author: Angel Cunill Camprubi <acunill86@gmail.com>
# Last update: 03 Mar 2022
# *************************************************************************

# SETTINGS ----------------------------------------------------------------

# Packages
library(raster)
library(sf)
library(stringr)
library(lubridate)
library(doParallel)
library(foreach)
source("src/functions/normalise.R")

# Others
options(stringsAsFactors = FALSE)
set.seed(24)

# LOAD DATABASE -----------------------------------------------------------

# Load database
db <- readRDS("./data/db01_field_lfmc.rds")

# Subset Mediterranean countries
db <- db[db$COUNTRY %in% c("Spain", "France", "Italy", "Tunisia"), ]

# DATE-TIME ATTRIBUTES ----------------------------------------------------

# Assign month and year fields
db$MONTH <- month(db$DATE)
db$YEAR <- year(db$DATE)
db$YDOY <- format(db$DATE, "%Y%j")

# Calculate normalised day-of-year
ndoy <- normalise(as.numeric(format(db$DATE, "%j")),
                  in_range = c(1, 366),
                  out_range = c(-pi, pi))
db$DOY_SIN <- round(sin(ndoy), 5)
db$DOY_COS <- round(cos(ndoy), 5)

# SPECTRAL INFORMATION ----------------------------------------------------

# Load data
si <- readRDS("./data/db02_point_nr.rds")
#si <- readRDS("./data/db02_area_nr.rds")

# Greenness related indices
si$NDVI <- with(si, (NR2 - NR1) / (NR2 + NR1))
si$EVI <- with(si, (2.5 * (NR2 - NR1)) / (NR2 + (6 * NR1) - (7.5 * NR3) + 1))
si$SAVI <- with(si, (1 + 0.5) * ((NR2 - NR1) / (NR2 + NR1 + 0.5)))
si$VARI <- with(si, (NR4 - NR1) / (NR4 + NR1 - NR3))
si$VIgreen <- with(si, (NR4 - NR1) / (NR4 + NR1))
si$Gratio <- with(si, NR4 / NR1)

# Water related indices
si$NDII6 <- with(si, (NR2 - NR6) / (NR2 + NR6))
si$NDII7 <- with(si, (NR2 - NR7) / (NR2 + NR7))
si$NDWI <- with(si, (NR2 - NR5) / (NR2 + NR5))
si$GVMI <- with(si, ((NR2 + 0.1) - (NR6 + 0.02)) / ((NR2 + 0.1) + (NR6 + 0.02)))
si$MSI <- with(si, NR6 / NR2)

# Dry-matter related indices
si$NDTI <- with(si, (NR6 - NR7) / (NR6 + NR7))
si$STI <- with(si, NR6 / NR7)

# Add SI to the database
db <- merge(db, si, by = "ID", all.x = TRUE)
remove(si)

# LAND SURFACE TEMPERATURE ------------------------------------------------

# MOD11A1
lst <- readRDS("data/db02_MOD11A1.rds")
db <- merge(db, lst, by = "ID", all.x = TRUE)
remove(lst)

# MOD11A2
lst2 <- readRDS("data/db02_MOD11A2.rds")
db <- merge(db, lst2, by = "ID", all.x = TRUE)
remove(lst2)

# MODIS LAND COVER --------------------------------------------------------

# Working directory
wdir <- "D:/Data/LFMC/MCD12Q1"

# Iteration counter
tyrs <- sort(unique(db$YEAR))

# Get IGBP LandCover classes
tmp <- lapply(tyrs, function(i) {
    # Subset db by year
    db_i <- db[db$YEAR == i, ]

    # Load raster file
    yr <- ifelse(i == 2000, 2001, i)
    r_file <- paste0("MCD12Q1.006_LC_Type1_doy", yr, "001_aid0001.tif")
    r <- raster(file.path(wdir, r_file))
    
    # Convert db_i to spatial object
    coordinates(db_i) <- ~ LON + LAT
    proj4string(db_i) <- CRS(SRS_string = "EPSG:4326")
    db_i <- spTransform(db_i, CRSobj = crs(r))
    
    # Extract landcover classes
    db_i$IGBP_LC <- extract(r, db_i)
    
    # Output
    db_i
})
lc_df <- as.data.frame(do.call(rbind, tmp))

# Add classes to the database
db <- merge(db, lc_df[, c("ID", "IGBP_LC")], by = "ID", all.x = TRUE)
remove(lc_df, tmp, tyrs)

# FUEL TYPE CLASSIFICATION ------------------------------------------------

# Load LUT
lut <- read.csv("data/landcover_fueltype_LUT.csv")

# Assing classes
db$IGBP_NAME <- lut[match(db$IGBP_LC, lut$IGBP_ID), "IGBP_NAME"]
db$FUEL <- lut[match(db$IGBP_LC, lut$IGBP_ID), "LC"]

# NDVIcv ------------------------------------------------------------------

# Load data
ndvi <- readRDS("./data/db02_ndvi.rds")

# Add NDVIcv to the database
db <- merge(db, ndvi, by = "ID", all.x = TRUE)
remove(ndvi)

# REMOVE UNREALISTIC OBSERVATIONS -----------------------------------------

# Explore unrealistic LFMC measurements
unreal <- db$LFMC_AVG > 20 & db$LFMC_AVG < 250
1 - sum(unreal) / nrow(db)
db[!unreal, c("IGBP_NAME", "SPECIES_COLLECTED")]
table(db[!unreal, "IGBP_NAME"])

# Remove unrealistic observations
db <- db[db$LFMC_AVG >= 20 & db$LFMC_AVG <= 250, ]

# AVERAGE SAMPLES LOCATED WHITHIN THE SAME PIXEL --------------------------

# Load raster template
r_tmpl <- raster("./data/modis_tmp.tif")

# Convert db to spatial object
db_sp <- db
coordinates(db_sp) <- ~ LON + LAT
proj4string(db_sp) <- CRS(SRS_string = "EPSG:4326")
db_sp_sinu <- spTransform(db_sp, CRSobj = crs(r_tmpl))

# Get pixel ID
db$cellID <- cellFromXY(r_tmpl, db_sp_sinu)

# Detect equal date-cellID samples
d1 <- which(duplicated(db[, c("DATE", "cellID")]))
d2 <- which(duplicated(db[, c("DATE", "cellID")], fromLast = TRUE))
dup <- sort(unique(c(d1, d2)))
ss_db <- db[dup, ]

# Assign removing-sample flags
db$PROX_FIL <- 0
ts <- unique(ss_db[, c("DATE", "cellID")])
for (i in seq_len(nrow(ts))) {
    # Current date-cellID samples
    ss <- ss_db[ss_db$DATE == ts[i, "DATE"] &
                    ss_db$cellID == ts[i, "cellID"], ]
    
    # Average LFMC measurements
    ss[, "LFMC_AVG"] <- mean(ss[, "LFMC_AVG"], na.rm = TRUE)
    
    # Random sample selection
    id <- sample(seq_len(nrow(ss)), size = 1)
    db[db$ID %in% ss[-id, "ID"], "PROX_FIL"] <- 1
}

# Remove assigned points which falls at the same pixel
db <- db[db$PROX_FIL == 0, ]
db$cellID <- NULL
db$PROX_FIL <- NULL

# GET LFMC VALUES FROM QUAN2021 -------------------------------------------

# Temporary directory
rtmdir <- "./data/rtm"
if (!dir.exists(rtmdir)) dir.create(rtmdir)

# Iteration to extract values
ncores <- detectCores() - 1
for (yr in 2000:2018) {
    
    # LFMC source directory
    qdir <- file.path("D:/Data/LFMC_Quan2021", yr)
    
    # List of available layers
    r_list <- list.files(qdir)
    r_ydoy <- str_sub(r_list, 1, 7)
    r_list <- r_list[r_ydoy %in% db$YDOY]
    
    # Iterate by layer
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    tmp <- foreach(lyr = r_list) %dopar% {
        library(raster)
        
        # Load raster data
        r <- raster(file.path(qdir, lyr))
        
        # Subset data by date
        lyr_ydoy <- stringr::str_sub(lyr, 1, 7)
        db_i <- db[db$YDOY == lyr_ydoy, ]
        
        # Convert db_i to spatial object
        coordinates(db_i) <- ~ LON + LAT
        proj4string(db_i) <- CRS(SRS_string = "EPSG:4326")
        db_i <- spTransform(db_i, CRSobj = crs(r))
        
        # Extract values
        db_i$RTM <- extract(r, db_i, method = "simple")
        as.data.frame(db_i)
    }
    stopCluster(cl)
    
    # Save output results for the current year
    yr_res <- do.call(rbind, tmp)
    saveRDS(yr_res, file.path(rtmdir, paste0("quan", yr, ".rds")))
}

# Add partial outputs to the database
tmp <- lapply(paste0("quan", 2000:2018, ".rds"), function(x) {
    readRDS(file.path(rtmdir, x))
})
lfmc_rtm <- do.call(rbind, tmp)
db <- merge(db, lfmc_rtm[, c("ID", "RTM")], by = "ID", all.x = TRUE)

# Remove temporary directory
unlink(rtmdir, recursive = TRUE)

# MANUALLY CORRECT ANOMALOUS FUEL TYPE CLASSIFICATIONS --------------------

# Explore land covers
table(db$LAND_COVER, db$FUEL)
table(db$FUEL)
features <- c("LAND_COVER", "IGBP_NAME", "FUEL", "SPECIES_COLLECTED", "REGION")

# Correct Croplands
tst <- db[db$FUEL == "Croplands", features]
unique(tst)
rule1 <- db$FUEL == "Croplands" & db$LAND_COVER == "Shrubland"
rule2<- db$FUEL == "Croplands" & db$SPECIES_COLLECTED == "Unknown grass"
db[rule1, "FUEL"] <- "3_Shrublands"
db[rule2, "FUEL"] <- "4_Grasslands"

# Correct Permanent wetlands
tst <- db[db$FUEL == "Permanent Wetlands", features]
unique(tst)
db[db$FUEL == "Permanent Wetlands", "FUEL"] <- "3_Shrublands"

# Correct Urban covers
tst <- db[db$FUEL == "Urban and Built-up Lands", features]
unique(tst)
db[db$FUEL == "Urban and Built-up Lands", "FUEL"] <- "3_Shrublands"

# Correct by source data
#db[db$SOURCE == "cat", "FUEL"] <- "3_Shrublands"

# FINAL SETTINGS ----------------------------------------------------------

# Format target variable
names(db)[names(db) == "LFMC_AVG"] <- "LFMC"

# SAVE COMPLETE DATABASE --------------------------------------------------

saveRDS(db, "./data/db03_lfmc.rds")

# END ---
