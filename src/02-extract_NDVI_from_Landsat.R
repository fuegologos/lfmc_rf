# *************************************************************************
# LFMC-RF - A remote sensing approach to predict LFMC at large scale using
# Random Forests
# Get values of NDVI from Landsat 5/7/8 from the extractions done 
# using Google Earth Engine and the scripts in src>gee folder.
# - Preparation of sample locations before running GEE
# - Get values from exported CSV files
# Author: Angel Cunill Camprubi <acunill86@gmail.com>
# Last update: 03 Mar 2022
# *************************************************************************

# SETTINGS ----------------------------------------------------------------

# Packages
library(raster)
library(sf)
library(stringr)
library(lubridate)

# Others
options(stringsAsFactors = FALSE)

# LOAD DATABASE -----------------------------------------------------------

# Load database
db <- readRDS("./data/db01_field_lfmc.rds")
summary(db)

# Subset Mediterranean countries
unique(db$COUNTRY)
db <- db[db$COUNTRY %in% c("Spain", "France", "Italy", "Tunisia"), ]

# Assign month and year fields
db$MONTH <- month(db$DATE)
db$YEAR <- year(db$DATE)

# CREATE INPUT SHAPEFILE FOR GEE (NOT RUN) --------------------------------

# Load raster template
r_tmpl <- raster("./data/modis_tmp.tif")

# Convert db to spatial object
db_sp <- db
coordinates(db_sp) <- ~ LON + LAT
proj4string(db_sp) <- CRS("+init=epsg:4326")
db_sp_sinu <- spTransform(db_sp, CRSobj = crs(r_tmpl))

# Get nearest cell coordinates
cellID <- cellFromXY(r_tmpl, db_sp_sinu)
db$cellLON <- xyFromCell(r_tmpl, cell = cellID)[, 1]
db$cellLAT <- xyFromCell(r_tmpl, cell = cellID)[, 2]

# Convert unique cell centroids to spatial points
db_sp2 <- db[!duplicated(db$SITE), ]
coordinates(db_sp2) <- ~ cellLON + cellLAT
proj4string(db_sp2) <- crs(r_tmpl)
db_sp3 <- spTransform(db_sp2, CRSobj = CRS("+init=epsg:4326"))

# Save cell points as shapefile
shapefile(db_sp3[, c("SITE")],
          filename = "./data/shp/plots_cellcoords.shp",
          overwrite = TRUE)

# GET VALUES FROM LANDSAT 5 -----------------------------------------------

# NDVI files path
wdir <- "D:/GIS/LANDSAT/LFMC_NDVI/landsat5"

# Create variables
db$L5_NDVIm <- NA
db$L5_NDVIsd <- NA
db$L5_totPixels <- NA
db$L5_unmaskPixels <- NA

# Iteration to assign values by year-month
for (yr in 2000:2011) {
    for (mth in 1:12) {
        cat(sprintf("\n%i-%i ", mth, yr))
        
        # Load NDVI database
        ndvi_file <- paste0("L5_month_", mth, "_year_", yr, ".csv")
        cat(ndvi_file)
        if (!file.exists(file.path(wdir, ndvi_file))) {
            cat(" <- Do not exist.")
            next
        }
        ndvi <- read.csv(file.path(wdir, ndvi_file))
        
        # Sample sites of current month-year
        ss <- which(db$YEAR == yr & db$MONTH == mth)
        
        # Match site IDs
        ids <- match(db[ss, "SITE"], ndvi$SITE)
        
        # Assign values
        db[ss, "L5_NDVIm"] <- ndvi[ids, "ndviMean"]
        db[ss, "L5_NDVIsd"] <- ndvi[ids, "ndviSD"]
        db[ss, "L5_totPixels"] <- ndvi[ids, "totPixels"]
        db[ss, "L5_unmaskPixels"] <- ndvi[ids, "unmaskedPixels"]
    }
}

# Validate NDVI values
db$L5_NDVIcv <- db$L5_NDVIsd / abs(db$L5_NDVIm)
db$L5_valid <- db$L5_unmaskPixels / db$L5_totPixels

# GET VALUES FROM LANDSAT 7 -----------------------------------------------

# NDVI files path
wdir <- "D:/GIS/LANDSAT/LFMC_NDVI/landsat7"

# Create variables
db$L7_NDVIm <- NA
db$L7_NDVIsd <- NA
db$L7_totPixels <- NA
db$L7_unmaskPixels <- NA

# Iteration to assign values by year-month
for (yr in 2000:2013) {
    for (mth in 1:12) {
        cat(sprintf("\n%i-%i ", mth, yr))
        
        # Load NDVI database
        ndvi_file <- paste0("L7_month_", mth, "_year_", yr, ".csv")
        cat(ndvi_file)
        if (!file.exists(file.path(wdir, ndvi_file))) {
            cat(" <- Do not exist.")
            next
        }
        ndvi <- read.csv(file.path(wdir, ndvi_file), header = TRUE)
        
        # Sample sites of current month-year
        ss <- which(db$YEAR == yr & db$MONTH == mth)
        
        # Match site IDs
        ids <- match(db[ss, "SITE"], ndvi$SITE)
        
        # Assign values
        db[ss, "L7_NDVIm"] <- ndvi[ids, "ndviMean"]
        db[ss, "L7_NDVIsd"] <- ndvi[ids, "ndviSD"]
        db[ss, "L7_totPixels"] <- ndvi[ids, "totPixels"]
        db[ss, "L7_unmaskPixels"] <- ndvi[ids, "unmaskedPixels"]
    }
}

# Validate NDVI values
db$L7_NDVIcv <- db$L7_NDVIsd / abs(db$L7_NDVIm)
db$L7_valid <- db$L7_unmaskPixels / db$L7_totPixels

# GET VALUES FROM LANDSAT 8 -----------------------------------------------

# NDVI files path
wdir <- "D:/GIS/LANDSAT/LFMC_NDVI/landsat8"

# Create variables
db$L8_NDVIm <- NA
db$L8_NDVIsd <- NA
db$L8_totPixels <- NA
db$L8_unmaskPixels <- NA

# Iteration to assign values by year-month
for (yr in 2013:2019) {
    for (mth in 1:12) {
        cat(sprintf("\n%i-%i ", mth, yr))
        
        # Load NDVI database
        ndvi_file <- paste0("L8_month_", mth, "_year_", yr, ".csv")
        cat(ndvi_file)
        if (!file.exists(file.path(wdir, ndvi_file))) {
            cat(" <- Do not exist.")
            next
        }
        ndvi <- read.csv(file.path(wdir, ndvi_file), header = TRUE)
        
        # Sample sites of current month-year
        ss <- which(db$YEAR == yr & db$MONTH == mth)
        
        # Match site IDs
        ids <- match(db[ss, "SITE"], ndvi$SITE)
        
        # Assign values
        db[ss, "L8_NDVIm"] <- ndvi[ids, "ndviMean"]
        db[ss, "L8_NDVIsd"] <- ndvi[ids, "ndviSD"]
        db[ss, "L8_totPixels"] <- ndvi[ids, "totPixels"]
        db[ss, "L8_unmaskPixels"] <- ndvi[ids, "unmaskedPixels"]
    }
}

# Validate NDVI values
db$L8_NDVIcv <- db$L8_NDVIsd / abs(db$L8_NDVIm)
db$L8_valid <- db$L8_unmaskPixels / db$L8_totPixels

# FILTERING RESULTS -------------------------------------------------------

str(db)
summary(db$DATE)

# Landsat 5
l5_ts <- seq(ymd("2000-03-01"), ymd("2011-10-31"), by = 1)
db[!db$DATE %in% l5_ts, "L5_valid"] <- NA
db[is.na(db$L5_valid), "L5_valid"] <- 0

# Landsat 7
l7_ts <- seq(ymd("2000-03-01"), ymd("2013-04-30"), by = 1)
db[!db$DATE %in% l7_ts, "L7_valid"] <- NA
db[is.na(db$L7_valid), "L7_valid"] <- 0

# Landsat 8
l8_ts <- seq(ymd("2013-05-01"), ymd("2020-01-01"), by = 1)
db[!db$DATE %in% l8_ts, "L8_valid"] <- NA
db[is.na(db$L8_valid), "L8_valid"] <- 0

# Assing NDVI coefficient of variation
db$NDVIcv <- NA
for (i in seq_len(nrow(db))) {
    if (db$L7_valid[i] >= 0.8) {
        db$NDVIcv[i] <- db$L7_NDVIcv[i]
    } 
    if (db$L5_valid[i] >= 0.8) {
        db$NDVIcv[i] <- db$L5_NDVIcv[i]
    }
    if (db$L8_valid[i] >= 0.8) {
        db$NDVIcv[i] <- db$L8_NDVIcv[i]
    }
}

# SAVE OUTPUT RESULTS -----------------------------------------------------

saveRDS(db[, c("ID", "NDVIcv")], "./data/db02_ndvi.rds")

# END ---
