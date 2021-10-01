#**************************************************************************
#* LFMC - A remote sensing approach to predict LFMC at large scale
#* Create database for LFMC modelling process
#*  - Target variable: LFMC_MED
#*  - Extraction protocol: Simple pixel extraction
#* Author: Angel Cunill Camprubi <acunill86@gmail.com>
#* Last update: 07 Sept 2021
#**************************************************************************

# SETTINGS ----------------------------------------------------------------

# Packages
library(raster)
library(sf)
library(stringr)
library(lubridate)
library(doParallel)
library(foreach)

# Others
options(stringsAsFactors = FALSE)
set.seed(1234)

# LOAD DATABASE -----------------------------------------------------------

# Load database
db <- readRDS("./data/tmp/db01_field_lfmc.rds")

# Subset Mediterranean countries
db <- db[db$LAT > 35 & db$LAT < 48, ]
db <- db[db$LON > -10 & db$LON < 35, ]

# Assign month and year fields
db$MONTH <- month(db$DATE)
db$YEAR <- year(db$DATE)
db$YDOY <- format(db$DATE, "%Y%j")

# CALCULATE SPECTRAL INDICES ----------------------------------------------

# Load data
si <- readRDS("./data/tmp/db02_point_nr.rds")

# Spectral indices
si$NDVI <- with(si, (NR2 - NR1) / (NR2 + NR1))
si$NDII6 <- with(si, (NR2 - NR6) / (NR2 + NR6))
si$NDII7 <- with(si, (NR2 - NR7) / (NR2 + NR7))
si$GVMI <- with(si, ((NR2 + 0.1) - (NR6 + 0.02)) / ((NR2 + 0.1) + (NR6 + 0.02)))
si$NDWI <- with(si, (NR2 - NR5) / (NR2 + NR5))
si$EVI <- with(si, (2.5 * (NR2 - NR1)) / (NR2 + (6 * NR1) - (7.5 * NR3) + 1))
si$SAVI <- with(si, (1 + 0.5) * ((NR2 - NR1) / (NR2 + NR1 + 0.5)))
si$VARI <- with(si, (NR4 - NR1) / (NR4 + NR1 - NR3))
si$VIgreen <- with(si, (NR4 - NR1) / (NR4 + NR1))
si$NDTI <- with(si, (NR6 - NR7) / (NR6 + NR7))
si$MSI <- with(si, NR6 / NR2)
si$Gratio <- with(si, NR4 / NR1)

# Add SI to the database
db <- merge(db, si, by = "ID", all.x = TRUE)
remove(si)

# HYBRID FOREST MASK ------------------------------------------------------

# Load hybrid forest mask
hybfor <- raster("D:/GIS/HybridForestMask/for2000_bg.tif")

# Convert db to spatial object
db_sp <- db
coordinates(db_sp) <- ~ LON + LAT
proj4string(db_sp) <- CRS(SRS_string = "EPSG:4326")
db_sp <- spTransform(db_sp, CRSobj = crs(hybfor))

# Extract values
hybfor_ft <- extract(hybfor, db_sp, method = "simple")
hybfor_ft[hybfor_ft == 128] <- 0
db$HYBFOR <- hybfor_ft
remove(hybfor, hybfor_ft, db_sp)

# GET LEAF AREA INDEX -----------------------------------------------------

# # Load data
# lai <- readRDS("./data/tmp/db02_point_lai.rds")
# 
# # Add LAI to the database
# db <- merge(db, lai, by = "ID", all.x = TRUE)
# remove(lai)

# FUEL TYPE CLASSIFICATION ------------------------------------------------

# Iteration counter
tyrs <- sort(unique(db$YEAR))

# Get IGBP LandCover classes
tmp <- lapply(tyrs, function(i) {
    # Subset db by year
    db_i <- db[db$YEAR == i, ]

    # Load raster file
    yr <- ifelse(i == 2000, 2001, i)
    r <- raster(paste0("D:/GIS/MODIS/MCD12Q1/MCD12Q1.006_LC_Type1_doy", yr, "001_aid0001.tif"))
    
    # Convert db_i to spatial object
    coordinates(db_i) <- ~ LON + LAT
    proj4string(db_i) <- CRS(SRS_string = "EPSG:4326")
    db_i <- spTransform(db_i, CRSobj = crs(r))
    
    # Extract landcover classes
    db_i$IGBP_LC <- extract(r, db_i)
    
    # Output
    db_i
})
lc <- as.data.frame(do.call(rbind, tmp))

# Add classes to the database
db <- merge(db, lc[, c("ID", "IGBP_LC")], by = "ID", all.x = TRUE)
remove(lc, tmp, tyrs)

# GET NDVIcv --------------------------------------------------------------

# Load data
ndvi <- readRDS("./data/tmp/db02_ndvi.rds")

# Add NDVIcv to the database
db <- merge(db, ndvi, by = "ID", all.x = TRUE)
remove(ndvi)

# TIME-SERIES FILTERS -----------------------------------------------------

# Load functions
source("src/functions/hampel.R")
source("src/functions/quan.R")

# Apply outlier identifiers
tmp <- lapply(unique(db$SITE), function(x) {
    # Subset site samples
    site_db <- db[db$SITE == x, ]

    # Create fields
    site_db$H3 <- 999
    site_db$H5 <- 999
    site_db$Q3 <- 999

    # 3-step windows identifiers
    if (nrow(site_db) >= 3) {
        site_db$H3 <- hampel(x = site_db$LFMC_MED, k = 1)
        site_db$Q3 <- quan(x = site_db$LFMC_MED, k = 1)
    }

    # 5-step windows identifier
    if (nrow(site_db) >= 5) {
        site_db$H5 <- hampel(x = site_db$LFMC_MED, k = 2)
    }

    # Output
    site_db[, c("ID", "H3", "H5", "Q3")]
})
tf <- do.call(rbind, tmp)

# Correct NaN values
tf[is.na(tf$H3), "H3"] <- 999
tf[tf$H3 == Inf, "H3"] <- 999

# Add identifiers to the database
db <- merge(db, tf, by = "ID", all.x = TRUE)
remove(tmp, tf)

# SPATIAL PROXIMITY FILTER ------------------------------------------------

# Load raster template
r <- raster("./data/tmp/modis_tmp.tif")

# Convert db to spatial object
db_sp <- db
coordinates(db_sp) <- ~ LON + LAT
proj4string(db_sp) <- CRS(SRS_string = "EPSG:4326")
db_sp_sinu <- spTransform(db_sp, CRSobj = crs(r))

# Get pixel ID
db$cellID <- cellFromXY(r, db_sp_sinu)

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
moddir <- "./data/tmp/modis"
if (!dir.exists(moddir)) dir.create(moddir)

# Iteration to extract values
ncores <- detectCores() - 1
for (yr in 2000:2018) {
    
    # LFMC source directory
    qdir <- file.path("D:/GIS/LFMC_QUAN2021", yr)
    
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
    saveRDS(yr_res, file.path(moddir, paste0("quan", yr, ".rds")))
}

# Add partial outputs to the database
tmp <- lapply(paste0("quan", 2000:2018, ".rds"), function(x) {
    readRDS(file.path(moddir, x))
})
lfmc_rtm <- do.call(rbind, tmp)
db <- merge(db, lfmc_rtm[, c("ID", "RTM")], by = "ID", all.x = TRUE)

# Remove temporary directory
unlink(moddir, recursive = TRUE)

# SAVE COMPLETE DATABASE --------------------------------------------------

# Remove missing values
si <- c('NDVI', 'NDII6', 'NDII7', 'GVMI', 'NDWI', 'EVI', 'SAVI', 'VARI', 'VIgreen', 'NDTI', 'MSI', 'Gratio')
bands <- c('NR1', 'NR2', 'NR3', 'NR4', 'NR5', 'NR6', 'NR7')
vars <- c('HYBFOR', si, bands)
db <- db[complete.cases(db[, vars]), ]

# Save database
saveRDS(db, "./data/tmp/db03_lfmc.rds")

# END ---
