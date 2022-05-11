# *************************************************************************
# LFMC-RF - A remote sensing approach to predict LFMC at large scale using
# Random Forests
# Extract values from MODIS Land Surface Temperature database.
#  - MOD11A2: Terra 8-day average
# Author: Angel Cunill Camprubi <acunill86@gmail.com>
# Last update: 03 Mar 2022
# *************************************************************************

# SETTINGS ----------------------------------------------------------------

# Packages
library(raster)
library(stringr)
library(lubridate)
library(sf)
library(RGISTools)
library(httr)
library(gdalUtils)
library(doParallel)
library(foreach)

# Product
MODIS_prod <- "MOD11A2"

# Working directories
wdir <- file.path("./data", MODIS_prod)
suppressWarnings(dir.create(wdir))
suppressWarnings(dir.create(file.path(wdir, "hdf")))
suppressWarnings(dir.create(file.path(wdir, "tif")))
suppressWarnings(dir.create(file.path(wdir, "tmp_point")))
suppressWarnings(dir.create(file.path(wdir, "out_point")))
outDir <- "data"

# NASA`s authentication
nasa <- file.path("C:/Users/ACC/nasa.netrc")

# Create log files for the process
chk_file <- file.path(wdir, "chk.txt")
if (!file.exists(chk_file)) {
    file.create(chk_file, showWarnings = FALSE)
}
date_file <- file.path(wdir, "MOD11A2_dates.txt")
if (!file.exists(date_file)) {
    file.create(date_file, showWarnings = FALSE)
}
cat(
    sprintf("#START SESSION: %s", as.character(Sys.time())),
    file = chk_file,
    sep = "\n",
    append = TRUE
)

# Turn off factors
options(stringsAsFactors = FALSE)

# MODIS native projection
sinu <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

# LOAD SAMPLE DATA --------------------------------------------------------

# Load database
db <- readRDS("./data/db01_field_lfmc.rds")
db <- db[db$COUNTRY %in% c("Spain", "France", "Italy", "Tunisia"), ]
db <- db[db$COUNTRY %in% c("Spain"), ]

# IDs not done
# notdone <- readRDS(file.path(wdir, "ids_notDone.rds"))
# db <- db[db$ID %in% notdone, ]
# remove(notdone)

d <- ymd("2001-06-22", "2001-06-24", "2002-06-08", "2002-06-09", "2002-06-10")
db <- db[db$DATE %in% d, ]

# Convert to spatial object
sp <- st_as_sf(
    x = db,
    crs = CRS(SRS_string = "EPSG:4326"),
    agr = "constant",
    coords = c("LON", "LAT"),
    remove = FALSE
)
remove(db)

# EXTRACT VALUES BY COUNTRY AND DATE-SITE ---------------------------------

#* Start iteration by country ----
countries <- unique(sp$COUNTRY)
for (country in countries) {
    
    # Subset country samples
    sp_country <- sp[sp$COUNTRY == country, ]
    cat(
        sprintf("**%s**", country),
        file = chk_file,
        sep = "\n",
        append = TRUE
    )
    
    #* Start iteration by date ----
    tseq <- unique(sp_country$DATE)
    for (i in seq_along(tseq)) {
        
        # Select current day sampling locations
        cday <- tseq[i]
        sp_country_cday <- sp_country[sp_country$DATE == cday, ]
        cat(sprintf(
            "\f\nCountry: %s (%i of %i)\nCurrent date: %s (%i of %i)",
            country,
            which(countries == country),
            length(countries),
            as.character(cday),
            i,
            length(tseq)
        ))
        
        # Remove existing files
        cat("\n(1) Remove existing files")
        file.remove(list.files(file.path(wdir, "hdf"), full.names = TRUE))
        file.remove(list.files(file.path(wdir, "tif"), full.names = TRUE))
        file.remove(list.files(file.path(wdir, "tmp_point"), full.names = TRUE))
        
        #** Searching MODIS tiles ----
        cat("\n(2) Searching tiles")
        search_files <- NULL
        try(search_files <- modSearch(product = MODIS_prod,
                                        collection = 6,
                                        dates = cday,
                                        extent = sp_country_cday)[["hdf"]]
        )
        
        # Checking if search connection fails
        if (is.null(search_files)) {
            cat(
                sprintf("- %s Search connection fails.", as.character(cday)),
                file = chk_file,
                sep = "\n",
                append = TRUE
            )
            Sys.sleep(30)
            next #skip current date iteration if connection fails
        }
        
        # Checking if list of files is empty
        empty_chk <- length(unlist(search_files))
        if (empty_chk == 0) {
            cat(
                sprintf("- %s No files found.", as.character(cday)),
                file = chk_file,
                sep = "\n",
                append = TRUE
            )
            next #skip current date iteration if no existing files
        }
        
        # Get MODIS date linked to the target date
        search_files_date <- do.call(
            c,
            lapply(str_split(search_files, "/"), function(x) {
                as.Date(x[8], "%Y.%m.%d")
            })
        )
        if (length(unique(search_files_date)) > 1) {
            stop("More than 1 file date.")
        }
        write(as.character(c(cday, search_files_date[1])),
              file = date_file,
              ncolumns = 2,
              append = TRUE,
              sep = ";")
        
        # Checking if target date exists
        if (length(search_files) == 0) {
            cat(
                sprintf("- %s Target date don't exist.", as.character(cday)),
                file = chk_file,
                sep = "\n",
                append = TRUE
            )
            next #skip current date iteration if no existing files
        }
        
        #** Download MODIS tiles from DataPool ----
        fail <- FALSE
        cat("\n(3) Start downloading files")
        for (fd in search_files) {
            # Keep original filename
            fd_filename <-  tail(str_split(fd, "/")[[1]], n = 1)
            
            # Write file to disk
            response <- NULL
            try(response <-
                    GET(
                        fd,
                        write_disk(file.path(wdir, "hdf", fd_filename), overwrite = TRUE),
                        config(netrc = TRUE, netrc_file = nasa),
                        #authenticate(user = nasa[1], password = nasa[2]),
                        set_cookies("LC" = "cookies")
                    )
            )
            
            # Check if the connection to the server fails
            if (is.null(response)) {
                fail <- TRUE
                Sys.sleep(120)
                break #exit downloading process if connection fails
            }
            
            # Check if file downloaded correctly
            if (response$status_code != 200) {
                fail <- TRUE
                break #exit downloading process if one file fails
            }
        } #end downloading for loop --
        
        # Check for complete files downloading
        if (fail) {
            cat(
                sprintf("- %s Download fails.", as.character(cday)),
                file = chk_file,
                sep = "\n",
                append = TRUE
            )
            next #skip current date iteration if one file fails downloading
        }
        
        # Get list of downloaded HDF files
        local_hdf_files <- list.files(file.path(wdir, "hdf"), pattern = ".hdf$")
        
        #** Export sub-datasets from HDF files to GeoTiff ----
        cat("\n(4) Export datasets from HDF to GeoTiff")
        ncores <- detectCores() - 1
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        foreach(lhdf = local_hdf_files) %dopar% {
            library(gdalUtils)
            library(stringr)
            
            # Get list of datasets
            lhdf_sds <- get_subdatasets(file.path(wdir, "hdf", lhdf))
            
            # Keep tile code id
            lhdf_tile <-
                paste0("_", format(cday, "%Y%j"), "_", str_split(lhdf, "\\.")[[1]][3])
            
            # Extract datasets to designed folder
            for (sds in lhdf_sds[1]) {
                lhdf_filename <-
                    paste0(tail(str_split(sds, ":")[[1]], 1), lhdf_tile, ".tif")
                gdal_translate(sds, dst_dataset = file.path(wdir, "tif", lhdf_filename))
            }
        }
        stopCluster(cl) #end export to GeoTiff for loop --
        
        #** Transform CRS of the locations to match MODIS native CRS ----
        sp_country_cday_points <- st_transform(sp_country_cday, sinu)
        remove(sp_country_cday)
        
        #** Extract values ----
        cat("\n(5) Extracting values")
        
        # Get list of LST tiles
        lst_files <- list.files(file.path(wdir, "tif"),
                                pattern = "LST_Day_1km",
                                full.names = TRUE)
        
        # Extract values for each tile in parallel
        ncores <- detectCores() - 1
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        foreach(tile = seq_along(lst_files)) %dopar% {
            library(raster)
            library(sf)
            library(stringr)
            
            # Keep current tile filename
            tile_filename <- str_replace(
                tail(str_split(lst_files[tile], "/")[[1]], 1),
                ".tif",
                ".csv"
            )
            
            # Load raster tiles
            lst0 <- raster(lst_files[tile])
            
            # Resample
            lst <- disaggregate(lst0, fact = 2, method = "bilinear")
            
            # Extract values by points
            pout_point <- extract(lst,
                                  sp_country_cday_points,
                                  method = "simple",
                                  sp = TRUE)
            pout_point <- pout_point@data[, c(1, ncol(pout_point))]
            
            # Save temporary results
            colnames(pout_point) <- c("ID", "LST8day")
            write.csv(
                pout_point,
                file.path(wdir, "tmp_point", tile_filename),
                row.names = FALSE
            )
        }
        stopCluster(cl) #end extraction of the values --
        
        # Load and merge point based temporary files
        tmp_point_files <- list.files(file.path(wdir, "tmp_point"),
                                      full.names = TRUE,
                                      pattern = ".csv$")
        tmp_point <- lapply(tmp_point_files, read.csv)
        tmp_point_values <- do.call(cbind, lapply(tmp_point, function(x) x[, 2]))
        tmp_point_values <- apply(tmp_point_values, 1, mean, na.rm = TRUE)
        out_point <- data.frame(
            "ID" = tmp_point[[1]][, "ID"],
            "value" = ifelse(is.nan(tmp_point_values), NA, tmp_point_values)
        )
        names(out_point) <- c("ID", "LST8day")
        
        # Save output results
        out_point_filename <- paste0("point_", format(cday, "%Y%j"), ".csv")
        write.csv(
            out_point,
            file.path(wdir, "out_point", out_point_filename),
            row.names = FALSE
        )
        
        cat("\nDONE")
    } #end current date iteration --
} #end country iteration --

# Log finish time
cat(
    sprintf("FINISH SESSION %s", as.character(Sys.time())),
    file = chk_file,
    sep = "\n",
    append = TRUE
)

# JOIN PARTIAL OUTPUTS ----------------------------------------------------

# Load database
db <- readRDS("./data/db01_field_lfmc.rds")
db <- db[db$COUNTRY %in% c("Spain", "France", "Italy", "Tunisia"), ]

# Get list of files
files <- list.files(file.path(wdir, "out_point"),
                    pattern = "csv$",
                    full.names = TRUE)
    
# Join partial outputs
tmp <- lapply(files, function(x) read.csv(x))    
lst_df <- do.call(rbind, tmp)

# Save output results
saveRDS(lst_df, file = file.path(outDir, paste0("db02_", MODIS_prod, ".rds")))

# CHECKINGS ---------------------------------------------------------------

# List not successful extractions
lst_point <- readRDS(file.path(outDir, paste0("db02_", MODIS_prod, ".rds")))
ids_point <- db[!db$ID %in% lst_point$ID, "ID"]
unique(db[!db$ID %in% lst_point$ID, "DATE"])
saveRDS(ids_point, file.path(wdir, "ids_notDone.rds"))

# Checking dates
d <- read.table(file.path(wdir, paste0(MODIS_prod, "_dates.txt")),
                sep = ";",
                header = FALSE)
d$V1 <- ymd(d$V1)
d$V2 <- ymd(d$V2)
d$dif <- as.integer(d$V1 - d$V2)
summary(d)

# REMOVE TEMPORARY FILES --------------------------------------------------

unlink(file.path(wdir, "hdf"))
unlink(file.path(wdir, "tif"))
#unlink(file.path(wdir, "tmp_point"))
#unlink(wdir, recursive = TRUE)

# END ---
