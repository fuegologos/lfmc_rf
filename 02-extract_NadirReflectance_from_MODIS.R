#**************************************************************************
#* LFMC - A remote sensing approach to predict LFMC at large scale
#* Extract values from MODIS MCD43A4 database (BRDF-Reflectances).
#* Author: Angel Cunill Camprubi <acunill86@gmail.com>
#* Last update: 09 Jul 2021
#**************************************************************************

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

# Working directories
wdir <- file.path("./data/tmp/MCD43A4")
suppressWarnings(dir.create(wdir))
suppressWarnings(dir.create(file.path(wdir, "hdf")))
suppressWarnings(dir.create(file.path(wdir, "tif")))
suppressWarnings(dir.create(file.path(wdir, "tmp_point")))
suppressWarnings(dir.create(file.path(wdir, "tmp_area")))
suppressWarnings(dir.create(file.path(wdir, "out_point")))
suppressWarnings(dir.create(file.path(wdir, "out_area")))
outDir <- "data/tmp"

# NASA`s authentication
nasa <- file.path("C:/Users/ACC/nasa.netrc")

# Create log file for the process
chk_file <- file.path(wdir, "chk.txt")
if (!file.exists(chk_file)) {
    file.create(chk_file, showWarnings = FALSE)
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
db <- readRDS("./data/tmp/db01_field_lfmc.rds")

# IDs not done
# notdone <- readRDS(file.path(outDir, "ids_notDone.rds"))
# db <- db[db$ID %in% notdone, ]
# remove(notdone)

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
        
        #** Searching MODIS tiles ----
        cat("\n(2) Searching tiles")
        search_results <- NULL
        try(search_results <- modSearch(product = "MCD43A4",
                                        collection = 6,
                                        dates = cday,
                                        extent = sp_country_cday)[["hdf"]]
        )
        
        # Checking if search connection fails
        if (is.null(search_results)) {
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
        empty_chk <- length(unlist(search_results))
        if (empty_chk == 0) {
            cat(
                sprintf("- %s No files found.", as.character(cday)),
                file = chk_file,
                sep = "\n",
                append = TRUE
            )
            next #skip current date iteration if no existing files
        }
        
        # Get list of files for the target date
        search_files_date <- do.call(
            c,
            lapply(str_split(search_results, "/"), function(x) {
                as.Date(x[8], "%Y.%m.%d")
            })
        )
        search_files_df <- data.frame(
            "hdf" = search_results,
            "filedate" = search_files_date
        )
        search_files <- search_files_df[search_files_df$filedate == cday, "hdf"]
        
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
            for (sds in lhdf_sds) {
                lhdf_filename <-
                    paste0(tail(str_split(sds, ":")[[1]], 1), lhdf_tile, ".tif")
                gdal_translate(sds, dst_dataset = file.path(wdir, "tif", lhdf_filename))
            }
        }
        stopCluster(cl) #end export to GeoTiff for loop --
        
        #** Transform CRS of the locations to match MODIS native CRS ----
        sp_country_cday_points <- st_transform(sp_country_cday, sinu)
        remove(sp_country_cday)
        
        #** Make square sample area surrounding sample points ----
        sp_country_cday_area <- st_buffer(
            sp_country_cday_points,
            #dist = ceiling(res(nr_masked)[1] / 2),
            dist = 232,
            nQuadSegs = 1,
            endCapStyle = "SQUARE",
            joinStyle = "MITRE"
        )
        
        #** Extract values for each band ----
        cat("\n(5) Extracting values")
        for (band in 1:7) {
            cat("\n- Band", band)
            
            # Remove existing temporary files
            file.remove(list.files(file.path(wdir, "tmp_point"), full.names = TRUE))
            file.remove(list.files(file.path(wdir, "tmp_area"), full.names = TRUE))
            
            # Output filename
            band_filename <-
                paste0(country, "_", format(cday, "%Y%j"), "_B", band)
            
            # Get list of nadir reflectance tiles
            nr_files_pattern <- paste0("Nadir_Reflectance_Band", band, ".+tif$")
            nr_files <- list.files(file.path(wdir, "tif"),
                                   pattern = nr_files_pattern,
                                   full.names = TRUE)
            
            # Get list of quality layer tiles
            qa_files_pattern <-
                paste0("BRDF_Albedo_Band_Mandatory_Quality_Band", band, ".+tif$")
            qa_files <- list.files(file.path(wdir, "tif"),
                                   pattern = qa_files_pattern,
                                   full.names = TRUE)
            
            # Extract values for each tile in parallel
            ncores <- detectCores() - 1
            cl <- makeCluster(ncores)
            registerDoParallel(cl)
            foreach(tile = seq_along(nr_files)) %dopar% {
                library(raster)
                library(sf)
                library(stringr)
                
                # Keep current tile filename
                tile_filename <- str_replace(
                    tail(str_split(nr_files[tile], "/")[[1]], 1),
                    ".tif",
                    ".csv"
                )
                
                # Load raster tiles
                nr <- raster(nr_files[tile])
                qa <- raster(qa_files[tile])
                
                # Mask by good quality pixels (0 = good quality, 1 = other QA)
                qa[qa == 1] <- NA
                nr_masked <- raster::mask(nr, mask = qa)
                remove(qa, nr)
                
                # Extract values by points
                pout_point <- extract(nr_masked,
                                      sp_country_cday_points,
                                      method = "simple",
                                      sp = TRUE)
                pout_point <- pout_point@data[, c(1, ncol(pout_point))]
                
                # Save temporary results
                colnames(pout_point) <- c("ID", paste0("B", band, "_", tile))
                write.csv(
                    pout_point,
                    file.path(wdir, "tmp_point", tile_filename),
                    row.names = FALSE
                )
                
                # Extract values by area
                pout_area <- extract(
                    nr_masked,
                    sp_country_cday_area,
                    fun = mean,
                    na.rm = TRUE,
                    weights = TRUE,
                    normalizeWeights = TRUE,
                    sp = TRUE
                )
                if (length(unlist(pout_area)) > 0) {
                    # Get data
                    pout_area <- pout_area@data[, c(1, ncol(pout_area))]
                    
                    # Check if all cells are missing values
                    missing_cells <- extract(
                        nr_masked,
                        sp_country_cday_area,
                        na.rm = FALSE,
                        fun = function(x, ...) {
                            as.integer(all(is.na(x)))
                        }
                    )[, 1]
                    missing_cells[is.na(missing_cells)] <- 1
                    pout_area[missing_cells == 1, 2] <- NA
                    
                    # Save temporary results
                    colnames(pout_area) <- c("ID", paste0("B", band, "_", tile))
                    write.csv(
                        pout_area,
                        file.path(wdir, "tmp_area", tile_filename),
                        row.names = FALSE
                    )
                }
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
            names(out_point) <- c("ID", paste0("B", band))
            
            # Save output results
            out_point_filename <- paste0("point_", band_filename, ".csv")
            write.csv(
                out_point,
                file.path(wdir, "out_point", out_point_filename),
                row.names = FALSE
            )
            
            # Load and merge area based temporary files
            tmp_area_files <- list.files(file.path(wdir, "tmp_area"),
                                         full.names = TRUE,
                                         pattern = ".csv$")
            tmp_area <- lapply(tmp_area_files, read.csv)
            tmp_area_values <- do.call(cbind, lapply(tmp_area, function(x) x[, 2]))
            tmp_area_values <- apply(tmp_area_values, 1, mean, na.rm = TRUE)
            out_area <- data.frame(
                "ID" = tmp_area[[1]][, "ID"],
                "value" = ifelse(is.nan(tmp_area_values), NA, tmp_area_values)
            )
            names(out_area) <- c("ID", paste0("B", band))
            
            # Save output results
            out_area_filename <- paste0("area_", band_filename, ".csv")
            write.csv(
                out_area,
                file.path(wdir, "out_area", out_area_filename),
                row.names = FALSE
            )
            
            cat(" DONE")
        } #end band iteration --
    } #end current date iteration --
} #end country iteration --

# Log finish time
cat(
    sprintf("FINISH SESSION %s", as.character(Sys.time())),
    file = chk_file,
    sep = "\n",
    append = TRUE
)

# REMOVE TEMPORARY FILES --------------------------------------------------

file.remove(list.files(file.path(wdir, "hdf"), full.names = TRUE))
file.remove(list.files(file.path(wdir, "tif"), full.names = TRUE))
file.remove(list.files(file.path(wdir, "tmp_point"), full.names = TRUE))
file.remove(list.files(file.path(wdir, "tmp_area"), full.names = TRUE))

# JOIN PARTIAL OUTPUTS ----------------------------------------------------

# Load database
db <- readRDS("./data/tmp/db01_field_lfmc.rds")

# Make joins by extraction method
method <- c("point", "area")
for (m in method) {
    # Get list of files
    files_path <- file.path(wdir, paste0("out_", m))
    files <- list.files(files_path, pattern = ".csv$")
    
    # Join partial outputs
    for (band in 1:7) {
        band_files <- files[str_which(files, paste0("B", band))]
        tmp <- lapply(band_files, function(x) {
            read.csv(file.path(files_path, x))
        })
        band_df <- do.call(rbind, tmp)
        if (band == 1) {
            nr_df <- band_df
        } else {
            nr_df <- merge(
                nr_df,
                band_df,
                by = "ID",
                all = TRUE,
                sort = TRUE
            )
        }    
    }
    
    # Save output results
    names(nr_df)[-1] <- str_replace(names(nr_df)[-1], "B", "NR")
    saveRDS(nr_df, file = file.path(outDir, paste0("db02_", m, "_nr.rds")))
}

# List not successful extractions
nr_point <- readRDS(file.path(outDir, "db02_point_nr.rds"))
nr_area <- readRDS(file.path(outDir, "db02_area_nr.rds"))
ids_point <- db[!db$ID %in% nr_point$ID, "ID"]
ids_area <- db[!db$ID %in% nr_area$ID, "ID"]
identical(ids_point, ids_area)
length(ids_area)
unique(db[!db$ID %in% nr_area$ID, "DATE"])
#saveRDS(ids_area, file.path(outDir, "ids_notDone.rds"))

# END ---
