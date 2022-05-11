# *************************************************************************
# LFMC-RF - A remote sensing approach to predict LFMC at large scale using
# Random Forests
# Preparation of field measurements
# Author: Angel Cunill Camprubi <acunill86@gmail.com>
# Last update: 03 Mar 2022
# *************************************************************************

# SETTINGS ----------------------------------------------------------------

# Packages
library(lubridate)
library(stringr)

# Directories
in_dir <- "./data/field_data"

# Others
options(stringsAsFactors = FALSE)

# Fields of interest
fields <- c(
    "ID",
    "SITE",
    "REGION",
    "COUNTRY",
    "LAT",
    "LON",
    "TIME",
    "DATE",
    "PROTOCOL",
    "LAND_COVER",
    "LFMC",
    "SPECIES_COLLECTED",
    "SOURCE"
)

# GLOBE-LFMC DATABASE -----------------------------------------------------

# Load data
df_globe <- read.csv(
    file.path(in_dir, "Globe-LFMC-v2.csv"),
    sep = ";"
)
df_globe$SOURCE <- "globe"

# Subset fields of interest
df_globe <- df_globe[, c(1, 3:9, 11:13, 19, 22)]
names(df_globe) <- fields

# Convert Date field
df_globe$DATE <-
    as.Date(as.character(df_globe$DATE), format = "%Y%m%d")

# CATALONIAN DATABASE -----------------------------------------------------

# Load data
df_cat <- read.csv(
    file.path(in_dir, "lfmc_catalonia_v1.csv"),
    sep = ";"
)

# Remove outliers
df_cat <- df_cat[is.na(df_cat$ManualOutlier), ]
df_cat <- df_cat[is.na(df_cat$AdditiveOutlier), ]

# Add sites metadata
sites_md <- read.csv(
    file.path(in_dir, "lfmc_catalonia_v1_sites.csv"),
    sep = ";"
)
df_cat$SamplingSiteCode <-
    do.call(rbind, str_split(df_cat$SiteSpCode, pattern = "s"))[, 2]
ids <- match(df_cat$SamplingSiteCode, sites_md$SamplingSiteCode)
df_cat <- cbind(df_cat,
                "REGION" = sites_md[ids, c("LocalityName")],
                "LON" = sites_md[ids, c("Longitude")],
                "LAT" = sites_md[ids, c("Latitude")])

# Unify sites with duplicated coordinates
which(duplicated(sites_md[, c("Longitude", "Latitude")], fromLast = FALSE))
ids <- which(df_cat$SamplingSiteCode == sites_md[10, "SamplingSiteCode"])
df_cat[ids, "SamplingSiteCode"] <- sites_md[2, "SamplingSiteCode"]
df_cat[ids, "REGION"] <- sites_md[2, "LocalityName"]

# Add species metadata
species_md <- read.csv(
    file.path(in_dir, "lfmc_catalonia_v1_species.csv"),
    sep = ";"
)
df_cat$SpeciesCode <-
    do.call(rbind, str_split(df_cat$SiteSpCode, pattern = "sp"))[, 2]
ids <- match(df_cat$SpeciesCode, species_md$SpeciesCode)
df_cat <- cbind(df_cat,
                "SPECIES_COLLECTED" = species_md[ids, "SpeciesName"])

# Adapt structure to Globe-LFMC db
df_cat$ID <- paste0("C", df_cat$SampleCode)
df_cat$SITE <- paste0("Cat", df_cat$SamplingSiteCode)
df_cat$COUNTRY <- "Spain"
df_cat$TIME <- "12:00"
df_cat$PROTOCOL <- 99
df_cat$LAND_COVER <- "Shrubland"
df_cat$DATE <- as.Date(df_cat$Date, "%d/%m/%Y")
df_cat$SOURCE <- "cat"
df_cat <- df_cat[, fields]

# MERGE TWO DATABASES -----------------------------------------------------

# Merge
df <- rbind(df_globe, df_cat)

# Filter
df <- df[complete.cases(df$LFMC), ]           # remove NAs
df <- df[df$DATE >= "2000-02-24", ]           # min date of MCD43A4 product

# ONE VALUE FOR EACH SITE-DATE --------------------------------------------

# Identify unique site-date samples
vars <- c("DATE", "LAT", "LON")
sitedate <- unique(df[, vars])

# Aggregate duplicated values
df_aggr <- lapply(seq_len(nrow(sitedate)), function(i) {
    d <- sitedate[i, ]
    ss <- df[df$DATE == d$DATE &
                 df$LON == d$LON &
                 df$LAT == d$LAT,]
    tmp <- ss[1, ]
    tmp$SPECIES_COLLECTED <- paste(ss$SPECIES_COLLECTED, collapse = ";")
    tmp$LFMC_AVG <- mean(ss$LFMC, na.rm = TRUE)
    tmp$LFMC_MED <- median(ss$LFMC, na.rm = TRUE)
    tmp
})
out <- do.call(rbind, df_aggr)

# OUTPUT ------------------------------------------------------------------

# Tidy output
out <- out[order(out$DATE, out$LON, out$LAT), ]
out$LFMC <- NULL
digits <- str_length(nrow(out))
out$ID <- sprintf(paste0("C%0", digits, "i"), seq_len(nrow(out)))

# Save results
saveRDS(out, "./data/db01_field_lfmc.rds")

# END ---
