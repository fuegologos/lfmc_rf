#**************************************************************************
#* LFMC - A remote sensing approach to predict LFMC at large scale
#* Preparation of field measurements
#* Author: Angel Cunill Camprubi <acunill86@gmail.com>
#* Last update: 09 Jun 2021
#**************************************************************************

# SETTINGS ----------------------------------------------------------------

# Packages
library(lubridate)
library(raster)
library(stringr)

# Directories
wdir <- "./data/raw"

# Turn off factors
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
    "LFMC",
    "SPECIES_COLLECTED"
)

# GLOBE-LFMC DATABASE -----------------------------------------------------

# Load data
df_globe <- read.csv(
    file.path(wdir, "Globe-LFMC-v2.csv"),
    sep = ";"
)

# Subset fields of interest
df_globe <- df_globe[, c(1, 3:9, 11, 13, 19)]
names(df_globe) <- fields

# Date field
df_globe$DATE <-
    as.Date(as.character(df_globe$DATE), format = "%Y%m%d")

# CATALONIAN DATABASE -----------------------------------------------------

# Load data
df_cat <- read.csv(
    file.path(wdir, "lfmc_catalonia_v1.csv"),
    sep = ";"
)

# Remove outliers
df_cat <- df_cat[is.na(df_cat$ManualOutlier), ]
df_cat <- df_cat[is.na(df_cat$AdditiveOutlier), ]

# Add sites metadata
sites_md <- read.csv(
    file.path(wdir, "lfmc_catalonia_v1_sites.csv"),
    sep = ";"
)
df_cat$SamplingSiteCode <-
    do.call(rbind, str_split(df_cat$SiteSpCode, pattern = "s"))[, 2]
df_cat <- cbind(df_cat, data.frame("REGION" = NA, "LON" = NA, "LAT" = NA))
for (i in unique(df_cat$SamplingSiteCode)) {
    ids <- which(df_cat$SamplingSiteCode == i)
    id_md <- which(sites_md$SamplingSiteCode == i)
    df_cat[ids, "REGION"] <- sites_md[id_md, "LocalityName"]
    df_cat[ids, "LON"] <- sites_md[id_md, "Longitude"]
    df_cat[ids, "LAT"] <- sites_md[id_md, "Latitude"]
}

# Unify sites with duplicated coordinates
which(duplicated(sites_md[, c("Longitude", "Latitude")], fromLast = FALSE))
ids <- which(df_cat$SamplingSiteCode == sites_md[10, "SamplingSiteCode"])
df_cat[ids, "SamplingSiteCode"] <- sites_md[2, "SamplingSiteCode"]
df_cat[ids, "REGION"] <- sites_md[2, "LocalityName"]

# Add species metadata
species_md <- read.csv(
    file.path(wdir, "lfmc_catalonia_v1_species.csv"),
    sep = ";"
)
df_cat$SpeciesCode <-
    do.call(rbind, str_split(df_cat$SiteSpCode, pattern = "sp"))[, 2]
df_cat$SPECIES_COLLECTED <- NA
for (i in unique(df_cat$SpeciesCode)) {
    df_cat[df_cat$SpeciesCode == i, "SPECIES_COLLECTED"] <-
        species_md[species_md$SpeciesCode == i, "SpeciesName"]
}

# Adapt structure to Globe-LFMC db
df_cat$ID <- paste0("C", df_cat$SampleCode)
df_cat$SITE <- paste0("Cat", df_cat$SamplingSiteCode)
df_cat$COUNTRY <- "Spain"
df_cat$TIME <- "12:00"
df_cat$PROTOCOL <- 99
df_cat$DATE <- as.Date(df_cat$Date, "%d/%m/%Y")
df_cat <- df_cat[, fields]

# MERGE TWO DATABASES -----------------------------------------------------

# Merge
df <- rbind(df_globe, df_cat)

# Filter
df <- df[complete.cases(df$LFMC), ]           # remove NAs
df <- df[df$DATE >= "2000-02-24", ]           # min date of MCD43A4 product

# ONE VALUE FOR EACH S44ITE-DATE --------------------------------------------

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
out[out$COUNTRY == "Republic of South Africa", "COUNTRY"] <- "RSA"
out$LFMC <- NULL
out$SPECIES_COLLECTED <- NULL
digits <- str_length(nrow(out))
out$ID <- sprintf(paste0("C%0", digits, "i"), seq_len(nrow(out)))

# Save results
saveRDS(out, "./data/tmp/db01_field_lfmc.rds")

# END ---
