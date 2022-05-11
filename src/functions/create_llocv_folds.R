# *************************************************************************
# LFMC - A remote sensing approach to predict LFMC at large scale
# Description: Stratifies spatial data into k-folds for cross-validation
#   1. Locations are clustered using k-means and their coordinates
#   2. Each cluster data split into k-folds
#   3. List of fold indices is obtained by merging those created for each cluster
# Arguments:
#   locations - column names that refer to spatial coordinates
#   id        - column name that refers to the index of each location
#   data      - input data matrix with a column ID identifying each row
#   kfold     - number of folds
#   nmeans    - number of centers for k-means clustering
#   seed      - seed for random numbers
#   nstart    - number of seeds (seed + sequence) for multiple initial cluster assignments
#   out_df    - output original dataframe
# Return value:
#   data_out  - locations w/ indices and cluster and fold assignment
#   fold_list - list of folds containing the incides of observations
# Author: Angel Cunill Camprubi <acunill86@gmail.com>
# Last update: 07 Jan 2021
# *************************************************************************

create_llocv_folds <- function(locations,
                               id,
                               data,
                               kfold = 5,
                               nmeans = kfold,
                               seed = sample(1:1000, 1),
                               nstart = 10,
                               out_df = FALSE) {

  # Unique sites
  n_samples <- aggregate(data$ID, list(data[, id]), length)
  colnames(n_samples) <- c("SITE", "N")
  loc <- data[!duplicated(data$SITE), c("SITE", "LON", "LAT")]
  df <- merge(n_samples, loc, by = "SITE")
  
  # Clustering spatial data
  seed_seq <- seed + seq(0, nstart - 1, 1)
  quality <- as.numeric(NA, length(seed_seq))
  min_size <- as.numeric(NA, length(seed_seq))
  for (i in seq_along(seed_seq)) {
    set.seed(seed_seq[i])
    km <- kmeans(df[, locations], centers = nmeans)
    quality[i] <- km$tot.withinss
    min_size[i] <- min(km$size)
  }
  min_size[quality != min(quality)] <- NA
  best_seed <- which.max(min_size)
  set.seed(seed_seq[best_seed])
  km <- kmeans(df[, locations], centers = nmeans)
  df$cluster <- km$cluster
  
  # Split clusters into k-folds
  set.seed(seed)
  clusters <- seq_along(unique(df$cluster))
  suppressWarnings(remove(out))
  for (j in clusters) {
    tmp <- df[df$cluster == j, ]
    tmp$fold <- NA
    fld <- caret::createFolds(seq_len(nrow(tmp)), k = kfold)
    for (i in seq_along(fld)) {
      fld_ids <- fld[[i]]
      tmp[fld_ids, "fold"] <- i
    }
    if (exists("out")) {
      out <- rbind(out, tmp)
    } else {
      out <- tmp
    }
  }
  df <- merge(df, out[, c(id, "fold")], by = id)
  
  # Assing fold label to each sample
  data$FOLD <- NA
  for (i in seq_len(nrow(df))) {
    ids <- which(data[, id] == df[i, "SITE"])
    data[ids, "FOLD"] <- df[i, "fold"]
  }
  
  # Output
  if (out_df) {
    data
  } else {
    data$FOLD 
  }
}

# END ---
