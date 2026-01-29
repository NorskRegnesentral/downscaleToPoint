library(data.table)
library(matrixStats)
library(geosphere)
library(parallel)
library(purrr)
library(downscaleToPoint)

# Define all necessary paths
# ------------------------------------------------------------------------------
data_dir = "/nr/samba/user/smvandeskog/projects/downscaleToPoint/data/"
variables = c("precipitation", "temperature")
model_dirs = file.path(data_dir, "models", variables)
out_dir = file.path(data_dir, "bivariate-mean-simulations")
local_fit_dirs = file.path(model_dirs, "local_fits")

if (!dir.exists(out_dir)) dir.create(out_dir)

meta_path = file.path(data_dir, "meta.rds")
global_fit_paths = file.path(model_dirs, "global.rds")

# Define other useful variables
# ------------------------------------------------------------------------------
n_cores = 8 # Number of cores to use for running code in parallel
n_sims = 1e3 # Number of ensembles to simulate during the downscaling

# Random seeds for reproducibility
set.seed(20260128)
base_seed = sample.int(1e8, 1)
seed_jump = sample.int(1e4, 1)

K = 15

overwrite = FALSE

# Load the meta data
# ------------------------------------------------------------------------------
station_meta = readRDS(meta_path)

# Remove stations with few precipitation observations
station_meta = station_meta[n_good_precip_flag > 200]
station_meta = station_meta[n_unique_precip > 40]
station_meta = station_meta[n_tmean > 200]

# ==============================================================================
# Simulate temperature and precipitation from the full models
# ==============================================================================

global_fits = readRDS(global_fit_paths[variables == "precipitation"])
global_fits$temperature = readRDS(global_fit_paths[variables == "temperature"])

start_time = Sys.time()
parallel::mclapply(
  X = seq_len(nrow(station_meta)),
  mc.cores = n_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {

    set.seed(base_seed + i * seed_jump)

    # Print our progress so far
    time_passed = Sys.time() - start_time
    message(
      "Starting on iter nr. ", i, " / ", nrow(station_meta), " with K=", K,
      ". Time passed: ", round(as.numeric(time_passed), 2), " ", attr(time_passed, "units")
    )

    out_path = file.path(out_dir, paste0(station_meta$id[i], ".rds"))
    if (!overwrite && file.exists(out_path)) return(TRUE)

    # Compute distances to all other weather stations
    dists = geosphere::distHaversine(
      p1 = station_meta[i, c(lon, lat)],
      p2 = station_meta[, cbind(lon, lat)]
    )

    # Locate and load the local models from the K nearest weather stations
    # to weather station nr. i
    local_fits = list()
    for (index in order(dists)[-1]) {
      paths = file.path(local_fit_dirs, paste0(station_meta$id[index], ".rds"))
      if (!all(file.exists(paths))) next
      fits = lapply(paths, readRDS)
      fits = merge(fits[[1]], fits[[2]], by = "id")
      fits$dist = dists[index]
      local_fits[[length(local_fits) + 1]] = fits
      if (length(local_fits) == K) break
    }
    local_fits = data.table::rbindlist(local_fits)

    # Load the data for the current weather station, and add necessary covariates
    data = load_station_data(
      meta = station_meta[i, ],
      data_dir = data_dir,
      rm_bad_flags = TRUE,
      rm_na = FALSE
    )
    data[, let(
      yday = yday(date),
      year = year(date),
      month = month(date),
      station_elevation = log(station_elevation + 1),
      era_log_precip = log(era_precip + 1),
      era_precip_bool = era_precip > 0,
      precip_bool = precip > 0,
      day_count = as.integer(date) - as.integer(min(date)) + 1L
    )]
    data[, let(
      next_era_precip_bool = c(tail(era_precip_bool, -1), NA),
      next_era_log_precip = c(tail(era_log_precip, -1), NA),
      next_era_tmean = c(tail(era_tmean, -1), NA)
    )]
    data[, let(
      era_log_precip_change = next_era_log_precip - era_log_precip,
      era_tmean_change = next_era_tmean - era_tmean
    )]

    # Add offset from the global GAMs
    offsets = list()
    for (name in names(global_fits)) {
      if (grepl("wet", name)) {
        offsets[[name]] = fast_mgcv_pred(global_fits[[name]], data[-nrow(data)])
      } else {
        offsets[[name]] = fast_mgcv_pred(global_fits[[name]], data)
      }
    }

    # Preallocate a list that will hold all simulations from all our different models
    sims = list()

    # Simulate temperature data using the local GAMs, but not the local ARMA models
    sims$tmean = simulate_tmean_with_donors(
      n_sims = n_sims,
      data = data,
      local_fits = local_fits,
      offset = offsets$temperature,
      time_dep = TRUE
    )

    # Check if any of the donor stations appear to be outliers and
    # Remove them if this is the case
    bad_tmean_donors = unique(c(
      get_bad_donor_index(sims$tmean, mean, K),
      get_bad_donor_index(sims$tmean, sd, K)
    ))
    if (length(bad_tmean_donors) > 0) {
      sims$tmean = simulate_tmean_with_donors(
        n_sims = n_sims,
        data = data,
        local_fits = local_fits[-bad_tmean_donors],
        offset = offsets$temperature,
        time_dep = TRUE
      )
    }

    sims$precip = simulate_precip_with_donors(
      n_sims = n_sims,
      data = data,
      local_fits = local_fits,
      offsets = offsets,
      occurrence_time_dep = TRUE,
      intensity_time_dep = TRUE
    )

    bad_precip_donors = get_bad_donor_index(sims$precip, mean, K)
    if (length(bad_precip_donors) > 0) {
      sims$precip = simulate_precip_with_donors(
        n_sims = n_sims,
        data = data,
        local_fits = local_fits[-bad_precip_donors, ],
        offsets = offsets,
        intensity_time_dep = TRUE,
        occurrence_time_dep = TRUE
      )
    }

    ensemble_means = sapply(sims, function(x) apply(x, 1, mean))

    na_index = data[, union(which(is.na(precip)), which(is.na(tmean)))]
    data = data[-na_index]
    ensemble_means = ensemble_means[-na_index, ]

    out = rbind(
      data[, .(date, precip, tmean, tag = "obs")],
      data.table(
        date = data$date,
        precip = ensemble_means[, "precip"],
        tmean = ensemble_means[, "tmean"],
        tag = "sim"
      )
    )

    saveRDS(out, out_path)
  })





# ==============================================================================
# Compute correlations between simulated and observed data
# ==============================================================================

corr_data = parallel::mclapply(
  X = seq_len(nrow(station_meta)),
  mc.cores = n_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {

    if (i %% 100 == 0) message(i, " / ", nrow(station_meta))

    path = file.path(out_dir, paste0(station_meta$id[i], ".rds"))
    if (!file.exists(path)) {
      message("File ", path, " does not seem to exist")
      return(NULL)
    }

    data = readRDS(path)

    correlations = data[, .(corr = cor(precip, tmean)), by = "tag"]
    correlations$id = station_meta$id[i]

    correlations
  })
corr_data = rbindlist(corr_data)
corr_data = dcast(corr_data, id ~ tag, value.var = "corr")

summary(corr_data$obs)
summary(corr_data$sim)
cor(corr_data$obs, corr_data$sim)
