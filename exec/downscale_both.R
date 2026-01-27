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
out_dir = file.path(data_dir, "schaake_shuffle")
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
# Simulate temperature and precipitation from the full model, then
# use the Schaake Shuffle to pair members
# ==============================================================================

global_fits = readRDS(global_fit_paths[variables == "precipitation"])
global_fits$temperature = readRDS(global_fit_paths[variables == "temperature"])

parallel::mclapply(
  X = seq_len(nrow(station_meta)),
  mc.cores = n_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {

    if (FALSE) K = 1

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

    # Load data from all neighbouring stations, to learn their Schaake Shuffle ranks
    schaake_ranks = list()
    for (j in seq_along(local_fits$id)) {
      dt = load_station_data(
        meta = station_meta[id == local_fits$id[j]],
        data_dir = data_dir,
        rm_bad_flags = TRUE
      )
      dt = dt[order(tmean)]
      precip_ranks = rank(dt$precip, ties.method = "random")
      schaake_ranks[[j]] = cbind(seq_len(nrow(dt)), precip_ranks)
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

    sims$precip = simulate_precip_with_donors(
      n_sims = n_sims,
      data = data,
      local_fits = local_fits,
      offsets = offsets,
      occurrence_time_dep = TRUE,
      intensity_time_dep = TRUE
    )

    sapply(sims, dim)

    # Check if any of the donor stations appear to be outliers and
    # Remove them if this is the case
    bad_local_donors = unique(c(
      get_bad_donor_index(sims$local, mean, K),
      get_bad_donor_index(sims$local, sd, K)
    ))
    if (length(bad_local_donors) > 0) {
      sims$local = simulate_tmean_with_donors(
        n_sims = n_sims,
        data = data,
        local_fits = local_fits[-bad_local_donors, ],
        offset = data$tmean_offset,
        time_dep = FALSE
      )
    }

    # Compute temperature values for the local and the global deterministic downscaling models
    sims$global_deterministic = matrix(
      rep(data$tmean_offset, n_sims),
      nrow = nrow(data),
      ncol = n_sims
    )
    local_deterministic_donors = seq_len(K)
    if (length(bad_local_donors) > 0) {
      local_deterministic_donors = local_deterministic_donors[-bad_local_donors]
    }
    sims$local_deterministic = sapply(
      X = local_deterministic_donors,
      FUN = function(j) {
        fast_mgcv_pred(local_fits$marginal_fit[[j]], data) + data$tmean_offset
      })
    sims$local_deterministic = apply(sims$local_deterministic, 1, mean)
    sims$local_deterministic = matrix(
      rep(sims$local_deterministic, n_sims),
      nrow = nrow(data),
      ncol = n_sims
    )

    # Simulate temperature data using the full model, including both local GAMs and ARMA models
    sims$full = simulate_tmean_with_donors(
      n_sims = n_sims,
      data = data,
      local_fits = local_fits,
      offset = data$tmean_offset,
      time_dep = TRUE
    )

    # Check if any of the donor stations appear to be outliers and
    # Remove them if this is the case
    bad_full_donors = unique(c(
      get_bad_donor_index(sims$full, mean, K),
      get_bad_donor_index(sims$full, sd, K)
    ))
    if (length(bad_full_donors) > 0) {
      sims$full = simulate_tmean_with_donors(
        n_sims = n_sims,
        data = data,
        local_fits = local_fits[-bad_full_donors, ],
        offset = data$tmean_offset,
        time_dep = TRUE
      )
    }

    # Simulate temperature means using the global model
    sims$global = simulate_tmean_notime(
      n = n_sims,
      fit = global_fit,
      data = data,
      offset = data$era_tmean
    )

    # Start working on the object that will contain information about all relevant
    # scoring function values for the current weather station location
    res = data.table(
      id = station_meta$id[i],
      K = K,
      neighbour_ids = list(local_fits$id),
      neighbour_dists = list(local_fits$dist),
      n_obs = sum(!is.na(data$tmean)),
      n_bad_local_donors = length(bad_local_donors),
      n_bad_full_donors = length(bad_full_donors)
    )

    # Compare the ensemble means and ERA5 with the observed data, using RMSE
    rmse = function(x, y, ...) sqrt(mean((x - y)^2, ...))
    mean_sim = lapply(sims, matrixStats::rowMeans2)
    sim_rmse = sapply(mean_sim, rmse, x = data$tmean)
    era_rmse = rmse(data$tmean, data$era_tmean)
    res$rmse = list(c(era = era_rmse, sim_rmse))

    # Compare the ensemble median and ERA5 with the observed data, using MAE
    mae = function(x, y, ...) mean(abs(x - y), ...)
    median_sim = lapply(sims, matrixStats::rowMedians)
    sim_mae = sapply(median_sim, mae, x = data$tmean)
    era_mae = mae(data$tmean, data$era_tmean)
    res$mae = list(c(era = era_mae, sim_mae))

    # Compare marginal distributions of all daily temperature means
    era_iqd = iqd(data$era_tmean, data$tmean)
    sims_iqd = sapply(sims, function(x) iqd(as.vector(x), y = data$tmean))
    res$iqd = list(c(era = era_iqd, sims_iqd))

    # Compute quantile scores
    quantile_score = function(prob, pred, obs) {
      q = quantile(pred, probs = prob)
      mean(2 * (as.numeric(obs <= q) - prob) * (q - obs))
    }
    era_quantile_score = sapply(
      X = threshold_probs,
      FUN = quantile_score,
      pred = data$era_tmean,
      obs = data$tmean
    )
    sims_quantile_score = sapply(
      X = sims,
      FUN = function(sim) {
        sapply(
          X = threshold_probs,
          FUN = quantile_score,
          pred = sim,
          obs = data$tmean
        )
      }
    )
    res$quantile_score = list(cbind(era = era_quantile_score, sims_quantile_score))

    # Compare marginal distributions for all n-day differences, with n in `diff_lengths`
    # This is easiest to do if we first expand `data` so it contains one row for every
    # single date within `range(data$date)`
    data_with_all_dates = local({
      max_rows = diff(range(data$day_count)) + 1
      if (max_rows == nrow(data)) return(data)
      df2 = data.table(day_count = seq_len(max_rows)[-data$day_count])
      merge(data, df2, by = "day_count", all.y = TRUE, all.x = TRUE)[order(day_count)]
    })
    # Compute all n-day differences for the observed data
    obs_diffs = lapply(
      diff_lengths,
      function(j) tail(data_with_all_dates$tmean, -j) - head(data_with_all_dates$tmean, -j)
    )
    # Compute all n-day differences for ERA
    era_diffs = lapply(
      diff_lengths,
      function(j) {
        tail(data_with_all_dates$era_tmean, -j) - head(data_with_all_dates$era_tmean, -j)
      }
    )
    # Compute all n-day differences for each simulated ensemble member
    non_na_index = which(!is.na(data_with_all_dates$tmean))
    tmp = rep(NA_real_, nrow(data_with_all_dates))
    sims_diffs = lapply(
      sims, function(x) {
        diffs = lapply(
          X = seq_len(ncol(x)),
          FUN = function(j) {
            tmp[non_na_index] = x[, j]
            lapply(diff_lengths, function(k) tail(tmp, -k) - head(tmp, -k))
          })
        lapply(
          X = seq_along(diff_lengths),
          FUN = function(j) do.call(cbind, lapply(diffs, `[[`, j))
        )
      })

    # Compare the marginal distributions of the n-day differences using IQD
    era_diff_iqd = sapply(seq_along(era_diffs), function(j) iqd(era_diffs[[j]], obs_diffs[[j]]))
    sims_diff_iqd = sapply(
      X = sims_diffs,
      FUN = function(x) {
        sapply(
          X = seq_along(diff_lengths),
          FUN = function(k) iqd(as.vector(x[[k]]), y = obs_diffs[[k]])
        )
      })
    res$diff_iqd = list(cbind(era = era_diff_iqd, sims_diff_iqd))

    # Compare marginal distributions for the weekly means and standard deviations of
    # temperature data from observations, ERA5 and simulated temperature.
    #
    # week_indices is a data.table, describing which rows of `data` that contain data from which
    # week/year combinations
    week_indices = data[, .(
      week = week(date), year = year, index = seq_len(.N)
    )][, .(
      index = list(index)
    ), by = c("week", "year")]
    # Preallocate the object that will contain all weekly means and standard deviations
    weekly_data = list(
      mean = list(obs = NULL, era = NULL, sim = list()),
      sd = list(obs = NULL, era = NULL, sim = list())
    )
    # Compute all the means and standard deviations
    for (j in seq_len(nrow(week_indices))) {
      weekly_data$mean$obs[j] = mean(data$tmean[week_indices$index[[j]]])
      weekly_data$mean$era[j] = mean(data$era_tmean[week_indices$index[[j]]])
      weekly_data$mean$sim[[j]] = lapply(
        sims,
        function(x) matrixStats::colMeans2(x[week_indices$index[[j]], , drop = FALSE]))
      weekly_data$sd$obs[j] = sd(data$tmean[week_indices$index[[j]]])
      weekly_data$sd$era[j] = sd(data$era_tmean[week_indices$index[[j]]])
      weekly_data$sd$sim[[j]] = lapply(
        sims,
        function(x) matrixStats::colSds(x[week_indices$index[[j]], , drop = FALSE]))
    }
    for (j in seq_along(weekly_data)) {
      weekly_data[[j]]$sim = purrr::transpose(weekly_data[[j]]$sim)
      weekly_data[[j]]$sim = lapply(
        weekly_data[[j]]$sim,
        function(x) do.call(rbind, x)
      )
    }

    # Compare the weekly means and standard deviations using IQD
    weekly_mean_era_iqd = iqd(weekly_data$mean$obs, weekly_data$mean$era)
    weekly_mean_sim_iqd = sapply(
      weekly_data$mean$sim,
      function(x) iqd(as.vector(x), x = weekly_data$mean$obs))
    res$weekly_mean_iqd = list(c(era = weekly_mean_era_iqd, weekly_mean_sim_iqd))

    weekly_sd_era_iqd = iqd(weekly_data$sd$obs, weekly_data$sd$era)
    weekly_sd_sim_iqd = sapply(
      weekly_data$sd$sim,
      function(x) iqd(as.vector(x), x = weekly_data$sd$obs))
    res$weekly_sd_iqd = list(c(era = weekly_sd_era_iqd, weekly_sd_sim_iqd))

    # Compare marginal distributions for the monthly means and standard deviations of temperature data
    # from observations, ERA5 and simulated temperature.
    #
    # month_indices is a data.table, describing which rows of `data` that contain data from which
    # month/year combinations
    month_indices = data[, .(
      month = month(date), year = year, index = seq_len(.N)
    )][, .(
      index = list(index)
    ), by = c("month", "year")]
    # Preallocate the object that will contain all weekly means and standard deviations
    monthly_data = list(
      mean = list(obs = NULL, era = NULL, sim = list()),
      sd = list(obs = NULL, era = NULL, sim = list())
    )
    # Compute all the means and standard deviations
    for (j in seq_len(nrow(month_indices))) {
      monthly_data$mean$obs[j] = mean(data$tmean[month_indices$index[[j]]])
      monthly_data$mean$era[j] = mean(data$era_tmean[month_indices$index[[j]]])
      monthly_data$mean$sim[[j]] = lapply(
        sims,
        function(x) matrixStats::colMeans2(x[month_indices$index[[j]], , drop = FALSE])
      )
      monthly_data$sd$obs[j] = sd(data$tmean[month_indices$index[[j]]])
      monthly_data$sd$era[j] = sd(data$era_tmean[month_indices$index[[j]]])
      monthly_data$sd$sim[[j]] = lapply(
        sims,
        function(x) matrixStats::colSds(x[month_indices$index[[j]], , drop = FALSE])
      )
    }
    for (j in seq_along(monthly_data)) {
      monthly_data[[j]]$sim = purrr::transpose(monthly_data[[j]]$sim)
      monthly_data[[j]]$sim = lapply(
        monthly_data[[j]]$sim,
        function(x) do.call(rbind, x)
      )
    }

    # Compare the monthly means and standard deviations using IQD
    monthly_mean_era_iqd = iqd(monthly_data$mean$obs, monthly_data$mean$era)
    monthly_mean_sim_iqd = sapply(
      monthly_data$mean$sim,
      function(x) iqd(as.vector(x), x = monthly_data$mean$obs)
    )
    res$monthly_mean_iqd = list(c(era = monthly_mean_era_iqd, monthly_mean_sim_iqd))

    monthly_sd_era_iqd = iqd(monthly_data$sd$obs, monthly_data$sd$era)
    monthly_sd_sim_iqd = sapply(
      monthly_data$sd$sim,
      function(x) iqd(as.vector(x), x = monthly_data$sd$obs)
    )
    res$monthly_sd_iqd = list(c(era = monthly_sd_era_iqd, monthly_sd_sim_iqd))

    # Save the results
    saveRDS(res, out_path)
  })
