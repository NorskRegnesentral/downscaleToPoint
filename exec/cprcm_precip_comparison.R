library(data.table)
library(ggplot2)
library(matrixStats)
library(geosphere)
library(sf)
library(rnaturalearth)
library(scico)
library(parallel)
library(purrr)
library(patchwork)
library(downscaleToPoint)

# Define all necessary paths
# ------------------------------------------------------------------------------
data_dir = "/nr/samba/user/smvandeskog/projects/downscaleToPoint/data/"
model_dir = file.path(data_dir, "models", "precipitation")
image_dir = file.path(data_dir, "images", "precipitation_cprcm")
local_fits_dir = file.path(model_dir, "local_fits")
cv_dir = file.path(data_dir, "cross-validation", "precipitation_cprcm")
cv_dir_old = file.path(data_dir, "cross-validation", "precipitation")

if (!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)
if (!dir.exists(image_dir)) dir.create(image_dir, recursive = TRUE)
if (!dir.exists(local_fits_dir)) dir.create(local_fits_dir)
if (!dir.exists(cv_dir)) dir.create(cv_dir, recursive = TRUE)

meta_path = file.path(data_dir, "meta.rds")
global_fit_path = file.path(model_dir, "global.rds")

# Define other useful variables
# ------------------------------------------------------------------------------
n_cores = 8 # Number of cores to use for running code in parallel
n_sims = 150 # Number of ensembles to simulate during the downscaling
diff_lengths = c(1, 3, 7) # Which n-day-differences to evaluate in the cross-validation
zero_thresholds = c(0, .1, .5, 1) # Different precipitation thresholds for defining a day as dry
B = 1000 # Number of bootstraps to use when bootstrapping
n_bins = 6 # Number of bins to use for hex-plots

# Random seeds for reproducibility
set.seed(20260126)
base_seed = sample.int(1e8, 1)
seed_jump = sample.int(1e4, 1)

# Thresholds for computing threshold weighted IQD scores during the cross-validation
#upper_threshold_probs = c(.8, .9, .95, .99)
#lower_threshold_probs = c(.2, .1, .05, .01)
threshold_probs = c(.001, .005, .01, .05, .1, .2, .5, .8, .9, .95, .99, .995, .999)

# This is a data.frame containing information about all the different scoring functions
# we want to evaluate in the cross-validation study
score_info = as.data.frame(t(do.call(cbind, list(
  c("rmse", "RMSE", 1),
  c("mae", "MAE", 1),
  c("zero_probs_se", "ZP01", 2),
  c("zero_probs_se", "ZP1", 4),
  c("iqd", "IQD", 1),
  c("quantile_score", "Q95", 3),
  c("quantile_score", "Q99", 4),
  c("weekly_mean_iqd", "WM", 1),
  c("weekly_sd_iqd", "WS", 1),
  c("monthly_mean_iqd", "MM", 1),
  c("monthly_sd_iqd", "MS", 1),
  c("diff_iqd", "D1", 1),
  c("diff_iqd", "D3", 2),
  c("diff_iqd", "D7", 3)
))))
names(score_info) = c("name", "shortname", "row_index")
score_info$row_index = as.integer(score_info$row_index)

# Load the meta data
# ------------------------------------------------------------------------------
station_meta = readRDS(meta_path)

# Remove stations with few precipitation observations
station_meta = station_meta[n_good_precip_flag > 200]
station_meta = station_meta[n_unique_precip > 40]

# ==============================================================================
# CPRCM stuff
# ==============================================================================

cprcm_lon_range = c(-1.7374, 18.6974)
cprcm_lat_range = c(39.083, 50.3347)

meta = station_meta |>
  _[lon < cprcm_lon_range[2] - 1] |>
  _[lon > cprcm_lon_range[1] + 1] |>
  _[lat < cprcm_lat_range[2] - 1] |>
  _[lat > cprcm_lat_range[1] + 1]

# ==============================================================================
# Perform model evaluation
# ==============================================================================

global_fit = readRDS(global_fit_path)

K = 15

out_dir = file.path(cv_dir, paste0(K, "_neighbours"))
if (!dir.exists(out_dir)) dir.create(out_dir)

overwrite = FALSE
start_time = Sys.time()

success = parallel::mclapply(
  X = seq_len(nrow(meta)),
  mc.cores = n_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {

    set.seed(base_seed + i * seed_jump)

    # Print our progress so far
    time_passed = Sys.time() - start_time
    message(
      "Starting on iter nr. ", i, " / ", nrow(meta), " with K=", K,
      ". Time passed: ", round(as.numeric(time_passed), 2), " ", attr(time_passed, "units")
    )

    out_path = file.path(out_dir, paste0(meta$id[i], ".rds"))
    if (!overwrite && file.exists(out_path)) return(NULL)

    # Load the data for the current weather station, and add necessary covariates
    data = load_station_data(
      meta = meta[i, ],
      data_dir = data_dir,
      rm_bad_flags = TRUE,
      rm_na = FALSE
    )
    data[, let(
      yday = yday(date),
      year = year(date),
      month = month(date),
      era_precip_bool = era_precip > 0,
      station_elevation = log(station_elevation + 1),
      precip_bool = precip > 0,
      era_log_precip = log(era_precip + 1)
    )]

    cprcm_data = readRDS(file.path(data_dir, "cprcm", paste0(meta$id[i], ".rds")))
    cprcm_data = cprcm_data[date %in% data$date]
    data = data[date %in% cprcm_data$date]
    if (sum(!is.na(data$precip)) < 300) return(NULL)
    data$cprcm_precip = cprcm_data$precipitation

    data[, let(
      day_count = as.integer(date) - as.integer(min(date)) + 1L,
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
    for (name in names(global_fit)) {
      if (grepl("wet", name)) {
        offsets[[name]] = fast_mgcv_pred(global_fit[[name]], data[-nrow(data)])
      } else {
        offsets[[name]] = fast_mgcv_pred(global_fit[[name]], data)
      }
    }

    # Compute distances to all other weather stations
    dists = geosphere::distHaversine(
      p1 = meta[i, c(lon, lat)],
      p2 = station_meta[, cbind(lon, lat)]
    )

    # Locate and load the local models from the K nearest
    # weather stations to weather station nr. i
    local_fits = list()
    for (index in order(dists)[-1]) {
      path = file.path(local_fits_dir, paste0(station_meta$id[index], ".rds"))
      if (!file.exists(path)) next
      fit = readRDS(path)
      fit$dist = dists[index]
      local_fits[[length(local_fits) + 1]] = fit
      if (length(local_fits) == K) break
    }
    local_fits = rbindlist(local_fits)

    # Preallocate a list that will hold all simulations from all our different models
    sims = list()

    # Simulate precipitation data using the local GAMs, without time dependence
    sims$local = simulate_precip_with_donors(
      n_sims = n_sims,
      data = data,
      local_fits = local_fits,
      offsets = offsets,
      intensity_time_dep = FALSE,
      occurrence_time_dep = FALSE
    )

    # Check if any of the donor stations appear to be outliers and
    # Remove them if this is the case
    bad_local_donors = get_bad_donor_index(sims$local, mean, K)
    if (length(bad_local_donors) > 0) {
      sims$local = simulate_precip_with_donors(
        n_sims = n_sims,
        data = data,
        local_fits = local_fits[-bad_local_donors, ],
        offsets = offsets,
        intensity_time_dep = FALSE,
        occurrence_time_dep = FALSE
      )
    }

    # Simulate precipitation data using the *old* full model,
    # which only has temporal dependence in the intensity model
    sims$old_full = simulate_precip_with_donors(
      n_sims = n_sims,
      data = data,
      local_fits = local_fits,
      offsets = offsets,
      intensity_time_dep = TRUE,
      occurrence_time_dep = FALSE
    )

    # Check if any of the donor stations appear to be outliers and
    # Remove them if this is the case
    bad_old_full_donors = get_bad_donor_index(sims$old_full, mean, K)
    if (length(bad_old_full_donors) > 0) {
      sims$old_full = simulate_precip_with_donors(
        n_sims = n_sims,
        data = data,
        local_fits = local_fits[-bad_old_full_donors],
        offsets = offsets,
        intensity_time_dep = TRUE,
        occurrence_time_dep = FALSE
      )
    }

    # This is super hacky!
    local_fits$dry_to_wet = local_fits$dry_to_wet2
    local_fits$wet_to_wet = local_fits$wet_to_wet2

    # Simulate precipitation data using the full model,
    # which has time dependence for both intensity and occurrence
    sims$full2 = simulate_precip_with_donors(
      n_sims = n_sims,
      data = data,
      local_fits = local_fits,
      offsets = offsets,
      intensity_time_dep = TRUE,
      occurrence_time_dep = TRUE
    )

    # Check if any of the donor stations appear to be outliers and
    # Remove them if this is the case
    bad_full_donors = get_bad_donor_index(sims$full2, mean, K)
    if (length(bad_full_donors) > 0) {
      sims$full2 = simulate_precip_with_donors(
        n_sims = n_sims,
        data = data,
        local_fits = local_fits[-bad_full_donors],
        offsets = offsets,
        intensity_time_dep = TRUE,
        occurrence_time_dep = TRUE
      )
    }

    # This is even more hacky
    local_fits$dry_to_wet = local_fits$dry_to_wet1
    local_fits$wet_to_wet = local_fits$wet_to_wet1
    offsets$dry_to_wet = tail(offsets$occurrence, -1)
    offsets$wet_to_wet = tail(offsets$occurrence, -1)
    local_fits$occurrence_prob = rep(
      global_fit$occurrence$family$linkinv(offsets$occurrence[1]),
      K
    )

    sims$full1 = simulate_precip_with_donors(
      n_sims = n_sims,
      data = data,
      local_fits = local_fits,
      offsets = offsets,
      intensity_time_dep = TRUE,
      occurrence_time_dep = TRUE
    )

    # Check if any of the donor stations appear to be outliers and
    # Remove them if this is the case
    bad_full_donors = get_bad_donor_index(sims$full1, mean, K)
    if (length(bad_full_donors) > 0) {
      sims$full1 = simulate_precip_with_donors(
        n_sims = n_sims,
        data = data,
        local_fits = local_fits[-bad_full_donors],
        offsets = offsets,
        intensity_time_dep = TRUE,
        occurrence_time_dep = TRUE
      )
    }

    # Simulate precipitation data using the global model
    sims$global = local({
      occurrence = simulate_occurrence_notime(
        n = n_sims,
        fit = global_fit$occurrence,
        data = data
      )
      intensity = simulate_intensity_notime(
        n = n_sims,
        fit = global_fit$intensity,
        data = data
      )
      intensity * occurrence
    })

    # Start working on the object that will contain information about all relevant
    # scoring function values for the current weather station location
    res = data.table(
      id = meta$id[i],
      K = K,
      neighbour_ids = list(local_fits$id),
      neighbour_dists = list(local_fits$dist),
      n_obs = sum(!is.na(data$precip)),
      n_bad_local_donors = length(bad_local_donors),
      n_bad_full_donors = length(bad_full_donors),
      n_bad_old_full_donors = length(bad_old_full_donors)
    )

    # Start evaluating the different model simulations
    # ------------------------------------------------

    # First, remove NAs
    for (j in seq_along(sims)) sims[[j]] = sims[[j]][!is.na(data$precip), ]
    data = data[!is.na(precip)]

    # Compare the ensemble means and ERA5 with the observed data, using RMSE
    rmse = function(x, y, ...) sqrt(mean((x - y)^2, ...))
    mean_sim = lapply(sims, matrixStats::rowMeans2)
    sim_rmse = sapply(mean_sim, rmse, x = data$precip)
    era_rmse = rmse(data$precip, data$era_precip)
    cprcm_rmse = rmse(data$precip, data$cprcm_precip)
    res$rmse = list(c(
      cprcm = cprcm_rmse,
      era = era_rmse,
      sim_rmse
    ))

    era_logrmse = rmse(log(data$precip + 1), log(data$era_precip + 1))
    cprcm_logrmse = rmse(log(data$precip + 1), log(data$cprcm_precip + 1))
    sim_logrmse = sapply(mean_sim, function(x) rmse(log(data$precip + 1), log(x + 1)))
    res$logrmse = list(c(
      cprcm = cprcm_logrmse,
      era = era_logrmse,
      sim_logrmse
    ))

    tw_index = which(data$precip < quantile(data$precip, .95))
    era_twrmse = rmse(data$precip[tw_index], data$era_precip[tw_index])
    cprcm_twrmse = rmse(data$precip[tw_index], data$cprcm_precip[tw_index])
    sim_twrmse = sapply(mean_sim, function(x) rmse(data$precip[tw_index], x[tw_index]))
    res$twrmse = list(c(
      cprcm = cprcm_twrmse,
      era = era_twrmse,
      sim_twrmse
    ))

    # Compute the probability of zero precipitation, using multiple different zero-thresholds
    sim_zero_probs = sapply(
      X = sims,
      FUN = function(x) sapply(zero_thresholds, function(y) mean(x <= y, na.rm = TRUE))
    )
    obs_zero_probs = sapply(zero_thresholds, function(x) mean(data$precip <= x, na.rm = TRUE))
    era_zero_probs = sapply(zero_thresholds, function(x) mean(data$era_precip <= x, na.rm = TRUE))
    cprcm_zero_probs = sapply(
      zero_thresholds,
      function(x) mean(data$cprcm_precip <= x, na.rm = TRUE)
    )
    res$zero_probs = list(cbind(
      obs = obs_zero_probs,
      cprcm = cprcm_zero_probs,
      era = era_zero_probs,
      sim_zero_probs
    ))

    # Compare the ensemble median and ERA5 with the observed data, using MAE
    mae = function(x, y, ...) mean(abs(x - y), ...)
    median_sim = lapply(sims, matrixStats::rowMedians)
    sim_mae = sapply(median_sim, mae, x = data$precip)
    cprcm_mae = mae(data$precip, data$cprcm_precip)
    era_mae = mae(data$precip, data$era_precip)
    res$mae = list(c(
      cprcm = cprcm_mae,
      era = era_mae,
      sim_mae
    ))

    era_logmae = mae(log(data$precip + 1), log(data$era_precip + 1))
    cprcm_logmae = mae(log(data$precip + 1), log(data$cprcm_precip + 1))
    sim_logmae = sapply(median_sim, function(x) mae(log(data$precip + 1), log(x + 1)))
    res$logmae = list(c(
      cprcm = cprcm_logmae,
      era = era_logmae,
      sim_logmae))

    tw_index = which(data$precip < quantile(data$precip, .95))
    era_twmae = mae(data$precip[tw_index], data$era_precip[tw_index])
    cprcm_twmae = mae(data$precip[tw_index], data$cprcm_precip[tw_index])
    sim_twmae = sapply(median_sim, function(x) mae(data$precip[tw_index], x[tw_index]))
    res$twmae = list(c(
      cprcm = cprcm_twmae,
      era = era_twmae,
      sim_twmae))

    # Compare marginal distributions of all non-zero precipitation intensities
    era_iqd = iqd(data$era_precip, data$precip, rm_zero = TRUE, w = function(x) x > 2 & x < 40)
    cprcm_iqd = iqd(data$cprcm_precip, data$precip, rm_zero = TRUE, w = function(x) x > 2 & x < 40)
    sims_iqd = sapply(
      X = sims,
      FUN = function(x) {
        iqd(as.vector(x), y = data$precip, rm_zero = TRUE, w = function(x) x > 2 & x < 40)
      }
    )
    res$iqd = list(c(
      cprcm = cprcm_iqd,
      era = era_iqd,
      sims_iqd))

    # Compute quantile scores
    quantile_score = function(prob, pred, obs, rm_zero = TRUE, zero_threshold = 1) {
      if (rm_zero) {
        pred = pred[pred > zero_threshold]
        obs = obs[obs > zero_threshold]
      }
      q = quantile(pred, probs = prob)
      mean(2 * (as.numeric(obs <= q) - prob) * (q - obs))
    }
    era_quantile_score = sapply(
      X = threshold_probs,
      FUN = quantile_score,
      pred = data$era_precip,
      obs = data$precip
    )
    cprcm_quantile_score = sapply(
      X = threshold_probs,
      FUN = quantile_score,
      pred = data$cprcm_precip,
      obs = data$precip
    )
    sims_quantile_score = sapply(
      X = sims,
      FUN = function(sim) {
        sapply(
          X = threshold_probs,
          FUN = quantile_score,
          pred = sim,
          obs = data$precip
        )
      }
    )
    res$quantile_score = list(cbind(
      cprcm = cprcm_quantile_score,
      era = era_quantile_score,
      sims_quantile_score))

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
      function(j) tail(data_with_all_dates$precip, -j) - head(data_with_all_dates$precip, -j)
    )
    # Compute all n-day differences for ERA
    era_diffs = lapply(
      diff_lengths,
      function(j) {
        tail(data_with_all_dates$era_precip, -j) - head(data_with_all_dates$era_precip, -j)
      }
    )
    cprcm_diffs = lapply(
      diff_lengths,
      function(j) {
        tail(data_with_all_dates$cprcm_precip, -j) - head(data_with_all_dates$cprcm_precip, -j)
      }
    )
    # Compute all n-day differences for each simulated ensemble member
    non_na_index = which(!is.na(data_with_all_dates$precip))
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
    era_diff_iqd = sapply(
      seq_along(era_diffs),
      function(j) iqd(era_diffs[[j]], obs_diffs[[j]])
    )
    cprcm_diff_iqd = sapply(
      seq_along(cprcm_diffs),
      function(j) iqd(cprcm_diffs[[j]], obs_diffs[[j]])
    )
    sims_diff_iqd = sapply(
      X = sims_diffs,
      FUN = function(x) {
        sapply(
          X = seq_along(diff_lengths),
          FUN = function(k) iqd(as.vector(x[[k]]), y = obs_diffs[[k]])
        )
      })
    res$diff_iqd = list(cbind(
      cprcm = cprcm_diff_iqd,
      era = era_diff_iqd,
      sims_diff_iqd))

    # Compare marginal distributions for the weekly means and standard deviations
    # of precipitation data from observations, ERA5 and simulated precipitation.
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
      mean = list(obs = NULL, era = NULL, cprcm = NULL, sim = list()),
      sd = list(obs = NULL, era = NULL, cprcm = NULL, sim = list())
    )
    # Compute all the means and standard deviations
    for (j in seq_len(nrow(week_indices))) {
      weekly_data$mean$obs[j] = mean(data$precip[week_indices$index[[j]]])
      weekly_data$mean$era[j] = mean(data$era_precip[week_indices$index[[j]]])
      weekly_data$mean$cprcm[j] = mean(data$cprcm_precip[week_indices$index[[j]]])
      weekly_data$mean$sim[[j]] = lapply(
        sims,
        function(x) matrixStats::colMeans2(x[week_indices$index[[j]], , drop = FALSE]))
      weekly_data$sd$obs[j] = sd(data$precip[week_indices$index[[j]]])
      weekly_data$sd$era[j] = sd(data$era_precip[week_indices$index[[j]]])
      weekly_data$sd$cprcm[j] = sd(data$cprcm_precip[week_indices$index[[j]]])
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
    weekly_mean_cprcm_iqd = iqd(weekly_data$mean$obs, weekly_data$mean$cprcm)
    weekly_mean_sim_iqd = sapply(
      weekly_data$mean$sim,
      function(x) iqd(as.vector(x), x = weekly_data$mean$obs))
    res$weekly_mean_iqd = list(c(
      cprcm = weekly_mean_cprcm_iqd,
      era = weekly_mean_era_iqd,
      weekly_mean_sim_iqd))

    weekly_sd_era_iqd = iqd(weekly_data$sd$obs, weekly_data$sd$era)
    weekly_sd_cprcm_iqd = iqd(weekly_data$sd$obs, weekly_data$sd$cprcm)
    weekly_sd_sim_iqd = sapply(
      weekly_data$sd$sim,
      function(x) iqd(as.vector(x), x = weekly_data$sd$obs))
    res$weekly_sd_iqd = list(c(
      cprcm = weekly_sd_cprcm_iqd,
      era = weekly_sd_era_iqd,
      weekly_sd_sim_iqd))

    # Compare marginal distributions for the monthly means and standard deviations
    # of precipitation data from observations, ERA5 and simulated precipitation.
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
      mean = list(obs = NULL, era = NULL, cprcm = NULL, sim = list()),
      sd = list(obs = NULL, era = NULL, cprcm = NULL, sim = list())
    )
    # Compute all the means and standard deviations
    for (j in seq_len(nrow(month_indices))) {
      monthly_data$mean$obs[j] = mean(data$precip[month_indices$index[[j]]])
      monthly_data$mean$cprcm[j] = mean(data$cprcm_precip[month_indices$index[[j]]])
      monthly_data$mean$era[j] = mean(data$era_precip[month_indices$index[[j]]])
      monthly_data$mean$sim[[j]] = lapply(
        sims,
        function(x) matrixStats::colMeans2(x[month_indices$index[[j]], , drop = FALSE])
      )
      monthly_data$sd$obs[j] = sd(data$precip[month_indices$index[[j]]])
      monthly_data$sd$era[j] = sd(data$era_precip[month_indices$index[[j]]])
      monthly_data$sd$cprcm[j] = sd(data$cprcm_precip[month_indices$index[[j]]])
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
    monthly_mean_cprcm_iqd = iqd(monthly_data$mean$obs, monthly_data$mean$cprcm)
    monthly_mean_sim_iqd = sapply(
      monthly_data$mean$sim,
      function(x) iqd(as.vector(x), x = monthly_data$mean$obs)
    )
    res$monthly_mean_iqd = list(c(
      cprcm = monthly_mean_cprcm_iqd,
      era = monthly_mean_era_iqd,
      monthly_mean_sim_iqd))

    monthly_sd_era_iqd = iqd(monthly_data$sd$obs, monthly_data$sd$era)
    monthly_sd_cprcm_iqd = iqd(monthly_data$sd$obs, monthly_data$sd$cprcm)
    monthly_sd_sim_iqd = sapply(
      monthly_data$sd$sim,
      function(x) iqd(as.vector(x), x = monthly_data$sd$obs)
    )
    res$monthly_sd_iqd = list(c(
      cprcm = monthly_sd_cprcm_iqd,
      era = monthly_sd_era_iqd,
      monthly_sd_sim_iqd))

    # Save the results
    saveRDS(res, out_path)
  })

# Load all eval data from the cross-validation experiment
eval_files = list.files(out_dir, full.names = TRUE)
eval = vector("list", length(eval_files))
pb = progress_bar(length(eval_files))
for (i in seq_along(eval_files)) {
  eval[[i]] = readRDS(eval_files[i])
  pb$tick()
}
pb$terminate()
eval = rbindlist(eval, fill = TRUE)

# Add the zero prob squared error score to eval
eval$zero_probs_se = lapply(eval$zero_probs, function(x) (x[, -1] - x[, 1])^2)

data_types = c("era", "local", "full", "global", "cprcm")

# Compute bootstrapped confidence intervals for all the skill scores of interest
bootstrap_data = list()
for (i in seq_len(nrow(score_info))) {
  set.seed(base_seed + i * seed_jump)
  bootstrap_data[[i]] = bootstrap_skillscores(
    data = eval,
    score_name = score_info$name[i],
    data_types = data_types,
    K_vals = K,
    row_index = score_info$row_index[i]
  )
  bootstrap_data[[i]]$score_name = score_info$shortname[i]
}
bootstrap_data = rbindlist(bootstrap_data)

plot = bootstrap_data |>
  _[data_type1 == "cprcm"] |>
  _[data_type0 %in% c("full", "era", "cprcm")] |>
  _[, let(
    data_type0 = factor(
      data_type0,
      levels = rev(c("full", "local", "global", "era", "cprcm")),
      labels = rev(c("Full", "Local", "Global", "ERA5", "CPRCM"))
    ),
    data_type1 = "Precipitation",
    score_name = factor(score_name, levels = score_info$shortname)
  )] |>
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_point(
    aes(x = score_name, y = truth, col = data_type0, group = data_type0),
    position = position_dodge(.3),
    size = rel(.8)
  ) +
  geom_errorbar(
    aes(x = score_name, ymin = lower, ymax = upper, col = data_type0, group = data_type0),
    position = position_dodge(.3)
  ) +
  facet_wrap(~data_type1, nrow = 1) +
  scale_y_continuous(breaks = seq(-10, 1, by = .2), limits = c(-1.4, .88), expand = c(0, 0)) +
  theme_light() +
  theme(
    strip.text = element_text(colour = "black", size = rel(1)),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
    axis.text.x = element_text(size = rel(1.1), angle = 70, vjust = .5),
    text = element_text(size = 15)
  ) +
  theme(legend.position = "top") +
  labs(x = "Scoring function", y = "$\\tilde S_{\\text{skill}}(S_1, S_0)$", col = "$S_0$:")

# Save the ggplot so we can combine it with a similar temperature plot
# in cprcm_temp_comparison.R
saveRDS(plot, file.path(image_dir, "precip_scores_cprcm.rds"))

# ==============================================================================
# Create (skill) score scatter plots
# ==============================================================================

# we can look at skill on the y-axis against distance-to-sea, elevation, climatology
# we can also plot scores for the full model on the y-axis and for ERA on the x-axis

data_types = c("full", "cprcm")

score_data = list()
for (i in seq_len(nrow(score_info))) {
  score_data[[i]] = get_scores(
    data = eval,
    score_name = score_info$name[i],
    data_types = data_types,
    K_vals = K,
    row_index = score_info$row_index[i]
  )
  score_data[[i]]$score_name = score_info$shortname[i]
}
score_data = rbindlist(score_data)
score_data = merge(score_data, meta[, .(id, lon, lat, elev, elev_mean)], by = "id")
score_data[, let(elev_diff = elev - elev_mean)]

# Add distance to the sea
# ------------------------------------------------------------------------------
land_mask_path = file.path(data_dir, "era-land-mask.nc")
nc = ncdf4::nc_open(land_mask_path)
land_mask = ncdf4::ncvar_get(nc, "lsm")
land_mask_lon = ncdf4::ncvar_get(nc, "longitude")
land_mask_lon[land_mask_lon > 180] = land_mask_lon[land_mask_lon > 180] - 360
land_mask_lat = ncdf4::ncvar_get(nc, "latitude")
ncdf4::nc_close(nc)

get_dist_to_sea = function(lon, lat) {
  max_dist = .5
  while (TRUE) {
    lon_index = which(abs(land_mask_lon - lon) < max_dist)
    lat_index = which(abs(land_mask_lat - lat) < max_dist)
    land_mask_selection = land_mask[lon_index, lat_index]
    if (any(land_mask_selection < .1)) break
    max_dist = max_dist * 2
  }
  land_mask_coords = expand.grid(lon = land_mask_lon[lon_index], lat = land_mask_lat[lat_index])
  distances = geosphere::distHaversine(c(lon, lat), land_mask_coords)
  distances_to_sea_cells = distances[land_mask_selection < .3]
  min(distances_to_sea_cells)
}

score_data[, let(dist_to_sea = get_dist_to_sea(lon[1], lat[1])), by = c("lon", "lat")]

# Add precipitation climatologies
# ------------------------------------------------------------------------------

precip_means = parallel::mclapply(
  X = seq_len(nrow(meta)),
  mc.cores = 8,
  mc.preschedule = FALSE,
  FUN = function(i) {
    if (i %% 100 == 0) message(i, " / ", nrow(meta))
    data = load_station_data(
      meta = meta[i],
      data_dir = data_dir,
      verbose = FALSE,
      rm_na = FALSE
    )
    data.table(
      id = meta$id[i],
      precip_mean = mean(data$era_precip, na.rm = TRUE)
    )
  }
)
precip_means = rbindlist(precip_means)

score_data = merge(score_data, precip_means, by = "id")
score_data$precip_yearly_sum = score_data$precip_mean * 365

# Create the actual scatter plots
# ------------------------------------------------------------------------------

# Skill scatter plots with different variables along the x-axis
x_vars = c("elev", "elev_diff", "dist_to_sea", "precip_mean")
for (x_var in x_vars) {
  plot_data = score_data |>
    dcast(... ~ data_type, value.var = "value") |>
    _[, let(skill = skill_score(full, cprcm))]
  if (x_var == "elev_diff") plot_data[[x_var]] = abs(plot_data[[x_var]])
  plot = ggplot(plot_data) +
    geom_point(aes(x = !!sym(x_var), y = skill), alpha = .2) +
    geom_smooth(aes(x = !!sym(x_var), y = skill)) +
    labs(x = x_var, y = "Skill") +
    facet_wrap(~score_name, scales = "free_y", ncol = 5) +
    theme_light() +
    lims(y = c(-1, 1)) +
    theme(
      strip.text = element_text(colour = "black"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")
    )
  png(
    filename = file.path(image_dir, paste0("cprcm-skill-vs-", x_var, ".png")),
    width = 10, height = 6, units = "in", res = 150
  )
  print(plot)
  dev.off()
}

# full model score scatter plots with different variables along the x-axis
x_vars = c("elev", "elev_diff", "dist_to_sea", "precip_mean")
for (x_var in x_vars) {
  plot_data = score_data |>
    dcast(... ~ data_type, value.var = "value") |>
    _[, let(upper = quantile(cprcm, .998)), by = "score_name"] |>
    _[cprcm <= upper]
  if (x_var == "elev_diff") plot_data[[x_var]] = abs(plot_data[[x_var]])
  plot = ggplot(plot_data) +
    geom_point(aes(x = !!sym(x_var), y = cprcm), alpha = .2) +
    geom_smooth(aes(x = !!sym(x_var), y = cprcm)) +
    labs(x = x_var, y = "Skill") +
    facet_wrap(~score_name, scales = "free_y", ncol = 5) +
    theme_light() +
    theme(
      strip.text = element_text(colour = "black"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")
    )
  png(
    filename = file.path(image_dir, paste0("cprcm-score-vs-", x_var, ".png")),
    width = 10, height = 6, units = "in", res = 150
  )
  print(plot)
  dev.off()
}


# Create a map plot for skill scores between the full model and ERA5
# ------------------------------------------------------------------------------

data_types = c("full", "cprcm")

score_data = list()
for (i in seq_len(nrow(score_info))) {
  score_data[[i]] = get_scores(
    data = eval,
    score_name = score_info$name[i],
    data_types = data_types,
    K_vals = K,
    row_index = score_info$row_index[i]
  )
  score_data[[i]]$score_name = score_info$shortname[i]
}
score_data = rbindlist(score_data)
score_data = merge(score_data, station_meta[, .(id, lon, lat)], by = "id")

score_data = dcast(score_data, ... ~ data_type, value.var = "value")

score_data[, let(scores = lapply(seq_len(.N), function(i) c(cprcm[i], full[i])))]

my_hex_func = function(x) {
  s1 = mean(sapply(x, `[[`, 1))
  s0 = mean(sapply(x, `[[`, 2))
  res = skill_score(s1 = s1, s0 = s0)
  res = max(res, -1)
  res
}

plot_data = st_as_sf(
  score_data,
  coords = c("lon", "lat"),
  crs = st_crs(4326)
)
plot_data$lon = st_coordinates(plot_data)[, 1]
plot_data$lat = st_coordinates(plot_data)[, 2]
plot_data$score_name = factor(plot_data$score_name, levels = score_info$shortname)

map = rnaturalearth::ne_countries(returnclass = "sf", scale = 110)

plot = ggplot() +
  geom_sf(data = map) +
  stat_summary_hex(
    fun = my_hex_func,
    data = plot_data,
    aes(x = lon, y = lat, z = scores),
    bins = n_bins
  ) +
  geom_sf(data = map, fill = NA) +
  scale_fill_scico(palette = "vik", limits = c(-1, 1), direction = -1) +
  coord_sf(
    xlim = st_bbox(plot_data)[c(1, 3)] + c(-.5, 1.5),
    ylim = st_bbox(plot_data)[c(2, 4)] + c(-1, 1),
    crs = st_crs(plot_data)
  ) +
  labs(x = "", y = "", fill = "Skill\nscore") +
  facet_wrap(~score_name, ncol = 4) +
  theme_light() +
  theme(
    strip.text = element_text(colour = "black", size = rel(.8)),
    text = element_text(size = 18),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")
  ) +
  scale_x_continuous(
    breaks = seq(0, 20, by = 5),
    labels = paste0("$", seq(0, 20, by = 5), "^\\circ$", c("", "E", "E", "E", "E"))
  ) +
  scale_y_continuous(
    breaks = seq(40, 50, by = 2),
    labels = paste0("$", seq(40, 50, by = 2), "^\\circ$N")
  )

plot_tikz(
  file = file.path(image_dir, "precip_map_scores_cprcm.pdf"),
  tex_engine = "lualatex",
  plot = plot,
  width = 11,
  height = 9
)

# Convert the plot to png, to reduce the size of the final paper
pdf_convert(
  in_path = file.path(image_dir, "precip_map_scores_cprcm.pdf"),
  out_paths = file.path(image_dir, "precip_map_scores_cprcm.png"),
  format = "png"
)

# Create a map with some selected scoring functions
# ------------------------------------------------------------------------------

data_types = c("full", "cprcm")

score_data = list()
for (i in seq_len(nrow(score_info))) {
  score_data[[i]] = get_scores(
    data = eval,
    score_name = score_info$name[i],
    data_types = data_types,
    K_vals = K,
    row_index = score_info$row_index[i]
  )
  score_data[[i]]$score_name = score_info$shortname[i]
}
score_data = rbindlist(score_data)
score_data = merge(score_data, station_meta[, .(id, lon, lat, elev, elev_mean)], by = "id")
score_data[, let(elev_diff = elev - elev_mean)]

precip_data = load_station_data(
  meta = meta,
  data_dir = data_dir,
  verbose = TRUE,
  rm_bad_flags = TRUE
)
precip_data = precip_data[, .(
  precip = mean(precip),
  era_precip = mean(era_precip)
), by = .(yday(date), id)]
precip_data = precip_data[, .(
  precip = sum(precip),
  era_precip = sum(era_precip)
), by = "id"]

score_data = merge(score_data, precip_data[, .(id, precip, era_precip)], by = "id")

score_data = dcast(score_data, ... ~ data_type, value.var = "value")
score_data[, let(scores = lapply(seq_len(.N), function(i) c(cprcm[i], full[i])))]

score_data = score_data[score_name %in% c("RMSE", "ZP1", "D1")]

plot_data = st_as_sf(
  score_data,
  coords = c("lon", "lat"),
  crs = st_crs(4326)
)
plot_data$lon = st_coordinates(plot_data)[, 1]
plot_data$lat = st_coordinates(plot_data)[, 2]
plot_data$score_name = factor(
  plot_data$score_name,
  levels = score_info$shortname,
  labels = paste("Precipitation", score_info$shortname)
)

my_hex_func = function(x) {
  s1 = mean(sapply(x, `[[`, 1))
  s0 = mean(sapply(x, `[[`, 2))
  res = skill_score(s1 = s1, s0 = s0)
  res = max(res, -1)
  res
}

map = rnaturalearth::ne_countries(returnclass = "sf", scale = 110)

plots = list()
plots[[1]] = ggplot() +
  geom_sf(data = map) +
  stat_summary_hex(
    fun = mean,
    data = dplyr::filter(plot_data, score_name == score_name[1]) |>
   dplyr::mutate(tag = "Mean annual precipitation"),
    aes(x = lon, y = lat, z = precip),
    bins = n_bins
  ) +
  geom_sf(data = map, fill = NA) +
  scale_fill_viridis_c() +
  labs(x = "", y = "", fill = "Precipitation") +
  facet_wrap(~tag)

plots[[2]] = ggplot() +
  geom_sf(data = map) +
  stat_summary_hex(
    fun = my_hex_func,
    data = plot_data,
    aes(x = lon, y = lat, z = scores),
    bins = n_bins
  ) +
  geom_sf(data = map, fill = NA) +
  scale_fill_scico(palette = "vik", limits = c(-1, 1), direction = -1) +
  facet_wrap(~score_name, nrow = 1) +
  guides(y = "none") +
  labs(x = "", y = "", fill = "Skill\nscore")

for (i in seq_along(plots)) {
  plots[[i]] = plots[[i]] +
    coord_sf(
      xlim = st_bbox(plot_data)[c(1, 3)] + c(-.5, 1.5),
      ylim = st_bbox(plot_data)[c(2, 4)] + c(-1, 1),
      crs = st_crs(plot_data)
    ) +
    scale_x_continuous(
      breaks = seq(0, 20, by = 5),
      labels = paste0("$", seq(0, 20, by = 5), "^\\circ$", c("", "E", "E", "E", "E"))
    ) +
    scale_y_continuous(
      breaks = seq(40, 50, by = 2),
      labels = paste0("$", seq(40, 50, by = 2), "^\\circ$N")
    ) +
    theme_light() +
    theme(
      strip.text = element_text(colour = "black", size = rel(.8)),
      text = element_text(size = 18),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      title = element_text(size = rel(.8)),
      legend.title = element_text(vjust = 1)
    )
}

plot = patchwork::wrap_plots(plots, nrow = 1, widths = c(1, 3), guides = "collect")

# Save the ggplot so we can combine it with a similar temperature plot
# in cprcm_temp_comparison.R
saveRDS(plot, file.path(image_dir, "precip_scores_cprcm_selected.rds"))
