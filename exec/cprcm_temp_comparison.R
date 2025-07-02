library(data.table)
library(mgcv)
library(here)
library(ggplot2)
library(forecast)
library(matrixStats)
library(geosphere)
library(sf)
library(rnaturalearth)
library(scico)
library(MASS)
library(parallel)
library(purrr)
library(patchwork)
library(downscaleToPoint)
library(ncdf4)

# Define all necessary paths
# ------------------------------------------------------------------------------
data_dir = file.path(here::here(), "raw_data")
model_dir = file.path(data_dir, "models", "temperature")
image_dir = file.path(data_dir, "images", "temperature_cprcm")
local_fits_dir = file.path(model_dir, "local_fits")
cv_dir = file.path(data_dir, "cross-validation", "temperature_cprcm")
cv_dir_old = file.path(data_dir, "cross-validation", "temperature")

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
B = 1000 # Number of bootstraps to use when bootstrapping
n_bins = 6 # Number of bins to use for hex-plots

# K is the number of neighbours to use for simulating temperature at an unknown location.
# K_vals is the vector of all values of K we will test during the cross-validation experiment
K_vals = c(5, 10, 15, 20, 25, 30)

# Thresholds for computing threshold weighted IQD scores during the cross-validation
#upper_threshold_probs = c(.8, .9, .95, .99)
#lower_threshold_probs = c(.2, .1, .05, .01)
threshold_probs = c(.001, .005, .01, .05, .1, .2, .5, .8, .9, .95, .99, .995, .999)

# This is a data.frame containing information about all the different scoring functions
# we want to evaluate in the cross-validation study
score_info = as.data.frame(t(do.call(cbind, list(
  c("rmse", "RMSE", 1),
  c("mae", "MAE", 1),
  c("iqd", "IQD", 1),
  c("quantile_score", "Q01", 3),
  c("quantile_score", "Q99", 11),
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

# Remove stations with few temperature observations
station_meta = station_meta[n_tmean > 200]

# ==============================================================================
# Functions
# ==============================================================================

simulate_tmean = function(n, marginal_fit, arma_fit, data, offset = 0) {
  # Compute the linear predictor
  linpred = fast_mgcv_pred(marginal_fit, data) + offset
  # Simulate Gaussian ARMA time series
  n_time = max(data$day_count) - min(data$day_count) + 1
  arma_sims = sapply(
    X = seq_len(n),
    FUN = function(i) {
      as.vector(arima.sim(n = n_time, model = arma_fit$model, sd = sqrt(arma_fit$sigma2)))
    })
  # Remove all ARMA simulations from dates where we have no observations, to ensure
  # that the observed and the simulated data correspond to each other
  arma_sims = arma_sims[data$day_count - min(data$day_count) + 1, , drop = FALSE]
  # Transform the ARMA simulations to have the same marginal distribution as the local fit
  res = arma_sims * marginal_fit$sig2 + linpred
  # Return the simulated data in a matrix with `n` columns
  matrix(res, nrow = nrow(arma_sims), ncol = n)
}

simulate_tmean_notime = function(n, fit, data, offset = 0) {
  # Compute the linear predictor
  linpred = fast_mgcv_pred(fit, data) + offset
  # Simulate the corresponding temperature means
  res = rnorm(n * length(linpred), mean = linpred, sd = fit$sig2)
  # Return the simulated data in a matrix with `n` columns
  matrix(res, nrow = length(linpred), ncol = n)
}

simulate_tmean_with_donors = function(n_sims,
                                      data,
                                      local_fits,
                                      offset,
                                      use_arma = TRUE) {
  K = nrow(local_fits)
  n_sims_per_local_fit = ceiling(n_sims / K)
  simulations = lapply(
    X = seq_len(K),
    FUN = function(i) {
      if (use_arma) {
        simulate_tmean(
          n = n_sims_per_local_fit,
          marginal_fit = local_fits$marginal_fit[[i]],
          arma_fit = local_fits$arma_fit[[i]],
          data = data,
          offset = offset
        )
      } else {
         simulate_tmean_notime(
          n = n_sims_per_local_fit,
          fit = local_fits$marginal_fit[[i]],
          data = data,
          offset = offset
        )
      }
    })
  do.call(cbind, simulations)
}

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

K = 10

overwrite = FALSE
start_time = Sys.time()

success = parallel::mclapply(
  X = seq_len(nrow(meta)),
  mc.cores = n_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {

    # Print our progress so far
    time_passed = Sys.time() - start_time
    message(
      "Starting on iter nr. ", i, " / ", nrow(meta), " with K=", K,
      ". Time passed: ", round(as.numeric(time_passed), 2), " ", attr(time_passed, "units")
    )

    out_path = file.path(cv_dir, paste0(meta$id[i], "_", K, "-neighbours.rds"))
    if (!overwrite && file.exists(out_path)) return(NULL)

    # Load the data for the current weather station, and add necessary covariates
    data = load_station_data(
      meta = meta[i, ],
      data_dir = data_dir
    )
    data[, let(
      yday = yday(date),
      year = year(date),
      station_elevation = log(station_elevation + 1),
      era_log_precip = log(era_precip + 1)
    )]

    cprcm_data = readRDS(file.path(data_dir, "cprcm", paste0(meta$id[i], ".rds")))

    cprcm_data = cprcm_data[cprcm_data$date %in% data$date]
    data = data[data$date %in% cprcm_data$date]
    stopifnot(nrow(cprcm_data) == nrow(data))
    if (nrow(data) < 300) return(NULL)
    data$cprcm_tmean = cprcm_data$temperature

    data[, let(day_count = as.integer(date) - as.integer(min(date)) + 1L), by = "id"]

    # Compute distances to all other weather stations
    dists = geosphere::distHaversine(
      p1 = meta[i, c(lon, lat)],
      p2 = station_meta[, cbind(lon, lat)]
    )

    # Locate and load the local models from the K nearest weather stations
    # to weather station nr. i
    nearest_index = order(dists)[-1][seq_len(K)]
    local_fits = lapply(
      X = seq_along(nearest_index),
      FUN = function(j) {
        path = file.path(local_fits_dir, paste0(station_meta$id[nearest_index[j]], ".rds"))
        fit = readRDS(path)
        fit$dist = dists[nearest_index[j]]
        fit
      })
    local_fits = rbindlist(local_fits)

    # Add the offset from the global GAM
    data$tmean_offset = data$era_tmean + fast_mgcv_pred(global_fit, data)

    # Preallocate a list that will hold all simulations from all our different models
    sims = list()

    # Simulate temperature data using the local GAMs, but not the local ARMA models
    set.seed(1)
    sims$local = simulate_tmean_with_donors(
      n_sims = n_sims,
      data = data,
      local_fits = local_fits,
      offset = data$tmean_offset,
      use_arma = FALSE
    )

    # Check if any of the donor stations appear to be outliers and
    # Remove them if this is the case
    bad_local_donors = unique(c(
      get_bad_donor_index(sims$local, mean, K),
      get_bad_donor_index(sims$local, sd, K)
    ))
    if (length(bad_local_donors) > 0) {
      set.seed(1)
      sims$local = simulate_tmean_with_donors(
        n_sims = n_sims,
        data = data,
        local_fits = local_fits[-bad_local_donors, ],
        offset = data$tmean_offset,
        use_arma = FALSE
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
    set.seed(1)
    sims$full = simulate_tmean_with_donors(
      n_sims = n_sims,
      data = data,
      local_fits = local_fits,
      offset = data$tmean_offset,
      use_arma = TRUE
    )

    # Check if any of the donor stations appear to be outliers and
    # Remove them if this is the case
    bad_full_donors = unique(c(
      get_bad_donor_index(sims$full, mean, K),
      get_bad_donor_index(sims$full, sd, K)
    ))
    if (length(bad_full_donors) > 0) {
      set.seed(1)
      sims$full = simulate_tmean_with_donors(
        n_sims = n_sims,
        data = data,
        local_fits = local_fits[-bad_full_donors, ],
        offset = data$tmean_offset,
        use_arma = TRUE
      )
    }

    # Simulate temperature means using the global model
    set.seed(1)
    sims$global = simulate_tmean_notime(
      n = n_sims,
      fit = global_fit,
      data = data,
      offset = data$era_tmean
    )

    # Start working on the object that will contain information about all relevant
    # scoring function values for the current weather station location
    res = data.table(
      id = meta$id[i],
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
    cprcm_rmse = rmse(data$tmean, data$cprcm_tmean)
    res$rmse = list(c(cprcm = cprcm_rmse, era = era_rmse, sim_rmse))

    # Compare the ensemble median and ERA5 with the observed data, using MAE
    mae = function(x, y, ...) mean(abs(x - y), ...)
    median_sim = lapply(sims, matrixStats::rowMedians)
    sim_mae = sapply(median_sim, mae, x = data$tmean)
    era_mae = mae(data$tmean, data$era_tmean)
    cprcm_mae = mae(data$cprcm_tmean, data$era_tmean)
    res$mae = list(c(cprcm = cprcm_mae, era = era_mae, sim_mae))

    # Compare marginal distributions of all daily temperature means
    era_iqd = iqd(data$era_tmean, data$tmean)
    cprcm_iqd = iqd(data$cprcm_tmean, data$tmean)
    sims_iqd = sapply(sims, function(x) iqd(as.vector(x), y = data$tmean))
    res$iqd = list(c(cprcm = cprcm_iqd, era = era_iqd, sims_iqd))

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
    cprcm_quantile_score = sapply(
      X = threshold_probs,
      FUN = quantile_score,
      pred = data$cprcm_tmean,
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
    res$quantile_score = list(cbind(
      cprcm = cprcm_quantile_score,
      era = era_quantile_score,
      sims_quantile_score
    ))

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
    cprcm_diffs = lapply(
      diff_lengths,
      function(j) {
        tail(data_with_all_dates$cprcm_tmean, -j) - head(data_with_all_dates$cprcm_tmean, -j)
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
    cprcm_diff_iqd = sapply(seq_along(cprcm_diffs), function(j) iqd(cprcm_diffs[[j]], obs_diffs[[j]]))
    sims_diff_iqd = sapply(
      X = sims_diffs,
      FUN = function(x) {
        sapply(
          X = seq_along(diff_lengths),
          FUN = function(k) iqd(as.vector(x[[k]]), y = obs_diffs[[k]])
        )
      })
    res$diff_iqd = list(cbind(cprcm = cprcm_diff_iqd, era = era_diff_iqd, sims_diff_iqd))

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
      mean = list(obs = NULL, era = NULL, cprcm = NULL, sim = list()),
      sd = list(obs = NULL, era = NULL, cprcm = NULL, sim = list())
    )
    # Compute all the means and standard deviations
    for (j in seq_len(nrow(week_indices))) {
      weekly_data$mean$obs[j] = mean(data$tmean[week_indices$index[[j]]])
      weekly_data$mean$era[j] = mean(data$era_tmean[week_indices$index[[j]]])
      weekly_data$mean$cprcm[j] = mean(data$cprcm_tmean[week_indices$index[[j]]])
      weekly_data$mean$sim[[j]] = lapply(
        sims,
        function(x) matrixStats::colMeans2(x[week_indices$index[[j]], , drop = FALSE]))
      weekly_data$sd$obs[j] = sd(data$tmean[week_indices$index[[j]]])
      weekly_data$sd$era[j] = sd(data$era_tmean[week_indices$index[[j]]])
      weekly_data$sd$cprcm[j] = sd(data$cprcm_tmean[week_indices$index[[j]]])
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
      weekly_mean_sim_iqd
    ))

    weekly_sd_era_iqd = iqd(weekly_data$sd$obs, weekly_data$sd$era)
    weekly_sd_cprcm_iqd = iqd(weekly_data$sd$obs, weekly_data$sd$cprcm)
    weekly_sd_sim_iqd = sapply(
      weekly_data$sd$sim,
      function(x) iqd(as.vector(x), x = weekly_data$sd$obs))
    res$weekly_sd_iqd = list(c(
      cprcm = weekly_sd_cprcm_iqd,
      era = weekly_sd_era_iqd,
      weekly_sd_sim_iqd
    ))

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
      mean = list(obs = NULL, era = NULL, cprcm = NULL, sim = list()),
      sd = list(obs = NULL, era = NULL, cprcm = NULL, sim = list())
    )
    # Compute all the means and standard deviations
    for (j in seq_len(nrow(month_indices))) {
      monthly_data$mean$obs[j] = mean(data$tmean[month_indices$index[[j]]])
      monthly_data$mean$era[j] = mean(data$era_tmean[month_indices$index[[j]]])
      monthly_data$mean$cprcm[j] = mean(data$cprcm_tmean[month_indices$index[[j]]])
      monthly_data$mean$sim[[j]] = lapply(
        sims,
        function(x) matrixStats::colMeans2(x[month_indices$index[[j]], , drop = FALSE])
      )
      monthly_data$sd$obs[j] = sd(data$tmean[month_indices$index[[j]]])
      monthly_data$sd$era[j] = sd(data$era_tmean[month_indices$index[[j]]])
      monthly_data$sd$cprcm[j] = sd(data$cprcm_tmean[month_indices$index[[j]]])
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
      monthly_mean_sim_iqd
    ))

    monthly_sd_era_iqd = iqd(monthly_data$sd$obs, monthly_data$sd$era)
    monthly_sd_cprcm_iqd = iqd(monthly_data$sd$obs, monthly_data$sd$cprcm)
    monthly_sd_sim_iqd = sapply(
      monthly_data$sd$sim,
      function(x) iqd(as.vector(x), x = monthly_data$sd$obs)
    )
    res$monthly_sd_iqd = list(c(
      cprcm = monthly_sd_cprcm_iqd,
      era = monthly_sd_era_iqd,
      monthly_sd_sim_iqd
    ))

    # Save the results
    saveRDS(res, out_path)
  })

# Load all eval data from the cross-validation experiment
eval_files = list.files(cv_dir, full.names = TRUE)
eval = vector("list", length(eval_files))
pb = progress_bar(length(eval_files))
for (i in seq_along(eval_files)) {
  eval[[i]] = readRDS(eval_files[i])
  pb$tick()
}
pb$terminate()
eval = rbindlist(eval, fill = TRUE)


chosen_K = 10
data_types = c("cprcm", "era", "local_deterministic", "full", "global_deterministic")

# Compute bootstrapped confidence intervals for all the skill scores of interest
set.seed(1)
bootstrap_data = list()
for (i in seq_len(nrow(score_info))) {
  bootstrap_data[[i]] = bootstrap_skillscores(
    data = eval,
    score_name = score_info$name[i],
    data_types = data_types,
    K_vals = chosen_K,
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
      levels = rev(c("full", "local_deterministic", "global_deterministic", "era", "cprcm")),
      labels = rev(c("Full", "Local deterministic", "Global deterministic", "ERA5", "CPRCM"))
    ),
    data_type1 = "Temperature",
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

precip_plot = readRDS(file.path(image_dir, "..", "precipitation_cprcm", "precip_scores_cprcm.rds"))

big_plot = patchwork::wrap_plots(
  plot,
  precip_plot,
  nrow = 1,
  axes = "collect",
  guides = "collect"
) &
  theme(legend.position = "top", legend.direction = "horizontal")

plot_tikz(
  file = file.path(image_dir, "..", "cprcm_scores.pdf"),
  plot = big_plot,
  width = 12,
  height = 5
)

# Create a map plot for skill scores between the full model and CPRCM
# ------------------------------------------------------------------------------

data_types = c("full", "cprcm")

score_data = list()
for (i in seq_len(nrow(score_info))) {
  score_data[[i]] = get_scores(
    data = eval,
    score_name = score_info$name[i],
    data_types = data_types,
    K_vals = chosen_K,
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
  facet_wrap(~score_name, nrow = 3) +
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
  file = file.path(image_dir, "temp_map_scores_cprcm.pdf"),
  tex_engine = "lualatex",
  plot = plot,
  width = 11,
  height = 7
)

# Convert the plot to png, to reduce the size of the final paper
pdf_convert(
  in_path = file.path(image_dir, "temp_map_scores_cprcm.pdf"),
  out_paths = file.path(image_dir, "temp_map_scores_cprcm.png"),
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
    K_vals = chosen_K,
    row_index = score_info$row_index[i]
  )
  score_data[[i]]$score_name = score_info$shortname[i]
}
score_data = rbindlist(score_data)
score_data = merge(score_data, station_meta[, .(id, lon, lat, elev, elev_mean)], by = "id")
score_data[, let(elev_diff = elev - elev_mean)]

score_data = dcast(score_data, ... ~ data_type, value.var = "value")
score_data[, let(scores = lapply(seq_len(.N), function(i) c(cprcm[i], full[i])))]

score_data = score_data[score_name %in% c("RMSE", "IQD", "MS")]

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
  labels = paste("Temperature", score_info$shortname)
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
    fun = function(x) {
      res = mean(abs(x))
      res[res < 10] = 10
      res
    },
    data = dplyr::filter(plot_data, score_name == score_name[1]) |>
      dplyr::mutate(tag = "Mean elevation difference"),
    aes(x = lon, y = lat, z = elev_diff),
    bins = n_bins
  ) +
  geom_sf(data = map, fill = NA) +
  labs(x = "", y = "", fill = "Elevation") +
  facet_wrap(~tag) +
  scale_fill_viridis_c(
    trans = "log",
    breaks = 10 * 2^(0:6),
    labels = c("$<10$", 10 * 2^(1:6))
  )

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

precip_plot = readRDS(file.path(
  image_dir, "..", "precipitation_cprcm", "precip_scores_cprcm_selected.rds"
))

big_plot = patchwork::wrap_plots(
  list(plot[[1]], plot[[2]], precip_plot[[1]], precip_plot[[2]]),
  nrow = 2,
  byrow = TRUE,
  widths = c(1, 3.05),
  guides = "collect"
)

plot_tikz(
  file = file.path(image_dir, "..", "map_scores_cprcm.pdf"),
  tex_engine = "lualatex",
  plot = big_plot,
  width = 15,
  height = 7
)

# Convert the plot to png, to reduce the size of the final paper
pdf_convert(
  in_path = file.path(image_dir, "..", "map_scores_cprcm.pdf"),
  out_paths = file.path(image_dir, "..", "map_scores_cprcm.png"),
  format = "png"
)

# Create a map plot for raw RMSE and MAE values
# ------------------------------------------------------------------------------

rmse_data = as.data.table(do.call(rbind, eval[K == chosen_K, rmse]))
mae_data = as.data.table(do.call(rbind, eval[K == chosen_K, mae]))
rmse_data = rmse_data[, .(full, cprcm)][, let(id = eval[K == chosen_K, id], tag = "rmse")]
mae_data = mae_data[, .(full, cprcm)][, let(id = eval[K == chosen_K, id], tag = "mae")]

plot_data = rbind(
  melt(rmse_data, id.vars = c("id", "tag")),
  melt(mae_data, id.vars = c("id", "tag"))
)
plot_data[, let(
  both = list(c(value[variable == "full"], value[variable == "cprcm"]))
), by = c("id", "tag")]
plot_data = merge(plot_data, station_meta[, .(id, lon, lat, elev, elev_mean)], by = "id")
plot_data[, let(elev_diff = elev - elev_mean)]

plot_data = st_as_sf(
  plot_data,
  coords = c("lon", "lat"),
  crs = st_crs(4326)
)
plot_data$lon = st_coordinates(plot_data)[, 1]
plot_data$lat = st_coordinates(plot_data)[, 2]
plot_data$tag = factor(plot_data$tag, levels = score_info$name, labels = score_info$shortname)
plot_data$variable = factor(
  plot_data$variable,
  levels = c("full", "cprcm"),
  labels = c("Full", "CPRCM")
)

map = rnaturalearth::ne_countries(returnclass = "sf", scale = 110)

pseudo_log_transform = scales::new_transform(
  name = "pseudo_log",
  transform = function(x) asinh(x * 10),
  inverse = function(x) sinh(x) / 10
)

# Function for computing the proper skill score value inside each hexagon
skill_score_hex_func = function(x) {
  s1 = mean(sapply(x, `[[`, 1))
  s0 = mean(sapply(x, `[[`, 2))
  res = skill_score(s1 = s1, s0 = s0)
  res = pmax(res, -1)
  res
}

plots = list()
for (t in unique(plot_data$tag)) {
  for (v in unique(plot_data$variable)) {
    plots[[length(plots) + 1]] = ggplot() +
      geom_sf(data = map) +
      stat_summary_hex(
        data = dplyr::filter(plot_data, tag == t, variable == v),
        aes(x = lon, y = lat, z = value),
        bins = n_bins
      ) +
      geom_sf(data = map, fill = NA) +
      labs(x = "", y = "", fill = t, title = paste0(v, ", ", t)) +
      scale_fill_viridis_c(
        option = if (t == "MAE") "D" else "D",
        limits = if (t == "MAE") c(.2, 10) else c(.2, 10),
        transform = "log",
        breaks = c(.5, 1, 2, 4, 8),
        labels = paste0("$", c(.5, 1, 2, 4, 8), "$")
      )
  }
  plots[[length(plots) + 1]] = ggplot() +
    geom_sf(data = map) +
    stat_summary_hex(
      fun = skill_score_hex_func,
      data = dplyr::filter(plot_data, tag == t, variable == v),
      aes(x = lon, y = lat, z = both),
      bins = n_bins
    ) +
    geom_sf(data = map, fill = NA) +
    labs(x = "", y = "", fill = "Skill", title = paste0(t, " skill score")) +
    scale_fill_scico(
      palette = "vik",
      limits = c(-1, 1),
      transform = pseudo_log_transform,
      direction = -1,
      breaks = c(-.5, -.2, 0, .2, .5)
    )
}

plots[[length(plots) + 1]] = ggplot() +
  geom_sf(data = map) +
  stat_summary_hex(
    fun = function(x) {
      res = mean(abs(x))
      res[res < 10] = 10
      res
    },
    data = dplyr::filter(plot_data, tag == tag[1], variable == variable[1]),
    aes(x = lon, y = lat, z = elev_diff),
    bins = n_bins
  ) +
  geom_sf(data = map, fill = NA) +
  labs(x = "", y = "", fill = "Elevation", title = "Mean elevation\ndifference") +
  scale_fill_viridis_c(
    trans = "log",
    breaks = 10 * 2^(0:6),
    labels = c("$<10$", 10 * 2^(1:6))
  )

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
    ) +
    if (i %in% c(3, 6, 7)) {
      theme(legend.text = element_text(size = rel(.5), angle = 70, vjust = .5))
    } else {
      theme(legend.text = element_text(size = rel(.5)))
    }
}

plot_design = "
##bbddff
aabbddff
aacceegg
##cceegg
"

plot = patchwork::wrap_plots(
  plots[c(7, 1, 4, 2, 5, 3, 6)],
  design = plot_design,
  guides = "collect"
)
plot = plot & theme(legend.position = "bottom")

plot_tikz(
  file = file.path(image_dir, "temp_map_scores_cprcm2.pdf"),
  tex_engine = "lualatex",
  plot = plot,
  width = 15,
  height = 7
)

# Convert the plot to png, to reduce the size of the final paper
pdf_convert(
  in_path = file.path(image_dir, "temp_map_scores_cprcm2.pdf"),
  out_paths = file.path(image_dir, "temp_map_scores_cprcm2.png"),
  format = "png"
)
