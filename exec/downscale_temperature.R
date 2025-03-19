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
library(downscaleToPoint)

# Define all necessary paths
# ------------------------------------------------------------------------------
data_dir = file.path(here::here(), "raw_data")
model_dir = file.path(data_dir, "models", "temperature")
result_dir = file.path(data_dir, "results")
image_dir = file.path(data_dir, "images")
local_fits_dir = file.path(model_dir, "local_fits")
cv_dir = file.path(data_dir, "cross-validation", "temperature")

meta_path = file.path(data_dir, "meta.rds")
global_fit_path = file.path(model_dir, "global.rds")

# Define other useful variables
# ------------------------------------------------------------------------------
n_cores = 8 # Number of cores to use for running code in parallel
n_sims = 150 # Number of ensembles to simulate during the downscaling
diff_lengths = c(1, 3, 7) # Which n-day-differences to evaluate in the cross-validation
B = 1000 # Number of bootstraps to use when bootstrapping

# K is the number of neighbours to use for simulating temperature at an unknown location.
# K_vals is the vector of all values of K we will test during the cross-validation experiment
K_vals = c(5, 10, 15, 20, 25, 30)

# This is a data.frame containing information about all the different scoring functions
# we want to evaluate in the cross-validation study
score_info = as.data.frame(t(do.call(cbind, list(
  c("rmse", "RMSE", 1),
  c("mae", "MAE", 1),
  c("iqd", "IQD", 1),
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
station_meta = readRDS(file.path(data_dir, "meta.rds"))

# Remove stations with few temperature observations
station_meta = station_meta[n_tmean > 200]

# ==============================================================================
# Fit the global GAM
# ==============================================================================

# Only fit the global model if we have not already done so
if (!file.exists(global_fit_path)) {

  # Load all data from all available weather stations
  data = load_station_data(
    meta = station_meta,
    data_dir = data_dir,
    verbose = TRUE,
    era_stats = TRUE
  )

  # Add extra variables for the modelling
  data[, let(
    yday = yday(date),
    era_log_precip = log(era_precip + 1)
  )]

  formula = tmean ~
    s(era_tmean) +
    s(era_log_precip) +
    s(station_elevation) +
    s(elevation_diff) +
    s(yday, bs = "cc") +
    s(lat, lon, bs = "sos", m = 2)

  # Fit the global model. We assume that era_tmean is already a pretty good predictor for tmean,
  # so we use this as a fixed offset in the model
  fit = bam(
    formula = formula,
    data = data,
    gc.level = 2,
    samfrac = .01,
    offset = era_tmean,
    discrete = TRUE,
    control = list(trace = TRUE)
  )

  # Remove unneccesary variables that take up a lot of memory
  unneccessary_vars = c("wt", "y", "prior.weights", "model", "offset", "weights",
    "residuals", "fitted.values", "linear.predictors")
  fit[unneccessary_vars] = NULL
  gc()

  # Save the results
  saveRDS(fit, global_fit_path)

  rm(data)
  gc()
}

# ==============================================================================
# Perform local modelling
# ==============================================================================

global_fit = readRDS(global_fit_path)

overwrite = FALSE

# Loop over all weather stations and fit local GAM/ARMA models
start_time = Sys.time()
fits = parallel::mclapply(
  X = seq_len(nrow(station_meta)),
  mc.preschedule = FALSE,
  mc.cores = n_cores,
  FUN = function(i) {

    out_path = file.path(local_fits_dir, paste0(station_meta$id[i], ".rds"))

    # Check if the fits have already been created
    if (!overwrite && file.exists(out_path)) return()

    # Load data from the station of interest
    data = load_station_data(
      meta = station_meta[i, ],
      data_dir = data_dir,
      rm_na = FALSE
    )

    # Add extra variables for the modelling
    data[, let(
      yday = yday(date),
      era_log_precip = log(era_precip + 1)
    )]

    # Fit the local GAM
    # ------------------------------------------------------------------------------

    # Compute the offset term from the global model
    data$tmean_offset = fast_mgcv_pred(global_fit, data) + data$era_tmean

    formula = tmean ~
      s(era_tmean) +
      era_log_precip +
      s(yday, bs = "cc")
    environment(formula) = NULL

    fit = bam(
      formula = formula,
      data = data,
      offset = tmean_offset,
      discrete = TRUE
    )

    # Remove unneccesary variables that take up a lot of memory
    unneccessary_vars = c("wt", "y", "prior.weights", "model", "offset", "weights",
                          "residuals", "fitted.values", "linear.predictors")
    fit[unneccessary_vars] = NULL

    # Fit the local ARMA model
    # ------------------------------------------------------------------------------

    # Transform the observed temperature means to approximately Gaussian
    # random variables, using the probability integral transform
    linpred = fast_mgcv_pred(fit, data) + data$tmean_offset
    tmean_pit = pnorm(data$tmean, mean = linpred, sd = fit$sig2)
    gaussian_tmean = qnorm(tmean_pit)

    # Since the fitted GAMs do not provide perfect model fits, it will sometimes
    # happen that an observation is equal to the epsilon-quantile or the (1 - epsilon)-quantile
    # of the fitted GAMs, where epsilon is a really small number. When this happens, the
    # corresponding transformed Gaussian variables might obtain values of -Inf or Inf, which
    # is problematic when we want to fit an ARMA model to the transformed Gaussian variables.
    # Therefore, we censor all variables that are outside the 99.9998% confidence interval of
    # a standardised Gaussian random variable
    gaussian_tmean[gaussian_tmean < qnorm(1e-6)] = qnorm(1e-6)
    gaussian_tmean[gaussian_tmean > qnorm(1 - 1e-6)] = qnorm(1 - 1e-6)

    # Fit an ARMA model to the approximately Gaussian random variables
    arma_fit = forecast::auto.arima(
      y = gaussian_tmean,
      stationary = TRUE,
      seasonal = FALSE,
      allowmean = FALSE
    )
    arma_fit = list(coef = arma_fit$coef, sigma2 = arma_fit$sigma2, arma = arma_fit$arma)
    arma_fit$model = list()
    if (arma_fit$arma[1] > 0) arma_fit$model$ar = head(arma_fit$coef, arma_fit$arma[1])
    if (arma_fit$arma[2] > 0) arma_fit$model$ma = tail(arma_fit$coef, arma_fit$arma[2])
    
    # Save the local model fits
    res = data.table(
      id = station_meta$id[i],
      marginal_fit = list(fit),
      arma_fit = list(arma_fit)
    )
    saveRDS(res, out_path)

    # Print the progress of the local model fitting
    time_passed = Sys.time() - start_time
    seconds_passed = as.numeric(time_passed, units = "secs")
    eta = Sys.time() + (nrow(station_meta) - i) * (seconds_passed / i)
    time_remaining = eta - Sys.time()
    message(
      "Done with iter nr. ", i, " / ", nrow(station_meta),
      ". Time passed: ", round(as.numeric(time_passed), 2), " ", attr(time_passed, "units"),
      ". ETA: ", round(as.numeric(time_remaining), 2), " ", attr(time_remaining, "units")
    )

  })

# ==============================================================================
# Evaluation
# ==============================================================================

# Define functions for simulating temperature from the downscaling models
# ------------------------------------------------------------------------------

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

# Simulate downscaled temperature at all locations, and compare properties of
# simulations from different models with each other and with ERA
# ------------------------------------------------------------------------------

global_fit = readRDS(global_fit_path)

# Loop over all weather stations for all values of K, simulate data and
# compute all scoring functions of interest. This takes a lot of time
overwrite = FALSE
start_time = Sys.time()
for (K in K_vals) {
  parallel::mclapply(
    X = seq_len(nrow(station_meta)),
    mc.cores = n_cores,
    mc.preschedule = FALSE,
    FUN = function(i) {

      # Print our progress so far
      time_passed = Sys.time() - start_time
      message(
        "Starting on iter nr. ", i, " / ", nrow(station_meta), " with K=", K,
        ". Time passed: ", round(as.numeric(time_passed), 2), " ", attr(time_passed, "units")
      )

      out_path = file.path(cv_dir, paste0(station_meta$id[i], "_", K, "-neighbours.rds"))
      if (!overwrite && file.exists(out_path)) return(TRUE)

      # Compute distances to all other weather stations
      dists = geosphere::distHaversine(
        p1 = station_meta[i, c(lon, lat)],
        p2 = station_meta[, cbind(lon, lat)]
      )

      # Locate and load the local models from the K nearest weather stations to weather station nr. i
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

      # Load the data for the current weather station, and add necessary covariates
      data = load_station_data(
        meta = station_meta[i, ],
        data_dir = data_dir
      )
      data[, let(
        yday = yday(date),
        year = year(date),
        era_log_precip = log(era_precip + 1)
      )]
      data[, let(day_count = as.integer(date) - as.integer(min(date)) + 1L), by = "id"]

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
      sims$global_deterministic = matrix(rep(data$tmean_offset, n_sims), nrow = nrow(data), ncol = n_sims)
      local_deterministic_donors = seq_len(K)
      if (length(bad_local_donors) > 0) local_deterministic_donors = local_deterministic_donors[-bad_local_donors]
      sims$local_deterministic = sapply(
        X = local_deterministic_donors,
        FUN = function(j) {
          fast_mgcv_pred(local_fits$marginal_fit[[j]], data) + data$tmean_offset
        })
      sims$local_deterministic = apply(sims$local_deterministic, 1, mean)
      sims$local_deterministic = matrix(rep(sims$local_deterministic, n_sims), nrow = nrow(data), ncol = n_sims)

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
        function(j) tail(data_with_all_dates$era_tmean, -j) - head(data_with_all_dates$era_tmean, -j)
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
      
      # Compare marginal distributions for the weekly means and standard deviations of temperature data
      # from observations, ERA5 and simulated temperature.
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
}

eval_files = list.files(cv_dir, full.names = TRUE)

eval = vector("list", length(eval_files))
pb = progress_bar(length(eval_files))
for (i in seq_along(eval_files)) {
  eval[[i]] = readRDS(eval_files[i])
  pb$tick()
}
pb$terminate()
eval = rbindlist(eval, fill = TRUE)

# Find the best value of K
# ------------------------------------------------------------------------------

# Which values of K and which models should we evaluate?
chosen_Ks = c(5, 10, 15, 20)
data_types = c("local_deterministic", "full")

# Compute bootstrapped confidence intervals for all the skill scores of interest
set.seed(1)
bootstrap_data = list()
for (i in seq_len(nrow(score_info))) {
  bootstrap_data[[i]] = bootstrap_skillscores(
    data = eval,
    score_name = score_info$name[i],
    data_types = data_types,
    K_vals = chosen_Ks,
    row_index = score_info$row_index[i],
    pairwise_data_types = FALSE
  )
  bootstrap_data[[i]]$score_name = score_info$shortname[i]
}
bootstrap_data = rbindlist(bootstrap_data)

# Plot all of the different skill scores, with 95% confidence intervals
plot = bootstrap_data[K1 > K0][, let(
  data_type = factor(
    data_type,
    levels = c("full", "local_deterministic"),
    labels = c("Full", "Local deterministic")
  ),
  K1 = factor(K1, levels = chosen_Ks, labels = paste0("$K = ", chosen_Ks, "$")),
  K0 = factor(K0, levels = chosen_Ks, labels = paste0("$S_0: K = ", chosen_Ks, "$")),
  score_name = factor(score_name, levels = score_info$shortname)
)] |>
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_point(
    aes(x = score_name, y = truth, col = K1, group = K1),
    position = position_dodge(.3),
    size = rel(.8)
  ) +
  geom_errorbar(
    aes(x = score_name, ymin = lower, ymax = upper, col = K1, group = K1),
    position = position_dodge(.3)
  ) +
  facet_grid(data_type~K0) +
  theme_light() +
  theme(
    strip.text = element_text(colour = "black", size = rel(1)),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
    axis.text.x = element_text(size = rel(1.1), angle = 70, vjust = .5),
    text = element_text(size = 15)
  ) +
  theme(legend.position = "top") +
  labs(x = "Scoring function", y = "$\\tilde S_{\\text{skill}}(S_1, S_0)$", col = "$S_1$:")

# Make a pretty plot, using the tikzDevice package
plot_tikz(
  file = file.path(image_dir, "temp_K_scores.pdf"),
  plot = plot,
  width = 11,
  height = 8
)

# Find the best downscaling model for the best value of K
# ------------------------------------------------------------------------------

chosen_K = 15
data_types = c("era", "local_deterministic", "full", "global_deterministic")

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

plot = bootstrap_data[
  !(data_type0 == "global_deterministic" & data_type1 == "era")
][
  !(data_type0 == "local_deterministic" & data_type1 %in% c("era", "global_deterministic"))
][
  !data_type0 == "full"
][, let(
  data_type1 = factor(
    data_type1,
    levels = rev(c("full", "local_deterministic", "global_deterministic", "era")),
    labels = rev(c("Full", "Local deterministic", "Global deterministic", "ERA5"))
  ),
  data_type0 = factor(
    data_type0,
    levels = rev(c("full", "local_deterministic", "global_deterministic", "era")),
    labels = rev(paste("$S_0$:", c("Full", "Local deterministic", "Global deterministic", "ERA5")))
  ),
  score_name = factor(score_name, levels = score_info$shortname)
)] |>
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_point(
    aes(x = score_name, y = truth, col = data_type1, group = data_type1),
    position = position_dodge(.3),
    size = rel(.8)
  ) +
  geom_errorbar(
    aes(x = score_name, ymin = lower, ymax = upper, col = data_type1, group = data_type1),
    position = position_dodge(.3)
  ) +
  facet_wrap(~data_type0, nrow = 1) +
  scale_y_continuous(breaks = seq(-10, 1, by = .2)) +
  theme_light() +
  theme(
    strip.text = element_text(colour = "black", size = rel(1)),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
    axis.text.x = element_text(size = rel(1.1), angle = 70, vjust = .5),
    text = element_text(size = 15)
  ) +
  theme(legend.position = "top") +
  labs(x = "Scoring function", y = "$\\tilde S_{\\text{skill}}(S_1, S_0)$", col = "$S_1$:")

plot_tikz(
  file = file.path(image_dir, "temp_scores.pdf"),
  plot = plot,
  width = 11,
  height = 5
)

# Create a map plot for skill scores between the full model and ERA5
# ------------------------------------------------------------------------------

data_types = c("full", "era")

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

score_data[, let(scores = lapply(seq_len(.N), function(i) c(full[i], era[i])))]

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
  stat_summary_hex(fun = my_hex_func, data = plot_data, aes(x = lon, y = lat, z = scores), bins = 15) +
  geom_sf(data = map, fill = NA) +
  scale_fill_scico(palette = "vik", limits = c(-1, 1)) +
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
    breaks = seq(-20, 40, by = 20),
    labels = paste0("$", c(20, 0, 20, 40), "^\\circ$", c("W", "", "E", "E")),
    minor_breaks = seq(-10, 30, by = 10)
  ) +
  scale_y_continuous(
    breaks = seq(40, 70, by = 10),
    labels = paste0("$", seq(40, 70, by = 10), "^\\circ$N"),
    minor_breaks = seq(35, 65, by = 10)
  )

plot_tikz(
  file = file.path(image_dir, "temp_map_scores.pdf"),
  tex_engine = "lualatex",
  plot = plot,
  width = 12,
  height = 8
)

# Convert the plot to png, to reduce the size of the final paper 
pdf_convert(
  in_path = file.path(image_dir, "temp_map_scores.pdf"),
  out_paths = file.path(image_dir, "temp_map_scores.png"),
  format = "png"
)

# ==============================================================================
# Find good and bad stations, and plot time series' from them
# ==============================================================================

bad_ids = c("135860-99999", "164600-99999", "062730-99999", "166820-99999")
good_ids = c("075910-99999", "163450-99999", "156150-99999", "109460-99999")

time_series_data = lapply(
  X = c(bad_ids, good_ids),
  FUN = function(my_id) {

    i = which(station_meta$id == my_id)
    dists = geosphere::distHaversine(
      p1 = station_meta[i, c(lon, lat)],
      p2 = station_meta[, cbind(lon, lat)]
    )

    # Locate and load the local models from the K nearest weather stations to weather station nr. i
    nearest_index = order(dists)[-1][seq_len(chosen_K)]
    local_fits = lapply(
      X = seq_along(nearest_index),
      FUN = function(j) {
        path = file.path(local_fits_dir, paste0(station_meta$id[nearest_index[j]], ".rds"))
        fit = readRDS(path)
        fit$dist = dists[nearest_index[j]]
        fit
      })
    local_fits = rbindlist(local_fits)

    # Load the data for the current weather station, and add necessary covariates
    data = load_station_data(
      meta = station_meta[i, ],
      data_dir = data_dir
    )
    data[, let(
      yday = yday(date),
      year = year(date),
      era_log_precip = log(era_precip + 1)
    )]
    data[, let(day_count = as.integer(date) - as.integer(min(date)) + 1L), by = "id"]

    # Only choose the year with the most observations
    data[, let(n_per_year = .N), by = .(year(date))]
    data = data[n_per_year == max(n_per_year)]
    n_years = length(unique(year(data$date)))
    if (n_years > 1) data = data[year(date) == max(year(date))]

    # Add the offset from the global GAM
    data$tmean_offset = data$era_tmean + fast_mgcv_pred(global_fit, data)

    # Simulate temperature data using the full model, including both local GAMs and ARMA models
    set.seed(1)
    sims = simulate_tmean_with_donors(
      n_sims = n_sims,
      data = data,
      local_fits = local_fits,
      offset = data$tmean_offset,
      use_arma = TRUE
    )

    # Check if any of the donor stations appear to be outliers and
    # Remove them if this is the case
    bad_full_donors = unique(c(
      get_bad_donor_index(sims, mean, chosen_K),
      get_bad_donor_index(sims, sd, chosen_K)
    ))
    if (length(bad_full_donors) > 0) {
      set.seed(1)
      sims = simulate_tmean_with_donors(
        n_sims = n_sims,
        data = data,
        local_fits = local_fits[-bad_full_donors, ],
        offset = data$tmean_offset,
        use_arma = TRUE
      )
    }

    out = cbind(data[, .(date, tmean, era_tmean)], sims)
    out = melt(out, id.vars = "date")
    out[, let(
      type = factor(
        ifelse(grepl("^V\\d+$", as.character(variable)), "sim", as.character(variable)),
        levels = c("tmean", "era_tmean", "sim"),
        labels = c("Observed", "ERA5", "Downscaled")
      ),
      id = my_id
    )]

    out
  }
)
time_series_data = rbindlist(time_series_data)

format_name = function(x) {
  x = tolower(x)
  x = sub("^(\\w)", "\\U\\1", x, perl = TRUE)
  x = gsub("([ -/\\(])(\\w)", "\\1\\U\\2", x, perl = TRUE)
  x = sub("(Ii+ *)$", "\\U\\1", x, perl = TRUE)
  x
}
plot_data = merge(
  time_series_data,
  station_meta[, .(
    id,
    name = format_name(name),
    country = format_name(country)
  )],
  all.x = TRUE)

plot_data[, let(id = factor(id, levels = unique(time_series_data$id)))]
plot_data = plot_data[order(id)]

plot_data[, let(tag = paste(name, country, year(date), sep = ", "))]
plot_data[, let(tag = factor(tag, levels = unique(tag), labels = paste0(seq_along(unique(tag)), ") ", unique(tag))))]
plot_data[, let(yday = yday(date))]
plot_data[, let(variable = factor(variable, levels = rev(levels(variable))))]

plot = ggplot(plot_data) +
  geom_line(aes(
    x = yday,
    y = value,
    alpha = type,
    linetype = type,
    group = variable,
    col = type,
    linewidth = type
  )) +
  scale_alpha_manual(values = c(1, 1, .1)) +
  scale_linetype_manual(values = c("solid", "solid", "solid")) +
  scale_color_manual(values = c("red", "blue", "black")) +
  scale_linewidth_manual(values = c(1, .8, .5)) +
  theme_light() +
  scale_x_continuous(
    breaks = cumsum(c(0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 32)),
    labels = c("1. Jan", "", "", "1. Apr", "", "", "1. Jul", "", "", "1. Oct", "", "", "31. Dec"),
    minor_breaks = NULL,
  ) +
  labs(
    alpha = "",
    linetype = "",
    col = "",
    linewidth = "",
    y = "Daily mean temperature [$^\\circ$C]",
    x = "Date"
  ) +
  facet_wrap(~tag, nrow = 2, scales = "free_y", dir = "v") +
  theme(
    strip.text = element_text(colour = "black"),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
    legend.position = "top"
  )

plot_tikz(
  file = file.path(image_dir, "temp_time_series.pdf"),
  tex_engine = "lualatex",
  plot = plot,
  width = 12,
  height = 5
)

# Convert pdf to png, to reduce the size of the final paper 
pdf_convert(
  in_path = file.path(image_dir, "temp_time_series.pdf"),
  out_paths = file.path(image_dir, "temp_time_series.png"),
  format = "png"
)

