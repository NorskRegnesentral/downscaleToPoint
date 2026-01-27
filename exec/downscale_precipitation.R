library(data.table)
library(mgcv)
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
library(patchwork)
library(dplyr)

"And I should look into Schaake Shuffling!"

"
We could maybe use the Schaake Shuffle from neighbouring stations to improve the matching of the temperature/precipitation ensemble members! But we use quite few simulations per station, and that might cause problems. But it is worth writing about in the discussion!
"

# Define all necessary paths
# ------------------------------------------------------------------------------
data_dir = "/nr/samba/user/smvandeskog/projects/downscaleToPoint/data/"
model_dir = file.path(data_dir, "models", "precipitation")
image_dir = file.path(data_dir, "images", "precipitation")
local_fits_dir = file.path(model_dir, "local_fits")
cv_dir = file.path(data_dir, "cross-validation", "precipitation")

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

# Random seeds for reproducibility
set.seed(20260126)
base_seed = sample.int(1e8, 1)
seed_jump = sample.int(1e4, 1)

# Should we overwrite already existing results?
overwrite_local_models = FALSE
overwrite_evaluation = FALSE
overwrite_global_model = FALSE

# Thresholds for computing threshold weighted IQD scores during the cross-validation
threshold_probs = c(.8, .9, .95, .99, .995, .999)

# K is the number of neighbours to use for simulating precipitation at an unknown location.
# K_vals is the vector of all values of K we will test during the cross-validation experiment
K_vals = c(5, 10, 15, 20, 25, 30)

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
# Fit the global GAMs
# ==============================================================================

# Only fit the global model if we have not already done so
if (overwrite_global_model || !file.exists(global_fit_path)) {

  # Load all data from all available weather stations
  data = load_station_data(
    meta = station_meta,
    data_dir = data_dir,
    verbose = TRUE,
    rm_bad_flags = TRUE
  )

  # Add extra variables for the modelling
  data[, let(
    yday = yday(date),
    era_precip_bool = era_precip > 0,
    precip_bool = precip > 0,
    station_elevation = log(station_elevation + 1),
    era_log_precip = log(era_precip + 1)
  )]

  # Fit the intensity model
  # ------------------------------------------------------------------------------
  formula = precip ~
    s(era_log_precip) +
    s(station_elevation) +
    s(elevation_diff) +
    s(grid_elevation_sd) +
    s(era_tmean) +
    s(yday, bs = "cc") +
    s(lon, lat, bs = "sos")

  intensity_fit = bam(
    formula = formula,
    family = Gamma(link = "log"),
    data = data[precip > 0],
    discrete = TRUE,
    samfrac = .1,
    control = list(trace = TRUE)
  )

  # Estimate the gamma shape parameter using MASS
  intensity_fit$shape = MASS::gamma.shape(intensity_fit)$alpha

  # Remove unneccesary variables that take up a lot of memory
  unneccessary_vars = c("wt", "y", "prior.weights", "model", "offset", "weights",
    "residuals", "fitted.values", "linear.predictors")
  intensity_fit[unneccessary_vars] = NULL
  gc()

  # Fit the time-independent occurrence model
  # ------------------------------------------------------------------------------
  occurrence_formula = precip_bool ~
    era_precip_bool +
    s(era_log_precip) +
    s(station_elevation) +
    s(elevation_diff) +
    s(grid_elevation_sd) +
    s(era_tmean) +
    s(yday, bs = "cc") +
    s(lon, lat, bs = "sos")

  occurrence_fit = bam(
    formula = occurrence_formula,
    family = binomial(),
    data = data,
    discrete = TRUE,
    samfrac = .1,
    control = list(trace = TRUE)
  )

  occurrence_fit[unneccessary_vars] = NULL
  gc()

  # Fit the time-dependent occurrence models
  # ------------------------------------------------------------------------------
  data[, let(time_diff = c(as.integer(diff(date)), NA_integer_)), by = "id"]
  data = data[time_diff == 1]

  data = data[, .(
    precip_bool_change = (tail(precip_bool, -1) != head(precip_bool, -1)),
    precip_bool = head(precip_bool, -1),
    era_precip_bool = head(era_precip_bool, -1),
    next_era_precip_bool = tail(era_precip_bool, -1),
    era_log_precip = head(era_log_precip, -1),
    next_era_log_precip = tail(era_log_precip, -1),
    era_tmean = head(era_tmean, -1),
    next_era_tmean = tail(era_tmean, -1),
    yday = head(yday, -1),
    lon = lon[1],
    lat = lat[1],
    station_elevation = station_elevation[1],
    grid_elevation_sd = grid_elevation_sd[1],
    elevation_diff = elevation_diff[1]
  ), by = "id"]

  data[, let(
    era_log_precip_change = next_era_log_precip - era_log_precip,
    era_tmean_change = next_era_tmean - era_tmean
  )]

  # Formula for the occurrence models
  formula = precip_bool_change ~
    era_precip_bool +
    next_era_precip_bool +
    s(era_log_precip) +
    s(era_log_precip_change) +
    s(era_tmean) +
    s(era_tmean_change) +
    s(station_elevation) +
    s(elevation_diff) +
    s(grid_elevation_sd) +
    s(yday, bs = "cc") +
    s(lon, lat, bs = "sos")

  # Fit the dry -> wet model
  dry_to_wet_fit = bam(
    formula = formula,
    family = binomial(),
    data = data[precip_bool == FALSE],
    discrete = TRUE,
    samfrac = .1,
    control = list(trace = TRUE)
  )

  # Fit the wet -> wet model
  data$precip_bool_nochange = !data$precip_bool_change
  formula = update(formula, precip_bool_nochange ~ .)
  wet_to_wet_fit = bam(
    formula = formula,
    family = binomial(),
    data = data[precip_bool == TRUE],
    discrete = TRUE,
    samfrac = .1,
    control = list(trace = TRUE)
  )

  # Remove unneccesary variables that take up a lot of memory
  dry_to_wet_fit[unneccessary_vars] = NULL
  wet_to_wet_fit[unneccessary_vars] = NULL
  gc()

  # Save the results
  saveRDS(
    object = list(
      intensity = intensity_fit,
      occurrence = occurrence_fit,
      dry_to_wet = dry_to_wet_fit,
      wet_to_wet = wet_to_wet_fit
    ),
    file = global_fit_path
  )

  rm(data)
  gc()
}

# ==============================================================================
# Plot the global fit
# ==============================================================================

global_fit = readRDS(global_fit_path)

plot_data = lapply(
  seq_along(global_fit),
  function(i) {
    out = fast_mgcv_plot_data(
      global_fit[[i]],
      lon_range = c(-30, 50),
      lat_range = c(30, 75)
    )
    for (j in seq_along(out)) out[[j]]$model_name = names(global_fit)[i]
    out
  })
plot_data = unlist(plot_data, recursive = FALSE)

plots = lapply(global_fit, function(x) list())
for (i in seq_along(plot_data)) {
  plot = ggplot() +
    theme_light() +
    theme(text = element_text(size = 15))
  if (all(c("lon", "lat") %in% names(plot_data[[i]]))) {
    map = rnaturalearth::ne_countries(returnclass = "sf", scale = 110)
    plot = plot +
      geom_sf(data = map) +
      geom_raster(data = plot_data[[i]], aes(x = lon, y = lat, fill = value), alpha = .8) +
      geom_sf(data = map, fill = NA) +
      coord_sf(
        xlim = range(plot_data[[i]]$lon),
        ylim = range(plot_data[[i]]$lat),
        crs = st_crs(4326)
      ) +
      scale_fill_viridis_c() +
      labs(title = "Spatial effect", x = "", y = "", fill = "$\\eta$") +
      scale_x_continuous(
        breaks = seq(-20, 40, by = 20),
        labels = paste0("$", c(20, 0, 20, 40), "^\\circ$", c("W", "", "E", "E")),
        minor_breaks = seq(-10, 30, by = 10),
        expand = c(0, 0)
      ) +
      scale_y_continuous(
        breaks = seq(40, 70, by = 10),
        labels = paste0("$", seq(40, 70, by = 10), "^\\circ$N"),
        minor_breaks = seq(35, 65, by = 10),
        expand = c(0, 0)
      )
  } else {
    title_name = factor(
      names(plot_data)[i],
      levels = c(
        "era_tmean", "era_tmean_change",
        "era_log_precip", "era_log_precip_change",
        "station_elevation",
        "elevation_diff", "yday", "grid_elevation_sd"
      ),
      labels = c(
        "ERA5 temperature", "$\\Delta$ ERA5 temperature",
        "ERA5 precipitation", "$\\Delta$ ERA5 precipitation",
        "Station elevation",
        "Elevation difference", "Seasonal effect", "ERA5 elevation SD"
      )
    )
    x_name = factor(
      names(plot_data)[i],
      levels = c(
        "era_tmean", "era_tmean_change",
        "era_log_precip", "era_log_precip_change",
        "station_elevation", "elevation_diff", "yday", "grid_elevation_sd"
      ),
      labels = c(
        "$^\\circ$C", "$^\\circ$C",
        "mm/day", "mm/day",
        "m.a.s.l.", "m", "Day of the year", "m"
      )
    )
    plot = plot +
      geom_line(data = plot_data[[i]], aes(x = x, y = y)) +
      labs(title = title_name, y = "$\\eta$", x = x_name) +
      theme(axis.title.y = element_text(angle = 0, vjust = .5))
    if (title_name == "ERA5 precipitation") {
      plot = plot + scale_x_continuous(
        breaks = log(c(0, 2^seq(0, 20, by = 2)) + 1),
        labels = c(0, 2^seq(0, 20, by = 2))
      )
    } else if (title_name == "Station elevation") {
      plot = plot + scale_x_continuous(
        breaks = log(c(0, 2^seq(0, 20, by = 3)) + 1),
        labels = c(0, 2^seq(0, 20, by = 3))
      )
    }
  }
  model_name = plot_data[[i]]$model_name[1]
  plots[[model_name]][[names(plot_data)[i]]] = plot
}

plot_designs = c("
aaccee##
aacceegg
bbddffgg
bbddff##
",
"
aacceegg##
aacceeggii
bbddffhhii
bbddffhh##
"
)

for (i in seq_along(plots)) {
  if (length(plots[[i]]) == 7) {
    plot_design = plot_designs[1]
    width = 14
    height = 6
  } else {
    plot_design = plot_designs[2]
    width = 17
    height = 6
  }
  plot = patchwork::wrap_plots(plots[[i]], design = plot_design)
  model_name = names(plots)[i]
  plot_tikz(
    file = file.path(image_dir, paste0("global_precip_", model_name, "_model.pdf")),
    plot = plot,
    width = width,
    height = height
  )
}

# ==============================================================================
# Perform local modelling
# ==============================================================================

global_fit = readRDS(global_fit_path)

# Loop over all weather stations and fit local GAM/ARMA models
start_time = Sys.time()
fits = parallel::mclapply(
  X = seq_len(nrow(station_meta)),
  mc.preschedule = FALSE,
  mc.cores = n_cores,
  FUN = function(i) {

    out_path = file.path(local_fits_dir, paste0(station_meta$id[i], ".rds"))

    # Check if the fits have already been created
    if (!overwrite_local_models && file.exists(out_path)) return()

    # Load data from the station of interest
    data = load_station_data(
      meta = station_meta[i, ],
      data_dir = data_dir,
      rm_na = FALSE,
      rm_bad_flags = TRUE
    )

    # Add extra variables for the modelling
    data[, let(
      yday = yday(date),
      era_precip_bool = era_precip > 0,
      precip_bool = precip > 0,
      station_elevation = log(station_elevation + 1),
      era_log_precip = log(era_precip + 1)
    )]
    data[, let(
      day_count = as.integer(date) - as.integer(min(date)) + 1L
    ),
    by = "id"]

    occurrence_prob = mean(data$precip_bool, na.rm = TRUE)

    # Fit the local intensity GAM
    # ------------------------------------------------------------------------------

    # Compute offset terms from the global models
    data$intensity_offset = fast_mgcv_pred(global_fit$intensity, data)

    intensity_formula = precip ~
      offset(intensity_offset) +
      s(era_log_precip) +
      s(era_tmean) +
      s(yday, bs = "cc")

    intensity_fit = bam(
      formula = intensity_formula,
      family = Gamma(link = "log"),
      data = data[!is.na(precip) & precip > 0],
      discrete = TRUE,
      control = list(trace = FALSE)
    )

    # Estimate the gamma shape parameter using MASS
    intensity_fit$shape = MASS::gamma.shape(intensity_fit)$alpha

    # Remove unneccesary variables that take up a lot of memory
    unneccessary_vars = c("wt", "y", "prior.weights", "model", "offset", "weights",
                          "residuals", "fitted.values", "linear.predictors")
    intensity_fit[unneccessary_vars] = NULL

    # Fit the local time-independent occurrence GAM
    # ------------------------------------------------------------------------------
    data$occurrence_offset = fast_mgcv_pred(global_fit$occurrence, data)

    occurrence_formula = precip_bool ~
      era_precip_bool +
      offset(occurrence_offset) +
      s(era_log_precip) +
      s(era_tmean) +
      s(yday, bs = "cc")

    occurrence_fit = bam(
      formula = occurrence_formula,
      family = binomial(),
      data = data[!is.na(precip)],
      discrete = TRUE,
      control = list(trace = FALSE)
    )

    # Remove unneccesary variables that take up a lot of memory
    occurrence_fit[unneccessary_vars] = NULL

    # Fit the local ARMA model
    # ------------------------------------------------------------------------------

    # Compute the estimated rate parameter for each row in `data`
    intensity_linpred = fast_mgcv_pred(intensity_fit, data) + data$intensity_offset
    intensity_mean = intensity_fit$family$linkinv(intensity_linpred)
    intensity_rate = intensity_fit$shape / intensity_mean

    # Transform the observed precipitation amounts to approximately Gaussian
    # random variables, using the probability integral transform
    intensity_pit = pgamma(
      q = data[, ifelse(precip > 0, precip, NA_real_)],
      shape = intensity_fit$shape,
      rate = intensity_rate
    )
    gaussian_precip = qnorm(intensity_pit)

    # Since the fitted GAMs do not provide perfect model fits, it will sometimes
    # happen that an observation is equal to the epsilon-quantile or the (1 - epsilon)-quantile
    # of the fitted GAMs, where epsilon is a really small number. When this happens, the
    # corresponding transformed Gaussian variables might obtain values of -Inf or Inf, which
    # is problematic when we want to fit an ARMA model to the transformed Gaussian variables.
    # Therefore, we censor all variables that are outside the 99.9998% confidence interval of
    # a standardised Gaussian random variable
    gaussian_precip[gaussian_precip < qnorm(1e-6)] = qnorm(1e-6)
    gaussian_precip[gaussian_precip > qnorm(1 - 1e-6)] = qnorm(1 - 1e-6)

    # Fit an ARMA model to the approximately Gaussian random variables
    arma_fit = forecast::auto.arima(
      y = gaussian_precip,
      stationary = TRUE,
      seasonal = FALSE,
      allowmean = FALSE
    )
    arma_fit = list(coef = arma_fit$coef, sigma2 = arma_fit$sigma2, arma = arma_fit$arma)
    arma_fit$model = list()
    if (arma_fit$arma[1] > 0) arma_fit$model$ar = head(arma_fit$coef, arma_fit$arma[1])
    if (arma_fit$arma[2] > 0) arma_fit$model$ma = tail(arma_fit$coef, arma_fit$arma[2])

    # Fit the local Markov occurrence models
    # ------------------------------------------------------------------------------

    data = data[!is.na(precip)]
    data[, let(time_diff = c(as.integer(diff(date)), NA_integer_)), by = "id"]
    data = data[time_diff == 1]

    data = data[, .(
      precip_bool_change = (tail(precip_bool, -1) != head(precip_bool, -1)),
      precip_bool_nochange = (tail(precip_bool, -1) == head(precip_bool, -1)),
      precip_bool = head(precip_bool, -1),
      era_precip_bool = head(era_precip_bool, -1),
      next_era_precip_bool = tail(era_precip_bool, -1),
      era_log_precip = head(era_log_precip, -1),
      next_era_log_precip = tail(era_log_precip, -1),
      era_tmean = head(era_tmean, -1),
      next_era_tmean = tail(era_tmean, -1),
      yday = head(yday, -1),
      lon = lon[1],
      lat = lat[1],
      station_elevation = station_elevation[1],
      grid_elevation_sd = grid_elevation_sd[1],
      elevation_diff = elevation_diff[1]
    )]

    data[, let(
      era_log_precip_change = next_era_log_precip - era_log_precip,
      era_tmean_change = next_era_tmean - era_tmean
    )]

    data$wet_to_wet_offset = fast_mgcv_pred(global_fit$wet_to_wet, data)
    data$dry_to_wet_offset = fast_mgcv_pred(global_fit$dry_to_wet, data)

    dry_to_wet_formula = precip_bool_change ~
      offset(dry_to_wet_offset) +
      era_precip_bool +
      next_era_precip_bool +
      s(era_log_precip) +
      s(era_log_precip_change) +
      s(era_tmean) +
      s(era_tmean_change) +
      s(yday, bs = "cc")

    # Fit the dry -> wet model
    dry_to_wet_fit = bam(
      formula = dry_to_wet_formula,
      family = binomial(),
      data = data[precip_bool == FALSE],
      discrete = TRUE,
      control = list(trace = FALSE)
    )

    wet_to_wet_formula = precip_bool_nochange ~
      offset(wet_to_wet_offset) +
      era_precip_bool +
      next_era_precip_bool +
      s(era_log_precip) +
      s(era_log_precip_change) +
      s(era_tmean) +
      s(era_tmean_change) +
      s(yday, bs = "cc")

    # Fit the wet -> wet model
    wet_to_wet_fit = bam(
      formula = wet_to_wet_formula,
      family = binomial(),
      data = data[precip_bool == TRUE],
      discrete = TRUE,
      control = list(trace = FALSE)
    )

    # Remove unneccesary variables that take up a lot of memory
    dry_to_wet_fit[unneccessary_vars] = NULL
    wet_to_wet_fit[unneccessary_vars] = NULL

    # Save the local model fits
    # ------------------------------------------------------------------------------
    res = data.table(
      id = station_meta$id[i],
      occurrence_prob = occurrence_prob,
      intensity = list(intensity_fit),
      occurrence = list(occurrence_fit),
      intensity_arma = list(arma_fit),
      dry_to_wet = list(dry_to_wet_fit),
      wet_to_wet = list(wet_to_wet_fit)
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

# Check how many of the ARMA models that are just equal to white noise, and the
# general summary of p and q in the models
# ------------------------------------------------------------------------------
arma_models = list()
local_fit_files = list.files(local_fits_dir, full.names = TRUE)
pb = progress_bar(length(local_fit_files))
for (i in seq_along(local_fit_files)) {
  fit = readRDS(local_fit_files[i])
  arma_models[[i]] = fit$intensity_arma[[1]]$arma[1:2]
  pb$tick()
}
pb$terminate()
arma_models = do.call(rbind, arma_models)
quantile(arma_models[, 1], seq(0, 1, by = .05)) # Number of AR terms (p)
quantile(arma_models[, 2], seq(0, 1, by = .05)) # Number of MA terms (q)
quantile(apply(arma_models, 1, sum), seq(0, 1, by = .05)) # p + q

# ==============================================================================
# Evaluation
# ==============================================================================

# Simulate downscaled precipitation at all locations, and compare properties of
# simulations from different models with each other and with ERA
# ------------------------------------------------------------------------------

global_fit = readRDS(global_fit_path)

# Loop over all weather stations for all values of K, simulate data and
# compute all scoring functions of interest. This takes a lot of time
start_time = Sys.time()
for (K in K_vals) {
  out_dir = file.path(cv_dir, paste0(K, "_neighbours"))
  if (!dir.exists(out_dir)) dir.create(out_dir)
  parallel::mclapply(
    X = seq_len(nrow(station_meta)),
    mc.cores = n_cores,
    mc.preschedule = FALSE,
    FUN = function(i) {

      set.seed(base_seed + K + i * seed_jump)

      # Print our progress so far
      time_passed = Sys.time() - start_time
      message(
        "Starting on iter nr. ", i, " / ", nrow(station_meta), " with K=", K,
        ". Time passed: ", round(as.numeric(time_passed), 2), " ", attr(time_passed, "units")
      )

      out_path = file.path(out_dir, paste0(station_meta$id[i], ".rds"))
      if (!overwrite_evaluation && file.exists(out_path)) return(TRUE)

      # Compute distances to all other weather stations
      dists = geosphere::distHaversine(
        p1 = station_meta[i, c(lon, lat)],
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
        era_precip_bool = era_precip > 0,
        station_elevation = log(station_elevation + 1),
        precip_bool = precip > 0,
        era_log_precip = log(era_precip + 1),
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
      for (name in names(global_fit)) {
        if (grepl("wet", name)) {
          offsets[[name]] = fast_mgcv_pred(global_fit[[name]], data[-nrow(data)])
        } else {
          offsets[[name]] = fast_mgcv_pred(global_fit[[name]], data)
        }
      }

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

      # Simulate precipitation data using the full model,
      # which has time dependence for both intensity and occurrence
      sims$full = simulate_precip_with_donors(
        n_sims = n_sims,
        data = data,
        local_fits = local_fits,
        offsets = offsets,
        intensity_time_dep = TRUE,
        occurrence_time_dep = TRUE
      )

      # Check if any of the donor stations appear to be outliers and
      # Remove them if this is the case
      bad_full_donors = get_bad_donor_index(sims$full, mean, K)
      if (length(bad_full_donors) > 0) {
        sims$full = simulate_precip_with_donors(
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
        id = station_meta$id[i],
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
      res$rmse = list(c(era = era_rmse, sim_rmse))

      era_logrmse = rmse(log(data$precip + 1), log(data$era_precip + 1))
      sim_logrmse = sapply(mean_sim, function(x) rmse(log(data$precip + 1), log(x + 1)))
      res$logrmse = list(c(era = era_logrmse, sim_logrmse))

      tw_index = which(data$precip < quantile(data$precip, .95))
      era_twrmse = rmse(data$precip[tw_index], data$era_precip[tw_index])
      sim_twrmse = sapply(mean_sim, function(x) rmse(data$precip[tw_index], x[tw_index]))
      res$twrmse = list(c(era = era_twrmse, sim_twrmse))

      # Compute the probability of zero precipitation, using multiple different zero-thresholds
      sim_zero_probs = sapply(
        X = sims,
        FUN = function(x) sapply(zero_thresholds, function(y) mean(x <= y, na.rm = TRUE))
      )
      obs_zero_probs = sapply(zero_thresholds, function(x) mean(data$precip <= x, na.rm = TRUE))
      era_zero_probs = sapply(zero_thresholds, function(x) mean(data$era_precip <= x, na.rm = TRUE))
      res$zero_probs = list(cbind(obs = obs_zero_probs, era = era_zero_probs, sim_zero_probs))

      # Compare the ensemble median and ERA5 with the observed data, using MAE
      mae = function(x, y, ...) mean(abs(x - y), ...)
      median_sim = lapply(sims, matrixStats::rowMedians)
      sim_mae = sapply(median_sim, mae, x = data$precip)
      era_mae = mae(data$precip, data$era_precip)
      res$mae = list(c(era = era_mae, sim_mae))

      era_logmae = mae(log(data$precip + 1), log(data$era_precip + 1))
      sim_logmae = sapply(median_sim, function(x) mae(log(data$precip + 1), log(x + 1)))
      res$logmae = list(c(era = era_logmae, sim_logmae))

      tw_index = which(data$precip < quantile(data$precip, .95))
      era_twmae = mae(data$precip[tw_index], data$era_precip[tw_index])
      sim_twmae = sapply(median_sim, function(x) mae(data$precip[tw_index], x[tw_index]))
      res$twmae = list(c(era = era_twmae, sim_twmae))

      # Compare marginal distributions of all non-zero precipitation intensities
      era_iqd = iqd(data$era_precip, data$precip, rm_zero = TRUE, w = function(x) x > 2 & x < 40)
      sims_iqd = sapply(
        X = sims,
        FUN = function(x) {
          iqd(as.vector(x), y = data$precip, rm_zero = TRUE, w = function(x) x > 2 & x < 40)
        }
      )
      res$iqd = list(c(era = era_iqd, sims_iqd))

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
        function(j) tail(data_with_all_dates$precip, -j) - head(data_with_all_dates$precip, -j)
      )
      # Compute all n-day differences for ERA
      era_diffs = lapply(
        diff_lengths,
        function(j) tail(data_with_all_dates$era_precip, -j) - head(data_with_all_dates$era_precip, -j)
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

      # Compare marginal distributions for the weekly means and
      # standard deviations of precipitation data
      # from observations, ERA5 and simulated precipitation.
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
        weekly_data$mean$obs[j] = mean(data$precip[week_indices$index[[j]]])
        weekly_data$mean$era[j] = mean(data$era_precip[week_indices$index[[j]]])
        weekly_data$mean$sim[[j]] = lapply(
          sims,
          function(x) matrixStats::colMeans2(x[week_indices$index[[j]], , drop = FALSE]))
        weekly_data$sd$obs[j] = sd(data$precip[week_indices$index[[j]]])
        weekly_data$sd$era[j] = sd(data$era_precip[week_indices$index[[j]]])
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

      # Compare marginal distributions for the monthly means and
      # standard deviations of precipitation data
      # from observations, ERA5 and simulated precipitation.
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
        monthly_data$mean$obs[j] = mean(data$precip[month_indices$index[[j]]])
        monthly_data$mean$era[j] = mean(data$era_precip[month_indices$index[[j]]])
        monthly_data$mean$sim[[j]] = lapply(
          sims,
          function(x) matrixStats::colMeans2(x[month_indices$index[[j]], , drop = FALSE])
        )
        monthly_data$sd$obs[j] = sd(data$precip[month_indices$index[[j]]])
        monthly_data$sd$era[j] = sd(data$era_precip[month_indices$index[[j]]])
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

# Load all eval data from the cross-validation experiment
eval_files = list.files(cv_dir, full.names = TRUE, recursive = TRUE)
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

# Find the best value of K
# ------------------------------------------------------------------------------

# Which values of K and which models should we evaluate?
#chosen_Ks = c(15, 20, 25, 30)
chosen_Ks = c(5, 10, 15, 20, 25, 30)
data_types = c("full")

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
plot = bootstrap_data[K1 > K0][K0 <= 20][K0 > 10][, let(
    data_type = factor(data_type, levels = c("full"), labels = c("Full")),
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
  #facet_grid(data_type ~ K0) +
  facet_wrap(~K0, nrow = 1) +
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
  file = file.path(image_dir, "precip_K_scores.pdf"),
  plot = plot,
  width = 10,
  height = 5
)

# Find the best downscaling model for the best value of K
# ------------------------------------------------------------------------------

#chosen_K = 25
chosen_K = 20
data_types = c("era", "local", "full", "global")

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
  !(data_type0 == "global" & data_type1 == "era")
][
  !(data_type0 == "local" & data_type1 %in% c("era", "global"))
][
  !data_type0 == "full"
][, let(
  data_type1 = factor(
    data_type1,
    levels = rev(c("full", "local", "global", "era")),
    labels = rev(c("Full", "Local", "Global", "ERA5"))
  ),
  data_type0 = factor(
    data_type0,
    levels = rev(c("full", "local", "global", "era")),
    labels = rev(paste("$S_0$:", c("Full", "Local", "Global", "ERA5")))
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
  file = file.path(image_dir, "precip_scores.pdf"),
  plot = plot,
  width = 12,
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
  stat_summary_hex(
    fun = my_hex_func,
    data = plot_data,
    aes(x = lon, y = lat, z = scores),
    bins = 15
  ) +
  geom_sf(data = map, fill = NA) +
  scale_fill_scico(palette = "vik", limits = c(-1, 1), direction = -1) +
  coord_sf(
    xlim = st_bbox(plot_data)[c(1, 3)] + c(-.5, 1.5),
    ylim = st_bbox(plot_data)[c(2, 4)] + c(-1, 1),
    crs = st_crs(plot_data)
  ) +
  labs(x = "", y = "", fill = "Skill\nscore") +
  facet_wrap(~score_name, nrow = 4) +
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
  file = file.path(image_dir, "precip_map_scores.pdf"),
  #file = "Rplots.pdf",
  tex_engine = "lualatex",
  plot = plot,
  width = 11,
  height = 11
)

# Convert the plot to png, to reduce the size of the final paper
pdf_convert(
  in_path = file.path(image_dir, "precip_map_scores.pdf"),
  out_paths = file.path(image_dir, "precip_map_scores.png"),
  format = "png"
)

# Create a map plot for raw RMSE and MAE values
# ------------------------------------------------------------------------------

rmse_data = as.data.table(do.call(rbind, eval[K == chosen_K, rmse]))
mae_data = as.data.table(do.call(rbind, eval[K == chosen_K, mae]))
rmse_data = rmse_data[, .(full, era)][, let(id = eval[K == chosen_K, id], tag = "rmse")]
mae_data = mae_data[, .(full, era)][, let(id = eval[K == chosen_K, id], tag = "mae")]

plot_data = rbind(
  melt(rmse_data, id.vars = c("id", "tag")),
  melt(mae_data, id.vars = c("id", "tag"))
)
plot_data[, let(
  both = list(c(value[variable == "full"], value[variable == "era"]))
), by = c("id", "tag")]
plot_data = merge(
  plot_data,
  station_meta[, .(id, lon, lat, elev, elev_mean, elev_sd)],
  by = "id"
)

precip_data = load_station_data(
  meta = station_meta,
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

plot_data = merge(plot_data, precip_data[, .(id, precip, era_precip)], by = "id")
plot_data[, let(
  precip_diff = precip - era_precip,
  elev_diff = elev - elev_mean
)]

if (FALSE) {

  plots = list()
  plots[[1]] = ggplot(plot_data) +
    geom_point(aes(x = precip, y = value), alpha = .4) +
    facet_grid(variable ~ tag)
  plots[[2]] = ggplot(plot_data) +
    geom_point(aes(x = precip_diff, y = value), alpha = .4) +
    facet_grid(variable ~ tag)
  plots[[3]] = ggplot(plot_data) +
    geom_point(aes(x = elev_diff, y = value), alpha = .4) +
    facet_grid(variable ~ tag)
  plots[[4]] = ggplot(plot_data) +
    geom_point(aes(x = abs(elev_diff), y = value), alpha = .4) +
    facet_grid(variable ~ tag)
  plots[[5]] = ggplot(plot_data) +
    geom_point(aes(x = elev_sd, y = value), alpha = .4) +
    facet_grid(variable ~ tag)

  pdf("Rplots.pdf")
  for (plot in plots) print(plot)
  dev.off()

}


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
  levels = c("full", "era"),
  labels = c("Full", "ERA5")
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
        bins = 15
      ) +
      geom_sf(data = map, fill = NA) +
      labs(x = "", y = "",
           fill = t,
           title = paste0(v, ", ", t)) +
      scale_fill_viridis_c(
        option = if (t == "MAE") "D" else "D",
        limits = if (t == "MAE") c(0, 5) else c(2, 13),
        breaks = seq(0, 20, by = 2)
      )
  }
  plots[[length(plots) + 1]] = ggplot() +
    geom_sf(data = map) +
    stat_summary_hex(
      fun = skill_score_hex_func,
      data = dplyr::filter(plot_data, tag == t, variable == v),
      aes(x = lon, y = lat, z = both),
      bins = 15
    ) +
    geom_sf(data = map, fill = NA) +
    labs(x = "", y = "",
         fill = "Skill",
         title = paste0(t, " skill score")) +
    scale_fill_scico(
      palette = "vik",
      limits = c(-.6, .6),
      transform = pseudo_log_transform,
      breaks = c(-.5, -.2, 0, .2, .5),
      direction = -1,
      labels = paste0("$", c(-.5, -.2, 0, .2, .5), "$")
    )
}

plots[[length(plots) + 1]] = ggplot() +
  geom_sf(data = map) +
  stat_summary_hex(
    fun = mean,
    data = dplyr::filter(plot_data, tag == tag[1], variable == variable[1]),
    aes(x = lon, y = lat, z = precip),
    bins = 15
  ) +
  geom_sf(data = map, fill = NA) +
  labs(x = "", y = "",
       fill = "Precipitation",
       title = "Mean annual\nprecipitation") +
  scale_fill_viridis_c(breaks = c(500, 1000, 1500))

plots[[length(plots) + 1]] = ggplot() +
  geom_sf(data = map) +
  stat_summary_hex(
    fun = mean,
    data = dplyr::filter(plot_data, tag == tag[1], variable == variable[1]),
    aes(x = lon, y = lat, z = elev),
    bins = 15
  ) +
  geom_sf(data = map, fill = NA) +
  labs(x = "", y = "",
       fill = "Elevation",
       title = "Station elevation"
       ) +
  scale_fill_viridis_c()

for (i in seq_along(plots)) {
  plots[[i]] = plots[[i]] +
    coord_sf(
      xlim = st_bbox(plot_data)[c(1, 3)] + c(-.5, 1.5),
      ylim = st_bbox(plot_data)[c(2, 4)] + c(-1, 1),
      crs = st_crs(plot_data)
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
    ) +
    theme_light() +
    theme(
      strip.text = element_text(colour = "black", size = rel(.8)),
      text = element_text(size = 18),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      title = element_text(size = rel(.8)),
      legend.text = element_text(size = rel(.5)),
      legend.title = element_text(vjust = 1)
    ) +
    if (i %in% c(3, 6)) {
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
  file = file.path(image_dir, "precip_map_scores2.pdf"),
  tex_engine = "lualatex",
  plot = plot,
  width = 13,
  height = 7
)

# Convert the plot to png, to reduce the size of the final paper
pdf_convert(
  in_path = file.path(image_dir, "precip_map_scores2.pdf"),
  out_paths = file.path(image_dir, "precip_map_scores2.png"),
  format = "png"
)

# ==============================================================================
# Find good and bad stations, and plot time series' from them
# ==============================================================================

bad_ids = c("013250-99999", "111700-99999", "133780-99999", "112120-99999")
good_ids = c("152730-99999", "160540-99999", "152640-99999", "067200-99999")

time_series_data = lapply(
  X = c(bad_ids, good_ids),
  FUN = function(my_id) {

    i = which(station_meta$id == my_id)
    dists = geosphere::distHaversine(
      p1 = station_meta[i, c(lon, lat)],
      p2 = station_meta[, cbind(lon, lat)]
    )

    local_fits = list()
    for (index in order(dists)[-1]) {
      path = file.path(local_fits_dir, paste0(station_meta$id[index], ".rds"))
      if (!file.exists(path)) next
      fit = readRDS(path)
      fit$dist = dists[index]
      local_fits[[length(local_fits) + 1]] = fit
      if (length(local_fits) == chosen_K) break
    }
    local_fits = rbindlist(local_fits)

    # Load the data for the current weather station, and add necessary covariates
    data = load_station_data(
      meta = station_meta[i, ],
      data_dir = data_dir,
      rm_bad_flags = TRUE
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
    data[, let(day_count = as.integer(date) - as.integer(min(date)) + 1L), by = "id"]

    # Only choose the year with the most observations
    data[, let(n_per_year = .N), by = .(year(date))]
    data = data[n_per_year == max(n_per_year)]
    n_years = length(unique(year(data$date)))
    if (n_years > 1) data = data[year(date) == max(year(date))]

    # Add offsets from the global GAMs
    data$occurrence_offset = fast_mgcv_pred(global_fit$occurrence, data)
    data$intensity_offset = fast_mgcv_pred(global_fit$intensity, data)

    # Simulate precipitation data using the full model, including both local GAMs and ARMA models
    set.seed(1)
    sims = simulate_precip_with_donors(
      n_sims = n_sims,
      data = data,
      local_fits = local_fits,
      occurrence_offset = data$occurrence_offset,
      intensity_offset = data$intensity_offset,
      time_dep = TRUE
    )

    # Check if any of the donor stations appear to be outliers and
    # Remove them if this is the case
    bad_full_donors = get_bad_donor_index(sims, mean, chosen_K)
    if (length(bad_full_donors) > 0) {
      set.seed(1)
      sims = simulate_precip_with_donors(
        n_sims = n_sims,
        data = data,
        local_fits = local_fits[-bad_full_donors],
        occurrence_offset = data$occurrence_offset,
        intensity_offset = data$intensity_offset,
        time_dep = TRUE
      )
    }

    out = cbind(data[, .(date, precip, era_precip)], sims)
    out = melt(out, id.vars = "date")
    out[, let(
      type = factor(
        ifelse(grepl("^V\\d+$", as.character(variable)), "sim", as.character(variable)),
        levels = c("precip", "era_precip", "sim"),
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
  x = sub("&", "\\\\&", x)
  x
}
plot_data = merge(
  copy(time_series_data)[, let(value = cumsum(value)), by = .(year(date), variable, id)],
  station_meta[, .(
    id,
    name = format_name(name),
    country = format_name(country)
  )],
  all.x = TRUE)

plot_data[, let(id = factor(id, levels = unique(time_series_data$id)))]
plot_data = plot_data[order(id)]

plot_data[, let(tag = paste(name, country, year(date), sep = ", "))]
plot_data[, let(tag = factor(
  tag,
  levels = unique(tag),
  labels = paste0(seq_along(unique(tag)), ") ", unique(tag))
))]
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
    minor_breaks = NULL
  ) +
  labs(
    alpha = "",
    linetype = "",
    col = "",
    linewidth = "",
    y = "Cumulative precipitation for all dates\nwith non-missing observations [mm]",
    x = "Date"
  ) +
  facet_wrap(~tag, nrow = 2, scales = "free_y", dir = "v") +
  theme(
    strip.text = element_text(colour = "black"),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
    legend.position = "top"
  )

plot_tikz(
  file = file.path(image_dir, "precip_time_series.pdf"),
  tex_engine = "lualatex",
  plot = plot,
  width = 12,
  height = 5
)

# Convert pdf to png, to reduce the size of the final paper
pdf_convert(
  in_path = file.path(image_dir, "precip_time_series.pdf"),
  out_paths = file.path(image_dir, "precip_time_series.png"),
  format = "png"
)
