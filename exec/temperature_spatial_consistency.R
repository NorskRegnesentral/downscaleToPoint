library(data.table)
library(mgcv)
library(ggplot2)
library(forecast)
library(matrixStats)
library(geosphere)
library(sf)
library(rnaturalearth)
library(scico)
library(parallel)
library(purrr)
library(patchwork)
library(abind)
library(downscaleToPoint)

# Define all necessary paths
# ------------------------------------------------------------------------------
data_dir = "/nr/samba/user/smvandeskog/projects/downscaleToPoint/data/"
model_dir = file.path(data_dir, "models", "temperature")
image_dir = file.path(data_dir, "images", "temperature")
local_fits_dir = file.path(model_dir, "local_fits")
out_dir = file.path(data_dir, "spatial_consistency", "temperature")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

meta_path = file.path(data_dir, "meta.rds")
global_fit_path = file.path(model_dir, "global.rds")

# Define other useful variables
# ------------------------------------------------------------------------------
n_cores = 8 # Number of cores to use for running code in parallel
n_sims = 150 # Number of ensembles to simulate during the downscaling

K = 10

n_neighbour_min = 8

overwrite = FALSE

# Thresholds for computing threshold weighted IQD scores during the cross-validation
threshold_probs = c(.9, .95, .99)

# Random seeds for reproducibility
set.seed(20260129)
base_seed = sample.int(1e8, 1)
seed_jump = sample.int(1e4, 1)

# Load the meta data
# ------------------------------------------------------------------------------
station_meta = readRDS(meta_path)

# Remove stations with few temperature observations
station_meta = station_meta[n_tmean > 200]

# ==============================================================================
# Evaluation
# ==============================================================================

global_fit = readRDS(global_fit_path)

start_time = Sys.time()
parallel::mclapply(
  X = seq_len(nrow(station_meta)),
  mc.cores = n_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {

    set.seed(base_seed + K + i * seed_jump)

    # Print our progress so far
    time_passed = Sys.time() - start_time
    message(
      "Starting on iter nr. ", i, " / ", nrow(station_meta),
      ". Time passed: ", round(as.numeric(time_passed), 2), " ", attr(time_passed, "units")
    )

    out_path = file.path(out_dir, paste0(station_meta$id[i], ".rds"))
    if (!overwrite && file.exists(out_path)) return(TRUE)

    # Compute distances to all other weather stations
    dists = geosphere::distHaversine(
      p1 = station_meta[i, c(lon, lat)],
      p2 = station_meta[, cbind(lon, lat)]
    )

    # Locate all stations within a distance of 100 km
    neighbouring_ids = station_meta$id[dists <= 100e3]
    n_neighbouring_ids = length(neighbouring_ids)

    # Load data for all these stations
    obs = load_station_data(
      meta = station_meta[id %in% neighbouring_ids],
      data_dir = data_dir
    )

    # Only keep data from dates where at least `n_neighbour_min` of the stations
    # actually have any data
    obs[, let(n_obs_per_date = .N), by = c("date")]
    obs = obs[n_obs_per_date >= n_neighbour_min]
    if (nrow(obs) == 0) return(FALSE)

    obs_ids = unique(obs$id)
    obs_dates = sort(unique(obs$date))

    # Fill in extra rows with NAs, so obs has the same number of rows/dates for all ids
    obs = merge(
      obs,
      data.table(
        id = rep(obs_ids, each = length(obs_dates)),
        date = rep(obs_dates, length(obs_ids))
      ),
      all = TRUE
    )

    # Add the necessary covariates for performing downscaling
    obs[, let(
      yday = yday(date),
      year = year(date),
      station_elevation = log(station_elevation + 1),
      era_log_precip = log(era_precip + 1)
    )]
    obs[, let(day_count = as.integer(date) - as.integer(min(date)) + 1L)]

    # Simulate temperature at all of the neighbouring IDs
    sims = list()
    for (j in seq_along(obs_ids)) {
      # Compute distances to all other weather stations
      dists = geosphere::distHaversine(
        p1 = station_meta[id == obs_ids[j], c(lon, lat)],
        p2 = station_meta[, cbind(lon, lat)]
      )

      # Locate and load the local models from the K nearest weather stations
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

      # Compute the offset from the global GAM
      tmean_offset = obs[
        id == obs_ids[j],
        era_tmean + fast_mgcv_pred(global_fit, .SD)
      ]

      sim = simulate_tmean_with_donors(
        n_sims = n_sims,
        data = obs[id == obs_ids[j]],
        local_fits = local_fits,
        offset = tmean_offset,
        time_dep = TRUE
      )

      # Check if any of the donor stations appear to be outliers and
      # Remove them if this is the case
      na_rows = apply(sim, 1, function(x) any(is.na(x)))
      bad_donors = unique(c(
        get_bad_donor_index(sim[!na_rows, , drop = FALSE], mean, K),
        get_bad_donor_index(sim[!na_rows, , drop = FALSE], sd, K)
      ))
      if (length(bad_donors) > 0) {
        sim = simulate_tmean_with_donors(
          n_sims = n_sims,
          data = obs[id == obs_ids[j]],
          local_fits = local_fits[-bad_donors, ],
          offset = tmean_offset,
          time_dep = TRUE
        )
      }

      sims[[obs_ids[j]]] = sim[, seq_len(n_sims)]
    }
    sims = do.call(abind::abind, list(sims, along = 3))

    # Compute different stats over the entire domain
    sim_stats = list(
      mean = apply(sims, 2, matrixStats::rowMeans2, na.rm = TRUE),
      median = apply(sims, 2, matrixStats::rowMedians, na.rm = TRUE),
      sd = apply(sims, 2, matrixStats::rowSds, na.rm = TRUE),
      min = apply(sims, 2, matrixStats::rowMins, na.rm = TRUE),
      max = apply(sims, 2, matrixStats::rowMaxs, na.rm = TRUE)
    )

    obs_stats = obs[, .(
      mean = mean(tmean, na.rm = TRUE),
      median = median(tmean, na.rm = TRUE),
      sd = sd(tmean, na.rm = TRUE),
      min = min(tmean, na.rm = TRUE),
      max = max(tmean, na.rm = TRUE),
      n_stations = sum(!is.na(tmean))
    ), by = "date"][order(date)]

    era_stats = obs[, .(
      mean = mean(era_tmean, na.rm = TRUE),
      median = median(era_tmean, na.rm = TRUE),
      sd = sd(era_tmean, na.rm = TRUE),
      min = min(era_tmean, na.rm = TRUE),
      max = max(era_tmean, na.rm = TRUE)
    ), by = "date"][order(date)]

    # Create ranking histograms and estimate covarate probabilities
    res = list()
    for (name in names(sim_stats)) {
      obs_rank = cbind(obs_stats[[name]], sim_stats[[name]]) |>
        matrixStats::rowRanks(ties.method = "random") |>
        _[, 1]
      rank_mean = mean((obs_rank - 1) / n_sims)
      rank_sd = sd((obs_rank - 1) / n_sims)
      sim_quantiles = matrixStats::rowQuantiles(sim_stats[[name]], probs = c(.025, .05, .95, .975))

      res[[name]] = data.table(
        id = station_meta$id[i],
        stat = name,
        n_obs = nrow(obs_stats),
        rank_mean = rank_mean,
        rank_sd = rank_sd,
        coverage_90 = mean((obs_stats[[name]] >= sim_quantiles[, 2])
                           & (obs_stats[[name]] <= sim_quantiles[, 3])),
        coverage_95 = mean((obs_stats[[name]] >= sim_quantiles[, 1])
                           & (obs_stats[[name]] <= sim_quantiles[, 4]))
      )

      # Compare the ensemble means and ERA5 with the observed data, using RMSE
      rmse = function(x, y, ...) sqrt(mean((x - y)^2, ...))
      mean_sim = matrixStats::rowMeans2(sim_stats[[name]])
      sim_rmse = rmse(obs_stats[[name]], mean_sim)
      era_rmse = rmse(obs_stats[[name]], era_stats[[name]])
      res[[name]]$rmse = list(c(era = era_rmse, sim = sim_rmse))

      # Compare the ensemble median and ERA5 with the observed data, using MAE
      mae = function(x, y, ...) mean(abs(x - y), ...)
      median_sim = matrixStats::rowMedians(sim_stats[[name]])
      sim_mae = mae(obs_stats[[name]], median_sim)
      era_mae = mae(obs_stats[[name]], era_stats[[name]])
      res[[name]]$mae = list(c(era = era_mae, sim = sim_mae))

      # Compare marginal distributions of all daily temperature means
      sims_iqd = iqd(as.vector(obs_stats[[name]]), sim_stats[[name]])
      era_iqd = iqd(obs_stats[[name]], era_stats[[name]])
      res[[name]]$iqd = list(c(era = era_iqd, sim = sims_iqd))

      # Compute quantile scores
      quantile_score = function(prob, pred, obs) {
        q = quantile(pred, probs = prob)
        mean(2 * (as.numeric(obs <= q) - prob) * (q - obs))
      }
      era_quantile_score = sapply(
        X = threshold_probs,
        FUN = quantile_score,
        pred = era_stats[[name]],
        obs = obs_stats[[name]]
      )
      sims_quantile_score = sapply(
        X = threshold_probs,
        FUN = quantile_score,
        pred = sim_stats[[name]],
        obs = obs_stats[[name]]
      )
      res[[name]]$quantile_score = list(cbind(era = era_quantile_score, sim = sims_quantile_score))
    }
    res = data.table::rbindlist(res)

    # Save the results
    saveRDS(res, out_path)
  })

# ==============================================================================
# Load and evaluate the results
# ==============================================================================

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

eval

rmse = eval[, .(
  stat,
  sim = sapply(rmse, `[`, "sim"),
  era = sapply(rmse, `[`, "era")
)]
rmse[, .(sim = mean(sim), era = mean(era)), by = "stat"]

mae = eval[, .(
  stat,
  sim = sapply(mae, `[`, "sim"),
  era = sapply(mae, `[`, "era")
)]
mae[, .(sim = mean(sim), era = mean(era)), by = "stat"]

iqd = eval[, .(
  stat,
  sim = sapply(iqd, `[`, "sim"),
  era = sapply(iqd, `[`, "era")
)]
iqd[, .(sim = mean(sim), era = mean(era)), by = "stat"]

eval[, .(
  coverage_90 = mean(coverage_90),
  coverage_95 = mean(coverage_95)
), by = "stat"]

eval[, .(
  e1 = mean(rank_mean - .5),
  e2 = mean(rank_sd - 1 / sqrt(12))
), by = "stat"]

