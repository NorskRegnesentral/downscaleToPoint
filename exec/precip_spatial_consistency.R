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
model_dir = file.path(data_dir, "models", "precipitation")
image_dir = file.path(data_dir, "images", "precipitation")
local_fits_dir = file.path(model_dir, "local_fits")
out_dir = file.path(data_dir, "spatial_consistency", "precipitation")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

meta_path = file.path(data_dir, "meta.rds")
global_fit_path = file.path(model_dir, "global.rds")

# Define other useful variables
# ------------------------------------------------------------------------------
n_cores = 8 # Number of cores to use for running code in parallel
n_sims = 150 # Number of ensembles to simulate during the downscaling
diff_lengths = c(1, 3, 7) # Which n-day-differences to evaluate in the cross-validation
B = 1000 # Number of bootstraps to use when bootstrapping

K = 15

n_neighbour_min = 8
neighbour_radius = 100e3

overwrite = FALSE

# Thresholds for computing threshold weighted IQD scores during the cross-validation
threshold_probs = c(.01, .05, .1, .9, .95, .99)

# Random seeds for reproducibility
set.seed(20260130)
base_seed = sample.int(1e8, 1)
seed_jump = sample.int(1e4, 1)

score_info = as.data.frame(t(do.call(cbind, list(
  c("rmse", "RMSE", 1),
  c("mae", "MAE", 1),
  c("iqd", "IQD", 1),
  c("quantile_score", "Q95", 5),
  c("quantile_score", "Q99", 6)
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
# Evaluation
# ==============================================================================

global_fit = readRDS(global_fit_path)

start_time = Sys.time()
success = parallel::mclapply(
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

    # Locate all stations that are closer than neighbour_radius
    neighbouring_ids = station_meta$id[dists <= neighbour_radius]
    n_neighbouring_ids = length(neighbouring_ids)

    # Load data for all these stations
    obs = load_station_data(
      meta = station_meta[id %in% neighbouring_ids],
      data_dir = data_dir,
      rm_na = FALSE
    )

    # Add the necessary covariates for performing downscaling
    obs[, let(
      yday = yday(date),
      year = year(date),
      month = month(date),
      station_elevation = log(station_elevation + 1),
      era_log_precip = log(era_precip + 1),
      era_precip_bool = era_precip > 0,
      precip_bool = precip > 0,
      day_count = as.integer(date) - as.integer(min(date)) + 1L
    )]
    obs[, let(
      next_era_precip_bool = c(tail(era_precip_bool, -1), NA),
      next_era_log_precip = c(tail(era_log_precip, -1), NA),
      next_era_tmean = c(tail(era_tmean, -1), NA)
    ), by = "id"]
    obs[, let(
      era_log_precip_change = next_era_log_precip - era_log_precip,
      era_tmean_change = next_era_tmean - era_tmean
    )]

    # Only keep data from dates where at least `n_neighbour_min` of the stations
    # actually have any data
    obs[, let(n_obs_per_date = sum(!is.na(precip))), by = c("date")]
    if (obs[n_obs_per_date >= n_neighbour_min, .N] == 0) return(FALSE)
    date_range = obs[n_obs_per_date >= n_neighbour_min, range(date)]
    obs = obs[date >= min(date_range)][date <= max(date_range)]

    n_obs_per_id = table(obs$id)
    if (any(n_obs_per_id < 50)) obs = obs[id %in% names(n_obs_per_id[n_obs_per_id >= 50])]

    obs_ids = unique(obs$id)
    obs_dates = sort(unique(obs$date))

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

      local_obs = obs[id == obs_ids[j]][order(date)]

      # Compute offsets from the global GAMs
      offsets = list()
      for (name in names(global_fit)) {
        if (grepl("wet", name)) {
          offsets[[name]] = fast_mgcv_pred(global_fit[[name]], local_obs[-nrow(local_obs)])
        } else {
          offsets[[name]] = fast_mgcv_pred(global_fit[[name]], local_obs)
        }
      }

      sim = simulate_precip_with_donors(
        n_sims = n_sims,
        data = local_obs,
        local_fits = local_fits,
        offsets = offsets,
        intensity_time_dep = TRUE,
        occurrence_time_dep = TRUE
      )

      # Check if any of the donor stations appear to be outliers and
      # Remove them if this is the case
      bad_donors = get_bad_donor_index(sim, mean, K)
      if (length(bad_donors) > 0) {
        sim = simulate_precip_with_donors(
          n_sims = n_sims,
          data = local_obs,
          local_fits = local_fits[-bad_donors, ],
          offsets = offsets,
          intensity_time_dep = TRUE,
          occurrence_time_dep = TRUE
        )
      }

      # Change sim to NA for all dates where local_obs$precip was NA
      sim[is.na(local_obs$precip), ] = NA

      # Pad the output with NAs, so that all simulations from all stations have the
      # same number of rows
      sims[[obs_ids[j]]] = matrix(nrow = length(obs_dates), ncol = n_sims)
      sims[[obs_ids[j]]][obs_dates %in% local_obs$date, ] = sim[, seq_len(n_sims)]
    }
    sims = do.call(abind::abind, list(sims, along = 3))

    obs_stats = suppressWarnings({
      obs[, .(
        mean = mean(precip, na.rm = TRUE),
        median = median(precip, na.rm = TRUE),
        sd = sd(precip, na.rm = TRUE),
        min = min(precip, na.rm = TRUE),
        max = max(precip, na.rm = TRUE),
        n_stations = sum(!is.na(precip))
      ), by = "date"][order(date)]
    })

    era_stats = obs[, .(
      mean = mean(era_precip, na.rm = TRUE),
      median = median(era_precip, na.rm = TRUE),
      sd = sd(era_precip, na.rm = TRUE),
      min = min(era_precip, na.rm = TRUE),
      max = max(era_precip, na.rm = TRUE),
      n_stations = sum(!is.na(precip))
    ), by = "date"][order(date)]

    # Remove rows corresponding to few observations
    # We had to do this after the simulations, because of the Markov process
    # used during the simulations, which doesn't handle NAs or missing dates
    sims = sims[which(obs_stats$n_stations >= n_neighbour_min), , ]
    era_stats = era_stats[n_stations >= n_neighbour_min]
    obs_stats = obs_stats[n_stations >= n_neighbour_min]

    # Compute different stats over the entire domain
    sim_stats = list(
      mean = apply(sims, 2, matrixStats::rowMeans2, na.rm = TRUE),
      median = apply(sims, 2, matrixStats::rowMedians, na.rm = TRUE),
      sd = apply(sims, 2, matrixStats::rowSds, na.rm = TRUE),
      min = apply(sims, 2, matrixStats::rowMins, na.rm = TRUE),
      max = apply(sims, 2, matrixStats::rowMaxs, na.rm = TRUE)
    )

    # Evaluate the performance of the spatial stat ensembles
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
        n_neighbour_max = max(obs_stats$n_stations),
        n_neighbour_mean = mean(obs_stats$n_stations),
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
    message("Done with i = ", i)

    TRUE
  })

unlist(success) |> table()

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
eval$K = K # This is stupid, but necessary for bootstrap_skillscores()

stat_names = unique(eval$stat)

data_types = c("era", "sim")
bootstrap_data = list()
for (i in seq_len(nrow(score_info))) {
  set.seed(base_seed + i * seed_jump)
  for (stat_name in stat_names) {
    tmp = bootstrap_skillscores(
      data = eval[stat == stat_name],
      score_name = score_info$name[i],
      data_types = data_types,
      K_vals = K,
      row_index = score_info$row_index[i]
    )
    tmp$score_name = score_info$shortname[i]
    tmp$stat = stat_name
    bootstrap_data[[length(bootstrap_data) + 1]] = tmp
  }
}
bootstrap_data = rbindlist(bootstrap_data)
bootstrap_data[, let(
  score_name = factor(score_name, levels = score_info$shortname),
  stat = factor(
    stat,
    levels = c("mean", "sd", "median", "min", "max"),
    labels = c("Mean", "SD", "Median", "Min", "Max")
  )
)]

plot = bootstrap_data[data_type1 == "sim"] |>
  copy() |>
  _[, let(
    truth = pmax(truth, -1),
    lower = pmax(lower, -1),
    upper = pmax(upper, -1)
  )] |>
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_point(
    aes(x = score_name, y = truth, col = stat, group = stat),
    position = position_dodge(.2),
    size = rel(.8)
  ) +
  geom_errorbar(
    aes(x = score_name, ymin = lower, ymax = upper, col = stat, group = stat),
    position = position_dodge(.2)
  ) +
  scale_y_continuous(breaks = seq(-10, 1, by = .2), limits = c(-1, 1)) +
  theme_light() +
  theme(
    strip.text = element_text(colour = "black", size = rel(1)),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
    axis.text.x = element_text(size = rel(1.1), angle = 70, vjust = .5),
    text = element_text(size = 15)
  ) +
  theme(legend.position = "top") +
  labs(x = "Scoring function", y = "$\\tilde S_{\\text{skill}}(S_1, S_0)$", col = "Statistic")

plot_tikz(
  file = file.path(image_dir, "spatial_consistency_scores.pdf"),
  plot = plot,
  width = 8,
  height = 5
)

"Maybe add a title, so we can combine this with the temperature scores in a nice way?"
