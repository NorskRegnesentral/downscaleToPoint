#' Load weather station data, ERA data and covariate data from a set of locations
#'
#' @description
#' Precipitation flag definition
#' (from https://www.ncei.noaa.gov/data/global-summary-of-the-day/doc/readme.txt)
#' ------------------------------------------------------------------------------
#' A = 1 report of 6-hour precipitation amount.
#' B = Summation of 2 reports of 6-hour precipitation amount.
#' C = Summation of 3 reports of 6-hour precipitation amount.
#' D = Summation of 4 reports of 6-hour precipitation amount.
#' E = 1 report of 12-hour precipitation amount.
#' F = Summation of 2 reports of 12-hour precipitation amount.
#' G = 1 report of 24-hour precipitation amount.
#' H = Station reported '0' as the amount for the day (eg, from 6-hour reports),
#'   but also reported at least one occurrence of precipitation in hourly observations. 
#'   This could indicate a trace occurred, but should be considered as incomplete 
#'   data for the day.
#' I = Station did not report any precipitation data for the day and did not report any
#'   occurrences of precipitation in its hourly observations. It's still possible that
#'   precipitation occurred but was not reported.
#' @export
load_station_data = function(meta,
                             data_dir,
                             verbose = FALSE,
                             rm_bad_flags = FALSE,
                             bad_flag_start = "H",
                             rm_na = TRUE,
                             era_stats = FALSE,
                             rm_leading_na = TRUE) {
  stopifnot(is.data.frame(meta))
  data = vector("list", nrow(meta))
  pb = progress_bar(nrow(meta))
  for (i in seq_len(nrow(meta))) {
    gsod_data = readRDS(file.path(data_dir, "gsod", paste0(meta$id[i], ".rds")))
    era_data = readRDS(file.path(data_dir, "era", paste0(meta$id[i], ".rds")))
    era_data = rbindlist(era_data)
    setnames(era_data, c("total_precipitation", "2m_temperature"), c("era_precip", "era_tmean"))
    gsod_data[, let(tmin = NULL, tmax = NULL)]
    # Change the precip units to mm/day and the tmean units to Celsius
    era_data[, let(
      era_tmean = era_tmean - 273.15, # Kelvin to Celsius
      era_precip = era_precip * 1000 # m/day to mm/day
    )]
    setkey(era_data, date)
    if (era_stats) {
      era_data[, let(
        era_tmean_mean = mean(era_tmean),
        era_tmean_sd = sd(era_tmean),
        era_prec_intensity_mean = mean(era_precip[era_precip > 0]),
        era_prec_intensity_sd = sd(era_precip[era_precip > 0])
      )]
    }
    # Possibly remove observations with bad precip flags
    if (rm_bad_flags) gsod_data[precip_flag >= bad_flag_start, let(precip = NA)]
    # Merge the GSOD data and ERA data
    setkey(gsod_data, date)
    data[[i]] = merge(gsod_data, era_data, by = "date", all.y = TRUE)
    if (rm_leading_na) {
      obs_range = range(which(!is.na(data[[i]]$precip)))
      data[[i]] = data[[i]][obs_range[1]:obs_range[2], ]
    }
    if (rm_na) data[[i]] = data[[i]][!is.na(precip), ]
    # Add covariates
    data[[i]][, let(
      lon = meta$lon[i],
      lat = meta$lat[i],
      station_elevation = meta$elev[i],
      grid_elevation_mean = meta$elev_mean[i],
      grid_elevation_sd = meta$elev_sd[i],
      elevation_diff = meta$elev[i] - meta$elev_mean[i]
    )]
    # Remove NA from the id column
    data[[i]]$id = meta$id[i]
    if (verbose) pb$tick()
  }
  pb$terminate()
  rbindlist(data)
}
