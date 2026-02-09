
#' @export
simulate_tmean_with_donors = function(n_sims,
                                      data,
                                      local_fits,
                                      offset,
                                      time_dep = TRUE) {
  K = nrow(local_fits)
  n_sims_per_local_fit = ceiling(n_sims / K)
  simulations = lapply(
    X = seq_len(K),
    FUN = function(i) {
      if (time_dep) {
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

#' @export
simulate_precip_with_donors = function(n_sims,
                                       data,
                                       local_fits,
                                       offsets,
                                       occurrence_time_dep = TRUE,
                                       intensity_time_dep = TRUE) {
  K = nrow(local_fits)
  n_sims_per_local_fit = ceiling(n_sims / K)
  simulations = lapply(
    X = seq_len(K),
    FUN = function(i) {
      if (occurrence_time_dep) {
        occurrence = simulate_occurrence(
          n = n_sims_per_local_fit,
          data = data,
          init_prob = local_fits$occurrence_prob[i],
          dry_to_wet_fit = local_fits$dry_to_wet[[i]],
          wet_to_wet_fit = local_fits$wet_to_wet[[i]],
          dry_to_wet_offset = offsets$dry_to_wet,
          wet_to_wet_offset = offsets$wet_to_wet
        )
      } else {
        occurrence = simulate_occurrence_notime(
          n = n_sims_per_local_fit,
          data = data,
          fit = local_fits$occurrence[[i]],
          offset = offsets$occurrence
        )
      }
      if (intensity_time_dep) {
        intensity = simulate_intensity(
          n = n_sims_per_local_fit,
          marginal_fit = local_fits$intensity[[i]],
          arma_fit = local_fits$intensity_arma[[i]],
          data = data,
          offset = offsets$intensity
        )
      } else {
        intensity = simulate_intensity_notime(
          n = n_sims_per_local_fit,
          data = data,
          fit = local_fits$intensity[[i]],
          offset = offsets$intensity
        )
      }
      intensity * occurrence
    })
  do.call(cbind, simulations)
}

#' @export
simulate_occurrence_notime = function(n, fit, data, offset = 0) {
  # Compute the linear predictor
  linpred = fast_mgcv_pred(fit, data) + offset
  # Compute the probability of precipitation occurrence
  p = fit$family$linkinv(linpred)
  # Simulate zeros and ones
  res = rbinom(n * length(p), 1, rep(p, n))
  # Return the simulated data in a matrix with `n` columns
  matrix(res, nrow = length(p), ncol = n)
}

#' @export
simulate_occurrence = function(n,
                               dry_to_wet_fit,
                               wet_to_wet_fit,
                               data,
                               init_prob,
                               dry_to_wet_offset = 0,
                               wet_to_wet_offset = 0) {
  # Compute the linear predictors
  dry_to_wet_linpred = fast_mgcv_pred(dry_to_wet_fit, data[-nrow(data)]) + dry_to_wet_offset
  wet_to_wet_linpred = fast_mgcv_pred(wet_to_wet_fit, data[-nrow(data)]) + wet_to_wet_offset
  # Compute the probability of precipitation occurrence
  dry_to_wet_p = dry_to_wet_fit$family$linkinv(dry_to_wet_linpred)
  wet_to_wet_p = wet_to_wet_fit$family$linkinv(wet_to_wet_linpred)
  # Preallocate a matrix of n precipitation occurrence time series
  n_time = length(dry_to_wet_p) + 1
  out = matrix(NA, ncol = n_time, nrow = n)
  # Simulate the n initial dry/wet states
  out[, 1] = rbinom(n, 1, init_prob)
  # I think it is too slow to interatively sample new occurrences after each time step,
  # based on the values of the previous time steps. Therefore, I start by sampling
  # occurrences under the assumption that the last time step always was dry, and under
  # the assumption that the last time step always was wet.
  # Then I just need to loop through each time step and select the simulations
  # that fit with the previous time step
  dry_to_wet = rbinom(n * (n_time - 1), 1, rep(dry_to_wet_p, each = n))
  dim(dry_to_wet) = c(n, n_time - 1)
  wet_to_wet = rbinom(n * (n_time - 1), 1, rep(wet_to_wet_p, each = n))
  dim(wet_to_wet) = c(n, n_time - 1)
  # This is the actual for loop, where we select simulations based on all the previous
  # time steps
  for (i in seq_len(n_time - 1)) {
    out[, i + 1] = out[, i] * wet_to_wet[, i] + (1 - out[, i]) * dry_to_wet[, i]
  }
  # Return the simulated data, with each row representing a time step, and each column
  # representing an ensemble member
  out = t(out)
  out
}

#' @export
simulate_intensity = function(n, marginal_fit, arma_fit, data, offset = 0) {
  # Compute the linear predictor
  linpred = fast_mgcv_pred(marginal_fit, data) + offset
  # Compute the shape and rate of the local GAM
  mu = marginal_fit$family$linkinv(linpred)
  shape = marginal_fit$shape
  rate = marginal_fit$shape / mu
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
  # Censor really large values, to avoind Inf problems when transforming the data
  arma_sims[arma_sims > qnorm(1 - 2 * .Machine$double.eps)] = qnorm(1 - 2 * .Machine$double.eps)
  # Transform the ARMA simulations to have gamma marginal distributions, using the PIT
  res = qgamma(pnorm(arma_sims), shape = shape, rate = rep(rate, n))
  # Return the simulated data in a matrix with `n` columns
  matrix(res, nrow = nrow(arma_sims), ncol = n)
}

#' @export
simulate_intensity_notime = function(n, fit, data, offset = 0) {
  # Compute the linear predictor
  linpred = fast_mgcv_pred(fit, data) + offset
  # Compute the shape and rate of the local GAM
  mu = fit$family$linkinv(linpred)
  shape = fit$shape
  rate = fit$shape / mu
  # Simulate the corresponding precipitation intensities
  res = rgamma(n * length(mu), shape = shape, rate = rate)
  # Return the simulated data in a matrix with `n` columns
  matrix(res, nrow = length(mu), ncol = n)
}


#' @export
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

#' @export
simulate_tmean_notime = function(n, fit, data, offset = 0) {
  # Compute the linear predictor
  linpred = fast_mgcv_pred(fit, data) + offset
  # Simulate the corresponding temperature means
  res = rnorm(n * length(linpred), mean = linpred, sd = fit$sig2)
  # Return the simulated data in a matrix with `n` columns
  matrix(res, nrow = length(linpred), ncol = n)
}
