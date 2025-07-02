
#' Estimate the (possibly weighted) IQD between x and y, using numerical integration
#' @export
iqd = function(x, y, w = function(x) 1, rm_zero = FALSE, zero_threshold = 0) {
  if (rm_zero) {
    x = x[x > zero_threshold]
    y = y[y > zero_threshold]
  }
  F = ecdf(x)
  G = ecdf(y)
  lower = min(c(x, y), na.rm = TRUE)
  upper = max(c(x, y), na.rm = TRUE)
  xx = seq(lower, upper, length.out = 500)
  delta_xx = xx[2] - xx[1]
  w_vals = w(xx)
  f_vals = F(xx)
  g_vals = G(xx)
  fg_vals = f_vals * g_vals * w_vals
  f_vals = f_vals^2 * w_vals
  g_vals = g_vals^2 * w_vals
  int1 = sum((head(f_vals, -1) + tail(f_vals, -1)) * delta_xx / 2)
  int2 = sum((head(g_vals, -1) + tail(g_vals, -1)) * delta_xx / 2)
  int3 = sum((head(fg_vals, -1) + tail(fg_vals, -1)) * delta_xx / 2)
  int1 + int2 - 2 * int3
}

#' @export
bootstrap_skillscores = function(data,
                                 score_name,
                                 data_types,
                                 K_vals,
                                 B = 1000,
                                 row_index = 1,
                                 pairwise_K = TRUE,
                                 pairwise_data_types = TRUE,
                                 probs = c(.025, .975)) {
  # Ensure that we only have one score_name
  stopifnot(length(score_name) == 1)
  # Ensure that we have more than one of either K_vals or data_types
  stopifnot(max(length(data_types), length(K_vals)) > 1)
  # If each value element of `data[[score_name]]` is a matrix instead of a vector,
  # then we need to choose only one of the rows in that matrix for doing our computations
  if (is.matrix(data[[score_name]][[1]])) {
    # Make a copy of data, so we don't overwrite the actual thing outside the function
    data = copy(data)
    data[[score_name]] = lapply(data[[score_name]], function(x) x[row_index, ])
  }
  # Reshape data so it has one column for each value of K
  data = dcast(data[K %in% K_vals], id ~ K, value.var = score_name)

  output = list() # Preallocate the output

  # Loop over different combinations of K_vals and data_types, and compute skill scores
  # with bootstrapped uncertainty intervals.
  # 
  # If we have multiple K_vals, then we need to loop over all pairwise combination of K_vals
  # for each value of data_types. Similarly, if we have multiple data_types, then we need to
  # loop over all pairwise combinations of data_types for each value of K_vals
  # --------------------------------------------------------------------------------------

  # Start by looping over all pairwise combinations of K_vals
  if (pairwise_K) {
    for (data_type in data_types) {
      for (i in seq_along(K_vals)) {
        for (j in seq_len(i - 1)) {
          K1 = K_vals[i]
          K0 = K_vals[j]
          s1_vals = sapply(data[[as.character(K1)]], `[[`, data_type)
          s0_vals = sapply(data[[as.character(K0)]], `[[`, data_type)
          # Compute the actual skill score
          truth = skill_score(s1 = mean(s1_vals), s0 = mean(s0_vals))
          
          # Compute the bootstrapped skill scores
          bootstraps = sapply(
            X = seq_len(B),
            FUN = function(b) {
              boot_index = sample.int(nrow(data), nrow(data), replace = TRUE)
              skill_score(s1 = mean(s1_vals[boot_index]), s0 = mean(s0_vals[boot_index]))
            }
          )
          # Return data.tables with the true skill score, together with quantiles of
          # the bootstrapped skill scores
          output[[length(output) + 1]] = data.table(
            truth = truth,
            lower = quantile(bootstraps, probs[1]),
            upper = quantile(bootstraps, probs[2]),
            K1 = K1,
            K0 = K0,
            data_type = data_type
          )
          # Return a similar data.table that tells us what the numbers would have been
          # if we switched the places of s and s0 in the skill score function
          output[[length(output) + 1]] = data.table(
            truth = opposite_skill_score(truth),
            lower = opposite_skill_score(quantile(bootstraps, probs[2])),
            upper = opposite_skill_score(quantile(bootstraps, probs[1])),
            K1 = K0,
            K0 = K1,
            data_type = data_type
          )
        }
      }
    }
  }

  # Then, loop over all pairwise combinations of data_types
  if (pairwise_data_types) {
    for (K in K_vals) {
      for (i in seq_along(data_types)) {
        for (j in seq_len(i - 1)) {
          s1_vals = sapply(data[[as.character(K)]], `[[`, data_types[i])
          s0_vals = sapply(data[[as.character(K)]], `[[`, data_types[j])
          # Compute the actual skill score
          truth = skill_score(s1 = mean(s1_vals), s0 = mean(s0_vals))
          # Compute the bootstrapped skill scores
          bootstraps = sapply(
            X = seq_len(B),
            FUN = function(b) {
              boot_index = sample.int(nrow(data), nrow(data), replace = TRUE)
              skill_score(s1 = mean(s1_vals[boot_index]), s0 = mean(s0_vals[boot_index]))
            }
          )
          # Return data.tables with the true skill score, together with quantiles of
          # the bootstrapped skill scores
          output[[length(output) + 1]] = data.table(
            truth = truth,
            lower = quantile(bootstraps, probs[1]),
            upper = quantile(bootstraps, probs[2]),
            K = K,
            data_type1 = data_types[i],
            data_type0 = data_types[j]
          )
          # Return a similar data.table that tells us what the numbers would have been
          # if we switched the places of s and s0 in the skill score function
          output[[length(output) + 1]] = data.table(
            truth = opposite_skill_score(truth),
            lower = opposite_skill_score(quantile(bootstraps, probs[2])),
            upper = opposite_skill_score(quantile(bootstraps, probs[1])),
            K = K,
            data_type1 = data_types[j],
            data_type0 = data_types[i]
          )
        }
      }
    }
  }

  output = rbindlist(output, fill = TRUE)
  output$score_name = score_name

  output
}


#' @export
get_scores = function(data,
                      score_name,
                      data_types,
                      K_vals,
                      row_index = 1) {
  # Ensure that we only have one score_name
  stopifnot(length(score_name) == 1)
  # Ensure that we have at least one K_val and one data_type
  stopifnot(min(length(data_types), length(K_vals)) > 0)
  # If each value element of `data[[score_name]]` is a matrix instead of a vector,
  # then we need to choose only one of the rows in that matrix for doing our computations
  if (is.matrix(data[[score_name]][[1]])) {
    # Make a copy of data, so we don't overwrite the actual thing outside the function
    data = copy(data)
    data[[score_name]] = lapply(data[[score_name]], function(x) x[row_index, ])
  }
  # Reshape data so it has one column for each value of K
  data = dcast(data[K %in% K_vals], id ~ K, value.var = score_name)

  output = list() # Preallocate the output

  # Loop over different combinations of K_vals and data_types, and compute scores
  # --------------------------------------------------------------------------------------

  for (data_type in data_types) {
    for (K in K_vals) {
      scores = sapply(data[[as.character(K)]], `[[`, data_type)
      output[[length(output) + 1]] = data.table(
        id = data$id,
        value = scores,
        K = K,
        data_type = data_type,
        score_name = score_name
      )
    }
  }
  output = rbindlist(output, fill = TRUE)

  output
}

#' Function for computing the skill score of s1 and s0
#' @export
skill_score = function(s1, s0) 1 - s1 / s0

#' Function for computing what the skill score of (s0, s1) would have been, when we
#' know what the skill score, s, of (s1, s0) is
opposite_skill_score = function(s) 1 - 1 / (1 - s)

