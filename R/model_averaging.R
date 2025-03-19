
#' @export
get_bad_donor_index = function(sims, func, n_neighbour) {
  donor_nr = rep(seq_len(n_neighbour), each = ncol(sims) / n_neighbour)
  mean_func_per_donor = sapply(
    X = seq_len(n_neighbour),
    FUN = function(i) mean(apply(sims[, which(donor_nr == i), drop = FALSE], 2, func))
  )
  mean_func_range = 1.5 * iqr(mean_func_per_donor)
  mean_func_lims = quantile(mean_func_per_donor, probs = c(.25, .75)) + mean_func_range * c(-1, 1)
  bad_donors = which(
    mean_func_per_donor < min(mean_func_lims)
    | mean_func_per_donor > max(mean_func_lims)
  )
  bad_donors
}

