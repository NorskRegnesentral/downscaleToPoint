#' Perform fast prediction for gam objects with smooths and an intercept
#' @export
fast_mgcv_pred = function(fit, data, detailed = FALSE) {
  stopifnot(is(fit, "gam"))
  if (!is.data.table(data)) data = as.data.table(data)
  terms = attr(fit[["terms"]], "term.labels")
  stopifnot(all(terms %in% names(data)))

  preds = list()

  if (names(fit[["coefficients"]])[1] == "(Intercept)") {
    preds$intercept = fit[["coefficients"]][1]
  }

  # Add smooth terms
  for (smooth in fit[["smooth"]]) {
    terms = smooth[["term"]]
    compressed_data = unique(data[, ..terms])
    X = mgcv::PredictMat(smooth, compressed_data)
    beta = fit[["coefficients"]][seq(smooth[["first.para"]], smooth[["last.para"]], by = 1)]
    compressed_data$pred = as.vector(X %*% beta)
    pred = merge(data[, ..terms], compressed_data, by = terms, sort = FALSE)[["pred"]]
    preds[[paste(terms, collapse = ":")]] = pred
  }

  # Add fixed terms
  fixed_terms = attr(fit[["pterms"]], "term.labels")
  for (term in fixed_terms) {
    if (is.numeric(data[[term]])) {
      coeff_index = which(names(fit[["coefficients"]]) == term)
      if (length(coeff_index) != 1) stop("Could not find correct fixed coefficient")
      pred = data[[term]] * fit[["coefficients"]][coeff_index]
    } else if (is.logical(data[[term]])) {
      coeff_index = which(names(fit[["coefficients"]]) == paste0(term, "TRUE"))
      if (length(coeff_index) != 1) stop("Could not find correct fixed coefficient")
      pred = data[[term]] * fit[["coefficients"]][coeff_index]
    } else {
      stop("Covariate data must be numeric or logical")
    }
    preds[[term]] = pred
  }

  preds$pred = Reduce(`+`, preds)

  if (detailed) {
    out = preds
  } else {
    out = preds$pred
  }

  out
}
