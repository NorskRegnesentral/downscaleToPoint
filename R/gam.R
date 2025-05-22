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


#' Create data for plotting all the smooth effects of an mgcv object
#' @export
fast_mgcv_plot_data = function(fit,
                               n = 1000,
                               coord_names = c("lon", "lat"),
                               lon_range = NULL,
                               lat_range = NULL) {
  stopifnot(is(fit, "gam"))
  terms = attr(fit[["terms"]], "term.labels")

  out = list()
  for (smooth in fit[["smooth"]]) {
    terms = smooth[["term"]]
    if (length(terms) == 1) {
      if (!is.null(smooth[["Xu"]])) {
        X_range = range(as.vector(smooth[["Xu"]]) + as.vector(smooth[["shift"]]))
      } else if (!is.null(smooth[["xp"]])) {
        X_range = range(smooth[["xp"]])
      } else {
        stop("I don't know what the X-range for covariate ", terms, " is.")
      }
      data = data.frame(value = seq(min(X_range), max(X_range), length.out = n))
      names(data) = terms
      X = mgcv::PredictMat(smooth, data)
      beta = fit[["coefficients"]][seq(smooth[["first.para"]], smooth[["last.para"]], by = 1)]
      pred = as.vector(X %*% beta)
      cov = smooth[["S"]][[1]] * smooth[["S.scale"]]
      out[[terms]] = data.frame(
        x = data[[terms]],
        y = pred,
        name = terms
      )
    } else {
      if (!all(terms %in% coord_names)) {
        stop(
          "the fit contains a ", length(terms), "-variate model term that\n",
          "I don't know how to deal with. Term names:\n", paste(terms, collapse = ", ")
        )
      }
      lon = if (is.null(lon_range)) {
        seq(0, 360, length.out = n)
      } else {
        seq(min(lon_range), max(lon_range), length.out = sqrt(n))
      }
      lat = if (is.null(lat_range)) {
        seq(-90, 90, length.out = n)
      } else {
        seq(min(lat_range), max(lat_range), length.out = sqrt(n))
      }
      data = as.data.frame(expand.grid(lon = lon, lat = lat))
      X = mgcv::PredictMat(smooth, data)
      beta = fit[["coefficients"]][seq(smooth[["first.para"]], smooth[["last.para"]], by = 1)]
      out[[paste(terms, collapse = ":")]] = data.frame(
        lon = data$lon,
        lat = data$lat,
        value = as.vector(X %*% beta),
        name = paste(terms, collapse = ":")
      )
    }
  }

  out
}
