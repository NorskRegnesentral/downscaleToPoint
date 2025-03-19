#' Create a progress bar for tracking long-running processes.
#' This is a thin wrapper around the progress package
#' @export
progress_bar = function(n) {
  pb = progress::progress_bar$new(
    format = ":percent [:bar] time elapsed: :elapsedfull, eta: :eta",
    total = n, width = 70, clear = FALSE)
  # pb$tick() throws an error if we tick too many times, which can potentially stop
  # a script that is otherwise working fine. We don't want that to happen,
  # so we change the tick function slightly
  res = list(
    terminate = pb$terminate,
    tick = function(...) tryCatch(pb$tick(...), error = function(e) NULL))
  res
}

