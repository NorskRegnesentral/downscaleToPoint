#' Turn an R plot into a beautiful pdf made by LaTeX and TikZ,
#' using the tikzDevice package
#' @export
plot_tikz = function(plot = NULL,
                     expression = NULL,
                     file = "Rplots.pdf",
                     extra_packages = NULL,
                     tex_engine = c("pdflatex", "lualatex"),
                     ...) {
  expression = substitute(expression)
  if (is.null(plot) && is.null(expression)) {
    stop("Either `plot` or `expression` must be non-NULL")
  }

  # Create a temporary file for the tikz-output
  tmp = tempfile(tmpdir = getwd())
  # Clean up after yourself on early interrupt
  on.exit(unlink(tmp))

  # I am nor sure about what is going on here...
  opt = options()
  on.exit(options(opt), add = TRUE) #Reset global options on exit
  tikzDevice::setTikzDefaults(overwrite = FALSE)

  # Extract default tex usepackages and extra packages
  tex_packages = c(options()$tikzLatexPackages, "\\usepackage[utf8]{inputenc}")
  extra_packages = c(extra_packages, "bm", "amsmath", "amssymb")
  extra_packages = paste0("\\usepackage{", extra_packages, "}\n")
  tex_packages = unique(c(tex_packages, extra_packages))

  # Open a device for creating a tex-file
  tikzDevice::tikz(tmp, standAlone = TRUE, packages = tex_packages, ...)
  # Call dev.off() on exit in case of interruptions
  current_device = dev.cur()
  on.exit({
    if (current_device %in% dev.list()) dev.off(current_device)
  }, add = TRUE)

  # Plot something into the tex-file
  if (!is.null(plot)) {
    if (any(class(plot) %in% c("gg", "ggplot", "patchwork"))) {
      print(plot)
    } else {
      for (p in plot) print(p)
    }
  } else {
    eval(expression)
  }

  # Finish the creation of the tex-file
  dev.off()

  # Compile to pdf using lualatex
  system2(tex_engine[1], shQuote(tmp))

  # Copy pdf file to final destination
  file.copy(paste0(tmp, ".pdf"), file, overwrite = TRUE)

  # Clean up all temporary files
  tmp_filename = tail(strsplit(tmp, "/")[[1]], 1)
  files_to_clean = grep(tmp_filename, list.files(full.names = TRUE), value = TRUE)
  unlink(files_to_clean)

  # Silently return a TRUE to show that the function completed successfully
  invisible(TRUE)
}

#' @export
pdf_convert = function(in_path,
                       out_paths,
                       format = c("png", "jpeg"),
                       pages = NULL,
                       dpi = 150,
                       ...) {
  suppressWarnings({
    pdftools::pdf_convert(
      pdf = in_path,
      filenames = out_paths,
      format = match.arg(format),
      dpi = dpi,
      ...
    )
  })
}

#' Change the axes of a map plot so it becomes possible to plot using tikzDevice
#' @export
latex_friendly_map_plot = function(x) {
  info = ggplot2::ggplot_build(x)$layout$panel_params[[1]]$graticule
  east_ticks = info$degree[info$type == "E" & info$y_start == 0]
  north_ticks = info$degree[info$type == "N" & info$x_start == 0]
  east_labels = paste0(abs(east_ticks), "$^\\circ$", ifelse(sign(east_ticks) == 1, "E", ifelse(east_ticks == 0, "", "W")))
  north_labels = paste0(abs(north_ticks), "$^\\circ$", ifelse(sign(north_ticks) == 1, "N", ifelse(north_ticks == 0, "", "S")))
  x +
    ggplot2::scale_x_continuous(breaks = east_ticks, labels = east_labels) +
    ggplot2::scale_y_continuous(breaks = north_ticks, labels = north_labels)
}

