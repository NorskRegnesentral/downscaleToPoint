% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_tikz}
\alias{plot_tikz}
\title{Turn an R plot into a beautiful pdf made by LaTeX and TikZ,
using the tikzDevice package}
\usage{
plot_tikz(
  plot = NULL,
  expression = NULL,
  file = "Rplots.pdf",
  extra_packages = NULL,
  tex_engine = c("pdflatex", "lualatex"),
  ...
)
}
\description{
Turn an R plot into a beautiful pdf made by LaTeX and TikZ,
using the tikzDevice package
}
