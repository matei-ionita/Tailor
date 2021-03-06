% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tailor_methods.R
\name{plot_kdes_global}
\alias{plot_kdes_global}
\title{plot_kdes_global}
\usage{
plot_kdes_global(data, data_overlay = NULL, params = NULL)
}
\arguments{
\item{data}{A flowSet, flowFrame or a matrix containing events
along the rows, markers along columns.}

\item{data_overlay}{Some different data for comparison.}

\item{params}{A list of markers to use; must be subset of
colnames(data).}
}
\value{
Plots of kernel density estimate for each chosen parameter.
}
\description{
Plot 1D kdes of a dataset for visual inspection.
}
\examples{
fileName <- system.file("extdata", "sampled_flowset_old.rda",
                         package = "Tailor")
load(fileName)
tailor_params <- flowCore::colnames(fs_old)[c(7:9, 11:22)]
plot_kdes_global(fs_old, params = tailor_params)
}
