% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tailor_methods.R
\name{get_1D_mixtures}
\alias{get_1D_mixtures}
\title{get_1D_mixtures}
\usage{
get_1D_mixtures(
  data,
  params,
  max_mixture = 3,
  use_ICL = FALSE,
  sample_fraction = 0.05,
  parallel = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{data}{A flowSet, flowFrame or a matrix containing
events along the rows, markers along columns.}

\item{params}{A list of markers to use; must be subset
of colnames(data).}

\item{max_mixture}{Will attempt to model each marker as
k mixture components, for
1 <= k <= max_mixture. The best k is chosen based on a
biased version of the Integrated Complete Likelihood (ICL).}

\item{sample_fraction}{A number between 0 and 1: the
fraction of data points used in the
calculation of 1D mixture components, to improve runtime.}

\item{parallel}{Boolean flag; if true, uses multithreading
to process markers in parallel.}

\item{verbose}{Boolean flag; if true, outputs timing and
milestone information.}

\item{bias_ICL}{Bias the ICL towards more mixture components.}
}
\value{
A named list of 1D mixture models, giving mixture
proportions, means and variances for each marker.
}
\description{
Computes 1D mixture model for each marker
separately, for use in binning step of tailor. It is
difficult to find settings which work for all datasets.
Therefore, it is recommended to inspect the results with
inspect_1D_mixtures, and run get_1D_mixtures_custom
for problematic markers.
}
\examples{
fileName <- system.file("extdata", "sampled_flowset_old.rda",
                        package = "Tailor")
load(fileName)
tailor_params <- flowCore::colnames(fs_old)[c(7:9, 11:22)]
mixtures_1D <- get_1D_mixtures(fs_old, tailor_params)
}
