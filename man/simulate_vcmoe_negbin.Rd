\name{simulate_vcmoe_negbin}
\alias{simulate_vcmoe_negbin}
\title{Simulate Negative-Binomial VCMoE count data}
\usage{
simulate_vcmoe_negbin(n = 300L, k = 2L, seed = NULL,
  separation = 1, u = NULL, scenario = "well_separated",
  size_factor = NULL, mean_count = 5)
}
\arguments{
\item{n}{Number of observations.}
\item{k}{Number of components. Values 2 through 10 are supported.}
\item{seed}{Optional random seed.}
\item{separation}{Controls expert separation.}
\item{u}{Optional numeric vector of index values.}
\item{scenario}{Simulation scenario: \code{"well_separated"},
\code{"moderate"}, \code{"near_overlap"}, \code{"crossing"}, or
\code{"imbalanced_gating"}.}
\item{size_factor}{Optional positive size factors. If \code{NULL}, log-normal
size factors are generated.}
\item{mean_count}{Baseline count scale.}
}
\value{
A list with \code{data} and \code{truth}. Expert truth is on the log mean
count scale. The data include \code{size_factor} and \code{log_size_factor}
for use with \code{offset(log_size_factor)}.
}
\description{
Generates Negative-Binomial VCMoE simulations for gene-expression count
examples and tests.
}
