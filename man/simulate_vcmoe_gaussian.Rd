\name{simulate_vcmoe_gaussian}
\alias{simulate_vcmoe_gaussian}
\title{Simulate Gaussian VCMoE data}
\usage{
simulate_vcmoe_gaussian(n = 300L, k = 2L, seed = NULL,
  separation = 1, u = NULL, scenario = "well_separated")
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
}
\value{
A list with \code{data} and \code{truth}. The \code{truth} entry includes
component coefficients, gating logits, probabilities, means, standard
deviations, and sampled class labels.
}
\description{
Generates a small Gaussian no-offset VCMoE simulation for tutorials and tests.
}
