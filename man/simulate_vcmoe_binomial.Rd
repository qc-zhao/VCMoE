\name{simulate_vcmoe_binomial}
\alias{simulate_vcmoe_binomial}
\title{Simulate Binomial VCMoE data}
\usage{
simulate_vcmoe_binomial(n = 300L, k = 2L, seed = NULL,
  separation = 1, u = NULL, scenario = "well_separated", trials = 1L)
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
\item{trials}{Binomial trial counts. Use \code{1} for Bernoulli data, or a
positive integer scalar/vector for grouped Binomial data.}
}
\value{
A list with \code{data} and \code{truth}. Expert truth is on the logit scale.
The \code{truth} entry includes component coefficients, gating logits,
component probabilities, component-specific success probabilities, sampled
class labels, and success/failure counts.
}
\description{
Generates Binomial VCMoE simulations for Bernoulli and grouped-count examples
and tests.
}
