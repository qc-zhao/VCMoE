\name{vcmoe_bootstrap}
\alias{vcmoe_bootstrap}
\title{Parametric bootstrap inference for a VCMoE fit}
\usage{
vcmoe_bootstrap(
  fit,
  data,
  u = NULL,
  B = 200L,
  coefficient_set = c("expert", "gating"),
  seed = NULL,
  control = list(),
  min_successful = max(20L, ceiling(0.5 * B)),
  keep_fits = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{fit}{A fitted \code{vcmoe} object. Bootstrap inference supports
\code{k = 2:10}.}

\item{data}{Original data frame used to fit \code{fit}. The function resamples
from \code{data[fit$rows_used, ]}.}

\item{u}{Optional original \code{u} values or column name. If \code{NULL}, the
stored \code{u} column from \code{fit} is reused when available.}

\item{B}{Number of parametric bootstrap replicates.}

\item{coefficient_set}{Coefficient sets to store and summarize:
\code{"expert"}, \code{"gating"}, or both.}

\item{seed}{Optional random seed.}

\item{control}{Named list passed to bootstrap refits. Bandwidth is not
reselected inside bootstrap v0.}

\item{min_successful}{Minimum number of successful replicates expected for
reliable inference. The object is returned when at least two replicates
succeed, but a warning is recorded below this threshold.}

\item{keep_fits}{Whether to store successful bootstrap fit objects.}

\item{verbose}{Whether to print replicate progress messages.}
}
\value{
An object of class \code{vcmoe_bootstrap} with fields \code{fit},
\code{replicates}, \code{replicate_summary}, \code{alignment_summary},
\code{settings}, \code{warnings}, and optionally \code{fits}.
}
\description{
Runs parametric bootstrap inference for a fitted Gaussian, Binomial, or
Negative-Binomial VCMoE model with \code{k = 2:10}.
}
\details{
For Gaussian fits, each bootstrap data set draws a latent component from the
fitted gating probabilities and then draws the response from the selected
component Normal distribution. For Binomial fits, each bootstrap data set draws
success counts from the selected component success probability. For
Negative-Binomial fits, each bootstrap data set draws counts from the selected
component mean and theta. Bernoulli and grouped \code{cbind(success, failure)}
response formats are preserved for Binomial fits.

Each bootstrap replicate is refit with the same formula, family, number of
components, bandwidth, \code{u_grid}, and label strategy as the reference fit.
After the usual within-grid label alignment, one global component permutation
matches the bootstrap coefficient paths back to the reference fit. Ambiguous
bootstrap-to-reference matches are recorded in \code{alignment_summary}.
Exact permutation matching is used for small \code{k}; assignment-based
matching is used when exhaustive permutation is infeasible.

Binomial expert coefficients and intervals are on the logit coefficient scale.
Negative-Binomial expert coefficients and intervals are on the log mean count
scale.
}
