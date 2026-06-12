\name{vcmoe_select_bandwidth}
\alias{vcmoe_select_bandwidth}
\title{Select a VCMoE bandwidth by K-fold cross-validation}
\usage{
vcmoe_select_bandwidth(
  formula,
  data,
  u,
  k = 2L,
  family = "gaussian",
  bandwidth_grid = NULL,
  folds = 5L,
  u_grid = NULL,
  control = list(),
  label = "align",
  parameterization = "a1_epanechnikov_scaled",
  u_scale = c("unit", "none"),
  seed = NULL,
  refit = TRUE
)
}
\arguments{
\item{formula}{A formula of the form \code{y ~ expert_terms | gating_terms}.}

\item{data}{A data frame.}

\item{u}{Continuous index column name or numeric vector.}

\item{k}{Number of mixture components. Values from \code{2} through \code{10}
are supported.}

\item{family}{Model family. \code{"gaussian"}, \code{"binomial"}, and
\code{"negative-binomial"} are supported.}

\item{bandwidth_grid}{Candidate bandwidth values. If \code{NULL}, uses
multiples of the default bandwidth.}

\item{folds}{Number of random cross-validation folds.}

\item{u_grid}{Grid where coefficient functions are estimated.}

\item{control}{Named list passed to \code{vcmoe_fit()}.}

\item{label}{Label strategy passed to \code{vcmoe_fit()}.}

\item{u_scale}{\code{u} scaling strategy passed to \code{vcmoe_fit()}. The
default selects bandwidths on the unit-scaled analysis domain.}

\item{parameterization}{Estimator convention passed to \code{vcmoe_fit()} for
each CV fold and the optional final refit.}

\item{seed}{Optional random seed for fold assignment and, when
\code{control$seed} is absent, deterministic CV refits.}

\item{refit}{Whether to refit the final model on all data using the selected
bandwidth.}
}
\value{
An object of class \code{vcmoe_bandwidth_selection} with fields
\code{best_bandwidth}, \code{cv_summary}, \code{cv_details}, \code{cv_folds},
\code{fit}, and \code{settings}.
}
\description{
Selects the kernel bandwidth for a VCMoE model using random K-fold held-out
predictive log-likelihood. The selected bandwidth is the candidate with the
largest held-out likelihood after ranking fully successful candidates ahead of
partial-failure candidates.
}
\details{
The default candidate grid is the current Silverman-style default bandwidth
multiplied by \code{c(0.5, 0.75, 1, 1.25, 1.5, 2)}. Fold assignment is made
only among complete rows that \code{vcmoe_fit()} would keep.

For Gaussian models, validation scoring uses
\code{log sum_c pi_c Normal(y | mu_c, sigma_c)}. For Binomial models, scoring
uses \code{log sum_c pi_c Binomial(success | trials, p_c)}. Binomial Bernoulli
and grouped \code{cbind(success, failure)} responses are supported.
For Negative-Binomial models, scoring uses
\code{log sum_c pi_c NB(y | mu_c, theta_c)}.

Bandwidth selection supports \code{k = 2:10} for Gaussian, Binomial, and
Negative-Binomial models. High-k candidates use the same held-out predictive
likelihood scoring and should be interpreted together with the returned fit
diagnostics.
The selected object records the fitting parameterization and \code{u} scaling
strategy in \code{settings$parameterization} and \code{settings$u_scale}.
}
