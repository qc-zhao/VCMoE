\name{vcmoe_fit}
\alias{vcmoe_fit}
\title{Fit a varying-coefficient mixture-of-experts model}
\usage{
vcmoe_fit(formula, data, u, k = 2L, family = "gaussian",
  bandwidth = NULL, u_grid = NULL, control = list(), label = "align",
  parameterization = "a1_epanechnikov_scaled", u_scale = c("unit", "none"))
}
\arguments{
\item{formula}{A formula of the form \code{y ~ expert_terms | gating_terms}.
For grouped Binomial data, use
\code{cbind(success, failure) ~ expert_terms | gating_terms}. For
Negative-Binomial count data, use expert-side
\code{offset(log_size_factor)} for library-size or size-factor offsets.}
\item{data}{A data frame.}
\item{u}{Continuous index column name or numeric vector.}
\item{k}{Number of mixture components. Values 2 through 10 are accepted.
High-k fits are candidate support and require diagnostics.}
\item{family}{Model family. \code{"gaussian"}, \code{"binomial"}, and
\code{"negative-binomial"} are implemented.}
\item{bandwidth}{Kernel bandwidth. If \code{NULL}, a default is used.}
\item{u_grid}{Grid where coefficient functions are estimated.}
\item{control}{Named list overriding EM and label-alignment settings.}
\item{label}{Label strategy. \code{"align"} uses exact global alignment for
\code{k <= 6} and sequential assignment for \code{k >= 7}. \code{"global"}
requests exact global alignment when feasible and falls back to the same
sequential assignment path for \code{k >= 7}. \code{"greedy"} keeps the older
one-step alignment.}
\item{u_scale}{How to transform \code{u} before fitting. The default
\code{"unit"} maps complete-row \code{u} values to \code{[0, 1]};
\code{"none"} leaves \code{u} on the supplied scale. \code{bandwidth} and
\code{u_grid} are interpreted on the transformed analysis scale.}
\item{parameterization}{Estimator convention. The public package uses
\code{"a1_epanechnikov_scaled"}: Epanechnikov density weights
\code{0.75 * (1 - t^2)_+ / h} with \code{t = (u - u0) / h} and stores
local-linear slopes on the scaled \code{(u - u0) / h} basis.}
}
\value{
An object of class \code{vcmoe}.
}
\description{
Fits a Gaussian, Binomial, or Negative-Binomial VCMoE model by local-linear EM
and aligns component labels across the coefficient-function grid.
}
\details{
Rows with missing or non-finite response, covariates, or \code{u} are removed
consistently before fitting, with a warning.

By default, local EM uses Epanechnikov density kernel weights and a scaled
local-linear basis. Complete-row \code{u} values are mapped to \code{[0, 1]}
before fitting unless
\code{u_scale = "none"}. Fitted objects store both the analysis-scale grid and
the original-scale grid metadata.

For \code{family = "gaussian"}, the expert mean and log-standard-deviation
blocks are both local-linear inside each fit. The variance block uses
\code{log(sigma_ic) = delta_c0(u0) + delta_c1(u0) * (u_i - u0) / h};
\code{coef(fit, "sigma")} returns \code{exp(delta_c0(u0))} and
\code{coef(fit, "sigma_slope")} returns \code{delta_c1(u0)}.

For \code{family = "binomial"}, Bernoulli responses must be 0/1 values and
grouped responses must use non-negative finite whole-number success/failure
counts with positive row totals. Binomial expert coefficients are on the logit
success-probability scale. Binomial \code{k >= 3} fits are accepted as
experimental stress support and should not be treated as stable inference-ready
support. Binomial expert logistic M-steps use \code{control$binomial_ridge}
with default value \code{1}, and \code{control$binomial_structured_starts =
TRUE} uses structured local starts to improve single-trial Bernoulli stability.
For single-trial Bernoulli responses only, the gating ridge default is also
strengthened to \code{control$ridge = 1} unless the user explicitly supplies
\code{control$ridge}; grouped Binomial responses keep the global default.

For \code{family = "negative-binomial"}, responses must be finite
non-negative whole-number counts. Expert coefficients are on the log mean count
scale. Use \code{offset(log_size_factor)} in the expert formula to account for
library size or size factors. Gating-side offsets are not supported in v0.
\code{coef(fit, "theta")} returns component-specific NB size parameters.

The default label strategy performs local EM at each grid point, then applies a
post-processing dynamic-programming alignment over the full grid. Transition
costs use the previous local-linear slope to predict the next coefficient value,
plus gating, posterior, and, for Gaussian fits, variance consistency terms.
This improves label path tracking without changing the local EM estimating
equations.

For \code{k <= 6}, \code{label = "align"} uses exact derivative-aware global
alignment. For \code{k >= 7}, \code{label = "align"} and \code{label =
"global"} use sequential pairwise assignment alignment to avoid factorial
permutation growth and record the best-vs-second-best assignment margin for
ambiguity diagnostics.
}
