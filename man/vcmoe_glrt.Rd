\name{vcmoe_glrt}
\alias{vcmoe_glrt}
\title{Coefficient-specific GLRT for VCMoE coefficient variation}
\usage{
vcmoe_glrt(
  fit,
  data,
  test = c("coefficient", "constant_all"),
  coefficient_set = c("expert", "gating", "sigma", "theta"),
  component = NULL,
  term = NULL,
  calibration = c("analytic_epanechnikov", "bootstrap", "both", "none",
    "parametric_bootstrap"),
  B = 200L,
  seed = NULL,
  control = list(),
  refit_control = list(),
  verbose = FALSE
)
}
\arguments{
\item{fit}{A \code{vcmoe} fit.}
\item{data}{Original data frame used to fit \code{fit}.}
\item{test}{Test type. \code{"coefficient"} tests one coefficient function;
\code{"constant_all"} tests all fitted coefficient functions jointly.}
\item{coefficient_set}{Coefficient block for coefficient-specific tests.}
\item{component}{Component label or index for coefficient-specific tests.}
\item{term}{Term name for coefficient-specific tests.}
\item{calibration}{Calibration method. \code{"analytic_epanechnikov"} uses
the Epanechnikov modified chi-square calibration; \code{"bootstrap"}
uses parametric bootstrap calibration; \code{"both"} reports both.
\code{"parametric_bootstrap"} is accepted as a backwards-compatible alias for
\code{"bootstrap"}.}
\item{B}{Number of bootstrap calibration replicates.}
\item{seed}{Optional random seed.}
\item{control}{Controls for the constrained null optimizer.}
\item{refit_control}{Controls overriding bootstrap full-model refits.}
\item{verbose}{Whether to message progress.}
}
\value{
A \code{vcmoe_glrt} object with the observed \code{lambda}, analytic
statistic, null fit, optional bootstrap replicate summary, and calibrated
p-value when available.
}
\description{
Fits a constrained null under the same local objective and computes a
generalized likelihood-ratio statistic against the supplied VCMoE fit.
}
\details{
The primary path is coefficient-specific. The selected coefficient function is
constrained to be constant in \code{u}: its local-linear slope is fixed to zero
and its intercept is shared across grid points. The null is re-optimized rather
than obtained by post-hoc averaging.

Analytic Epanechnikov calibration requires Epanechnikov density weights,
scaled local-linear basis, and unit-scaled
\code{u}. It reports \code{lambda = ell_full - ell_null},
\code{analytic_statistic = rK * lambda}, and a modified chi-square
\code{analytic_p_value}. The diagnostic \code{lrt_statistic = 2 * lambda} is
reported separately and is not used for the analytic p-value. Ridge
penalties may be used to stabilize fitting but are excluded from the GLRT
likelihood ratio.

\code{vcmoe_glrt()} supports fitted \code{k = 2:10} models for both
coefficient-specific and \code{"constant_all"} nulls. For \code{k > 2}, gating
coefficient tests use identifiable baseline contrasts, e.g.
\code{component3_vs_component1}. Reported p-values should be interpreted
together with fit convergence, label ambiguity, component proportion, and
null-optimizer diagnostics.
}
