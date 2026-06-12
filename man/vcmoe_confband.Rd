\name{vcmoe_confband}
\alias{vcmoe_confband}
\title{Analytic-style confidence bands for a VCMoE fit}
\usage{
vcmoe_confband(
  fit,
  data = NULL,
  level = 0.95,
  type = c("pointwise", "simultaneous"),
  coefficient_set = c("expert", "gating", "sigma", "theta"),
  strict = TRUE,
  control = list()
)
}
\arguments{
\item{fit}{A \code{vcmoe} fit.}
\item{data}{Optional original data frame. The implementation uses
\code{fit$fitted}; refit with \code{control = list(keep_data = TRUE)} if
needed.}
\item{level}{Confidence level.}
\item{type}{Whether the convenience \code{lower} and \code{upper} columns use
pointwise intervals or simultaneous bands.}
\item{coefficient_set}{Coefficient blocks to return.}
\item{strict}{Whether weak local fits should return blocked intervals.}
\item{control}{Optional inference controls. HC0 is the only active covariance
adjustment.}
}
\value{
A \code{vcmoe_confband} object with \code{intervals}, \code{diagnostics}, and
\code{settings}.
}
\description{
Computes HC0 analytic-style Epanechnikov path confidence bands for
VCMoE fits with \code{k = 2:10} using the fitted Epanechnikov/scaled
parameterization. High-k intervals are diagnostic-gated and should be
interpreted with the returned block and Hessian diagnostics.
}
\details{
The returned interval table includes pointwise and simultaneous columns,
diagnostic status, block reasons, Hessian condition, effective local sample
size, and SCB metadata. Binomial expert intervals are on the logit coefficient
scale, Negative-Binomial expert intervals are on the log mean scale, and
Negative-Binomial \code{theta} intervals are nuisance diagnostics.
}
