\name{confint.vcmoe_bootstrap}
\alias{confint.vcmoe_bootstrap}
\title{Bootstrap confidence intervals for VCMoE coefficients}
\usage{
\method{confint}{vcmoe_bootstrap}(
  object,
  parm = c("expert", "gating"),
  level = 0.95,
  type = c("pointwise", "simultaneous"),
  ...
)
}
\arguments{
\item{object}{A \code{vcmoe_bootstrap} object.}

\item{parm}{Coefficient set to summarize: \code{"expert"}, \code{"gating"},
or both.}

\item{level}{Confidence level.}

\item{type}{Interval type. \code{"pointwise"} uses percentile bootstrap
intervals. \code{"simultaneous"} uses a max standardized bootstrap-deviation
band over \code{u} for each component and term.}

\item{...}{Unused.}
}
\value{
A tidy data frame with columns \code{coefficient_set}, \code{term},
\code{component}, \code{u}, \code{estimate}, \code{se}, \code{lower},
\code{upper}, \code{type}, \code{level}, and \code{n_successful}.
}
\description{
Summarizes pointwise or simultaneous bootstrap intervals for expert or gating
coefficient functions.
}
\details{
Pointwise intervals are percentile intervals at each grid point. Simultaneous
bands compute bootstrap standard errors and use the empirical quantile of the
maximum standardized absolute deviation over the \code{u_grid}. Near-zero
standard errors are floored internally to avoid division by zero.
}
