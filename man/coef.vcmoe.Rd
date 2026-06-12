\name{coef.vcmoe}
\alias{coef.vcmoe}
\title{Extract VCMoE coefficients}
\usage{
\method{coef}{vcmoe}(object, type = c("all", "expert", "gating", "sigma", "sigma_slope", "theta"), ...)
}
\arguments{
\item{object}{A \code{vcmoe} object.}
\item{type}{Coefficient block to return.}
\item{...}{Unused.}
}
\value{
A list or array of fitted coefficient functions.
}
\description{
Extracts expert coefficients, gating coefficients, Gaussian variance
intercepts, Gaussian log-sigma local-linear slopes, Negative-Binomial theta,
or all fitted coefficient blocks from a VCMoE fit.
}
\details{
For Gaussian fits, \code{coef(fit, "sigma")} returns the component-specific
standard deviation function at each \code{u_grid} point, and
\code{coef(fit, "sigma_slope")} returns the scaled local-linear slope of
\code{log(sigma)} on the \code{(u - u0) / h} basis.

For Binomial fits, expert coefficients are on the logit success-probability
scale and \code{coef(fit, "sigma")} returns \code{NULL}. For
Negative-Binomial fits, expert coefficients are on the log mean count scale,
\code{coef(fit, "theta")} returns the component-specific size parameter, and
\code{coef(fit, "sigma")} returns \code{NULL}.
}
