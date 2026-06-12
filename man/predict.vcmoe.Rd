\name{predict.vcmoe}
\alias{predict.vcmoe}
\title{Predict from a VCMoE fit}
\usage{
\method{predict}{vcmoe}(object, newdata = NULL, u = NULL,
  type = c("mean", "posterior", "component", "prior"), ...)
}
\arguments{
\item{object}{A \code{vcmoe} object.}
\item{newdata}{Optional data frame.}
\item{u}{Optional index values for \code{newdata}.}
\item{type}{Prediction type.}
\item{...}{Unused.}
}
\value{
A vector or matrix depending on \code{type}.
}
\description{
Returns fitted means, component-specific means, posterior probabilities, or
gating probabilities.
}
\details{
For Gaussian fits, \code{type = "component"} returns component-specific means
and \code{type = "mean"} returns the posterior-weighted fitted mean. For
Binomial fits, \code{type = "component"} returns component-specific success
probabilities and \code{type = "mean"} returns the marginal success
probability. For Negative-Binomial fits, \code{type = "component"} returns
component-specific mean counts and \code{type = "mean"} returns the marginal
mean count.
}
