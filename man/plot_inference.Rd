\name{plot_inference}
\alias{plot_inference}
\title{Plot bootstrap inference intervals}
\usage{
plot_inference(
  object,
  coefficient_set = "expert",
  type = c("pointwise", "simultaneous"),
  level = 0.95
)
}
\arguments{
\item{object}{A \code{vcmoe_bootstrap} object.}

\item{coefficient_set}{Coefficient set to plot: \code{"expert"} or
\code{"gating"}.}

\item{type}{Interval type passed to \code{confint()}.}

\item{level}{Confidence level.}
}
\value{
A \code{ggplot} object.
}
\description{
Plots fitted coefficient functions with bootstrap pointwise intervals or
simultaneous bands.
}
