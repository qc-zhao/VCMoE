\name{plot_coefficients}
\alias{plot_coefficients}
\title{Plot fitted coefficient functions}
\usage{
plot_coefficients(object, type = c("expert", "gating"))
}
\arguments{
\item{object}{A \code{vcmoe} object.}
\item{type}{\code{"expert"} or \code{"gating"}.}
}
\value{
A \code{ggplot} object.
}
\description{
Plots fitted expert or gating coefficient functions over the VCMoE grid.
}
