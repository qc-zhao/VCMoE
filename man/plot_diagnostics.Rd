\name{plot_diagnostics}
\alias{plot_diagnostics}
\title{Plot VCMoE fit diagnostics}
\usage{
plot_diagnostics(object)
}
\arguments{
\item{object}{A \code{vcmoe} object.}
}
\value{
A \code{ggplot} object.
}
\description{
Plots convergence, posterior entropy, component proportions, effective local
sample size, and label ambiguity flags over the coefficient grid.
}
\details{
This plot is intended as a first real-data sanity check before interpreting
coefficient functions. Ambiguity or non-convergence at many grid points should
be treated as evidence that the fitted component labels or coefficient paths
need closer review.
}
