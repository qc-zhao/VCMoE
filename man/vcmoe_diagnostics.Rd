\name{vcmoe_diagnostics}
\alias{vcmoe_diagnostics}
\title{Summarize VCMoE fit diagnostics}
\usage{
vcmoe_diagnostics(object)
}
\arguments{
\item{object}{A \code{vcmoe} object.}
}
\value{
A data frame with one row per coefficient grid point.
}
\description{
Returns a compact diagnostic table for reviewing whether a fitted VCMoE model
is reliable enough to interpret.
}
\details{
The table includes convergence status, iterations, local log-likelihood,
local-weighted posterior entropy, label ambiguity flags, alignment margin,
effective local sample size, local-weighted component posterior proportions,
and Binomial expert optimizer diagnostics when available.

Posterior entropy and component proportions use the same local kernel weights
as the fitted grid point when the fit retains training data. If the fit was
created with \code{control$keep_data = FALSE}, component proportions fall back
to unweighted posterior means and effective local sample size is \code{NA}.
}
