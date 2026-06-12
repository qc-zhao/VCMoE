\name{vcmoe_parameterization}
\alias{vcmoe_parameterization}
\title{Inspect VCMoE Parameterization Metadata}
\usage{
vcmoe_parameterization(object)
}
\arguments{
\item{object}{A \code{vcmoe} object.}
}
\value{
A named list describing the estimator convention used by the fit, including
kernel weights, local-linear basis scale, gating-logit storage, dispersion
block, label-alignment method, and optimization controls.
}
\description{
Inspect VCMoE parameterization metadata.
}
\details{
The package default is \code{"a1_epanechnikov_scaled"}: Epanechnikov density
weights \code{0.75 * (1 - t^2)_+ / h} with \code{t = (u - u0) / h}, scaled
local-linear slope storage, and centered gating logits. This helper reports
those conventions for reproducible model summaries.
}
