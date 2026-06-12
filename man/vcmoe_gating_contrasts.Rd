\name{vcmoe_gating_contrasts}
\alias{vcmoe_gating_contrasts}
\title{Report Identifiable VCMoE Gating Contrasts}
\usage{
vcmoe_gating_contrasts(object, baseline = NULL, scaled = FALSE)
}
\arguments{
\item{object}{A \code{vcmoe} object.}
\item{baseline}{Component used as the contrast baseline. By default,
\code{k = 2} uses component 2 and \code{k > 2} uses component 1.}
\item{scaled}{If \code{TRUE}, slope contrasts are converted to the scaled
local-linear basis \code{(u - u0) / h}.}
}
\value{
A data frame with one row per grid point, contrast, term, and block.
}
\description{
Report identifiable gating contrasts.
}
\details{
VCMoE stores gating coefficients as centered logits, so the absolute level of
all component logits is not identifiable. Interpretable gating effects are
component contrasts such as component 1 versus component 2 for \code{k = 2},
or component 1 and component 2 versus component 3 for \code{k = 3} comparisons
when \code{baseline = 3}.
}
