\name{vcmoe_scaled_slopes}
\alias{vcmoe_scaled_slopes}
\title{Inspect VCMoE Local-Linear Slopes On The Scaled Basis}
\usage{
vcmoe_scaled_slopes(object, type = c("expert", "gating"), bandwidth = NULL)
}
\arguments{
\item{object}{A \code{vcmoe} object.}
\item{type}{Coefficient block, either \code{"expert"} or \code{"gating"}.}
\item{bandwidth}{Optional bandwidth recorded in the returned attributes.
Defaults to the fitted bandwidth.}
}
\value{
An array with the same dimensions as the stored slope block.
}
\description{
Inspect VCMoE local-linear slopes on the scaled basis.
}
\details{
VCMoE stores slopes on the scaled local-linear basis \code{(u - u0) / h}.
This helper returns the stored scaled-basis slope block.
}
