\name{vcmoe_scaled_slopes}
\alias{vcmoe_scaled_slopes}
\title{Convert VCMoE Slopes To A Scaled Local-Linear Basis}
\usage{
vcmoe_scaled_slopes(object, type = c("expert", "gating"), bandwidth = NULL)
}
\arguments{
\item{object}{A \code{vcmoe} object.}
\item{type}{Coefficient block, either \code{"expert"} or \code{"gating"}.}
\item{bandwidth}{Optional bandwidth used for conversion. Defaults to the fitted bandwidth.}
}
\value{
An array with the same dimensions as the stored slope block.
}
\description{
Convert VCMoE local-linear slopes to scaled-basis slopes.
}
\details{
VCMoE stores slopes directly for the scaled basis
\code{(u - u0) / h}, so this helper returns stored slopes unchanged. Legacy
Gaussian/raw fits store slopes for \code{u - u0}; for those fits the equivalent
scaled-basis slope is \code{h * slope_raw}. This helper performs the needed
conversion without changing the fitted object.
}
