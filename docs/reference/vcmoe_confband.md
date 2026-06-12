<div id="main" class="col-md-9" role="main">

# Analytic-style confidence bands for a VCMoE fit

<div class="ref-description section level2">

Computes HC0 analytic-style Epanechnikov path confidence bands for VCMoE
fits with `k = 2:10` using the fitted Epanechnikov/scaled
parameterization. High-k intervals are diagnostic-gated and should be
interpreted with the returned block and Hessian diagnostics.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_confband(
  fit,
  data = NULL,
  level = 0.95,
  type = c("pointwise", "simultaneous"),
  coefficient_set = c("expert", "gating", "sigma", "theta"),
  strict = TRUE,
  control = list()
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   fit:

    A `vcmoe` fit.

-   data:

    Optional original data frame. The implementation uses `fit$fitted`;
    refit with `control = list(keep_data = TRUE)` if needed.

-   level:

    Confidence level.

-   type:

    Whether the convenience `lower` and `upper` columns use pointwise
    intervals or simultaneous bands.

-   coefficient\_set:

    Coefficient blocks to return.

-   strict:

    Whether weak local fits should return blocked intervals.

-   control:

    Optional inference controls. HC0 is the only active covariance
    adjustment.

</div>

<div class="section level2">

## Value

A `vcmoe_confband` object with `intervals`, `diagnostics`, and
`settings`.

</div>

<div class="section level2">

## Details

The returned interval table includes pointwise and simultaneous columns,
diagnostic status, block reasons, Hessian condition, effective local
sample size, and SCB metadata. Binomial expert intervals are on the
logit coefficient scale, Negative-Binomial expert intervals are on the
log mean scale, and Negative-Binomial `theta` intervals are nuisance
diagnostics.

</div>

</div>
