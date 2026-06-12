<div id="main" class="col-md-9" role="main">

# Extract VCMoE coefficients

<div class="ref-description section level2">

Extracts expert coefficients, gating coefficients, Gaussian variance
intercepts, Gaussian log-sigma local-linear slopes, Negative-Binomial
theta, or all fitted coefficient blocks from a VCMoE fit.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
# S3 method for class 'vcmoe'
coef(object, type = c("all", "expert", "gating", "sigma", "sigma_slope", "theta"), ...)
```

</div>

</div>

<div class="section level2">

## Arguments

-   object:

    A `vcmoe` object.

-   type:

    Coefficient block to return.

-   ...:

    Unused.

</div>

<div class="section level2">

## Value

A list or array of fitted coefficient functions.

</div>

<div class="section level2">

## Details

For Gaussian fits, `coef(fit, "sigma")` returns the component-specific
standard deviation function at each `u_grid` point, and
`coef(fit, "sigma_slope")` returns the scaled local-linear slope of
`log(sigma)` on the `(u - u0) / h` basis.

For Binomial fits, expert coefficients are on the logit
success-probability scale and `coef(fit, "sigma")` returns `NULL`. For
Negative-Binomial fits, expert coefficients are on the log mean count
scale, `coef(fit, "theta")` returns the component-specific size
parameter, and `coef(fit, "sigma")` returns `NULL`.

</div>

</div>
