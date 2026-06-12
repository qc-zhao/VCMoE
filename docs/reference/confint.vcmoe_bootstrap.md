<div id="main" class="col-md-9" role="main">

# Bootstrap confidence intervals for VCMoE coefficients

<div class="ref-description section level2">

Summarizes pointwise or simultaneous bootstrap intervals for expert or
gating coefficient functions.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
# S3 method for class 'vcmoe_bootstrap'
confint(
  object,
  parm = c("expert", "gating"),
  level = 0.95,
  type = c("pointwise", "simultaneous"),
  ...
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   object:

    A `vcmoe_bootstrap` object.

-   parm:

    Coefficient set to summarize: `"expert"`, `"gating"`, or both.

-   level:

    Confidence level.

-   type:

    Interval type. `"pointwise"` uses percentile bootstrap intervals.
    `"simultaneous"` uses a max standardized bootstrap-deviation band
    over `u` for each component and term.

-   ...:

    Unused.

</div>

<div class="section level2">

## Value

A tidy data frame with columns `coefficient_set`, `term`, `component`,
`u`, `estimate`, `se`, `lower`, `upper`, `type`, `level`, and
`n_successful`.

</div>

<div class="section level2">

## Details

Pointwise intervals are percentile intervals at each grid point.
Simultaneous bands compute bootstrap standard errors and use the
empirical quantile of the maximum standardized absolute deviation over
the `u_grid`. Near-zero standard errors are floored internally to avoid
division by zero.

</div>

</div>
