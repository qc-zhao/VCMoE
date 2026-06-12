<div id="main" class="col-md-9" role="main">

# Plot bootstrap inference intervals

<div class="ref-description section level2">

Plots fitted coefficient functions with bootstrap pointwise intervals or
simultaneous bands.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
plot_inference(
  object,
  coefficient_set = "expert",
  type = c("pointwise", "simultaneous"),
  level = 0.95
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   object:

    A `vcmoe_bootstrap` object.

-   coefficient\_set:

    Coefficient set to plot: `"expert"` or `"gating"`.

-   type:

    Interval type passed to `confint()`.

-   level:

    Confidence level.

</div>

<div class="section level2">

## Value

A `ggplot` object.

</div>

</div>
