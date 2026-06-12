<div id="main" class="col-md-9" role="main">

# Plot VCMoE fit diagnostics

<div class="ref-description section level2">

Plots convergence, posterior entropy, component proportions, effective
local sample size, and label ambiguity flags over the coefficient grid.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
plot_diagnostics(object)
```

</div>

</div>

<div class="section level2">

## Arguments

-   object:

    A `vcmoe` object.

</div>

<div class="section level2">

## Value

A `ggplot` object.

</div>

<div class="section level2">

## Details

This plot is intended as a first real-data sanity check before
interpreting coefficient functions. Ambiguity or non-convergence at many
grid points should be treated as evidence that the fitted component
labels or coefficient paths need closer review.

</div>

</div>
