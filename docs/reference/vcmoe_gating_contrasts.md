<div id="main" class="col-md-9" role="main">

# Report Identifiable VCMoE Gating Contrasts

<div class="ref-description section level2">

Report identifiable gating contrasts.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_gating_contrasts(object, baseline = NULL, scaled = FALSE)
```

</div>

</div>

<div class="section level2">

## Arguments

-   object:

    A `vcmoe` object.

-   baseline:

    Component used as the contrast baseline. By default, `k = 2` uses
    component 2 and `k > 2` uses component 1.

-   scaled:

    If `TRUE`, slope contrasts are converted to the scaled local-linear
    basis `(u - u0) / h`.

</div>

<div class="section level2">

## Value

A data frame with one row per grid point, contrast, term, and block.

</div>

<div class="section level2">

## Details

VCMoE stores gating coefficients as centered logits, so the absolute
level of all component logits is not identifiable. Interpretable gating
effects are component contrasts such as component 1 versus component 2
for `k = 2`, or component 1 and component 2 versus component 3 for
`k = 3` comparisons when `baseline = 3`.

</div>

</div>
