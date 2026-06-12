<div id="main" class="col-md-9" role="main">

# Inspect VCMoE Local-Linear Slopes On The Scaled Basis

<div class="ref-description section level2">

Inspect VCMoE local-linear slopes on the scaled basis.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_scaled_slopes(object, type = c("expert", "gating"), bandwidth = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   object:

    A `vcmoe` object.

-   type:

    Coefficient block, either `"expert"` or `"gating"`.

-   bandwidth:

    Optional bandwidth recorded in the returned attributes. Defaults to
    the fitted bandwidth.

</div>

<div class="section level2">

## Value

An array with the same dimensions as the stored slope block.

</div>

<div class="section level2">

## Details

VCMoE stores slopes on the scaled local-linear basis `(u - u0) / h`.
This helper returns the stored scaled-basis slope block.

</div>

</div>
