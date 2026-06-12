<div id="main" class="col-md-9" role="main">

# Convert VCMoE Slopes To A Scaled Local-Linear Basis

<div class="ref-description section level2">

Convert VCMoE local-linear slopes to scaled-basis slopes.

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

    Optional bandwidth used for conversion. Defaults to the fitted
    bandwidth.

</div>

<div class="section level2">

## Value

An array with the same dimensions as the stored slope block.

</div>

<div class="section level2">

## Details

VCMoE stores slopes directly for the scaled basis `(u - u0) / h`, so
this helper returns stored slopes unchanged. Legacy Gaussian/raw fits
store slopes for `u - u0`; for those fits the equivalent scaled-basis
slope is `h * slope_raw`. This helper performs the needed conversion
without changing the fitted object.

</div>

</div>
