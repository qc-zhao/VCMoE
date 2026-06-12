<div id="main" class="col-md-9" role="main">

# Inspect VCMoE Parameterization Metadata

<div class="ref-description section level2">

Inspect VCMoE parameterization metadata.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_parameterization(object)
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

A named list describing the estimator convention used by the fit,
including kernel weights, local-linear basis scale, gating-logit
storage, dispersion block, label-alignment method, and optimization
controls.

</div>

<div class="section level2">

## Details

The package default is `"a1_epanechnikov_scaled"`: Epanechnikov density
weights `0.75 * (1 - t^2)_+ / h` with `t = (u - u0) / h`, scaled
local-linear slope storage, and centered gating logits. This helper
reports those conventions for reproducible model summaries.

</div>

</div>
