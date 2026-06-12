<div id="main" class="col-md-9" role="main">

# Summarize VCMoE fit diagnostics

<div class="ref-description section level2">

Returns a compact diagnostic table for reviewing whether a fitted VCMoE
model is reliable enough to interpret.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_diagnostics(object)
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

A data frame with one row per coefficient grid point.

</div>

<div class="section level2">

## Details

The table includes convergence status, iterations, local log-likelihood,
local-weighted posterior entropy, label ambiguity flags, alignment
margin, effective local sample size, local-weighted component posterior
proportions, and Binomial expert optimizer diagnostics when available.

Posterior entropy and component proportions use the same local kernel
weights as the fitted grid point when the fit retains training data. If
the fit was created with `control$keep_data = FALSE`, component
proportions fall back to unweighted posterior means and effective local
sample size is `NA`.

</div>

</div>
