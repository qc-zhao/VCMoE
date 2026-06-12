<div id="main" class="col-md-9" role="main">

# Coefficient-specific GLRT for VCMoE coefficient variation

<div class="ref-description section level2">

Fits a constrained null under the same local objective and computes a
generalized likelihood-ratio statistic against the supplied VCMoE fit.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_glrt(
  fit,
  data,
  test = c("coefficient", "constant_all"),
  coefficient_set = c("expert", "gating", "sigma", "theta"),
  component = NULL,
  term = NULL,
  calibration = c("analytic_epanechnikov", "bootstrap", "both", "none",
    "parametric_bootstrap"),
  B = 200L,
  seed = NULL,
  control = list(),
  refit_control = list(),
  verbose = FALSE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   fit:

    A `vcmoe` fit.

-   data:

    Original data frame used to fit `fit`.

-   test:

    Test type. `"coefficient"` tests one coefficient function;
    `"constant_all"` tests all fitted coefficient functions jointly.

-   coefficient\_set:

    Coefficient block for coefficient-specific tests.

-   component:

    Component label or index for coefficient-specific tests.

-   term:

    Term name for coefficient-specific tests.

-   calibration:

    Calibration method. `"analytic_epanechnikov"` uses the Epanechnikov
    modified chi-square calibration; `"bootstrap"` uses parametric
    bootstrap calibration; `"both"` reports both.
    `"parametric_bootstrap"` is accepted as a backwards-compatible alias
    for `"bootstrap"`.

-   B:

    Number of bootstrap calibration replicates.

-   seed:

    Optional random seed.

-   control:

    Controls for the constrained null optimizer.

-   refit\_control:

    Controls overriding bootstrap full-model refits.

-   verbose:

    Whether to message progress.

</div>

<div class="section level2">

## Value

A `vcmoe_glrt` object with the observed `lambda`, analytic statistic,
null fit, optional bootstrap replicate summary, and calibrated p-value
when available.

</div>

<div class="section level2">

## Details

The primary path is coefficient-specific. The selected coefficient
function is constrained to be constant in `u`: its local-linear slope is
fixed to zero and its intercept is shared across grid points. The null
is re-optimized rather than obtained by post-hoc averaging.

Analytic Epanechnikov calibration requires Epanechnikov density weights,
scaled local-linear basis, and unit-scaled `u`. It reports
`lambda = ell_full - ell_null`, `analytic_statistic = rK * lambda`, and
a modified chi-square `analytic_p_value`. The diagnostic
`lrt_statistic = 2 * lambda` is reported separately and is not used for
the analytic p-value. Ridge penalties may be used to stabilize fitting
but are excluded from the GLRT likelihood ratio.

`vcmoe_glrt()` supports fitted `k = 2:10` models for both
coefficient-specific and `"constant_all"` nulls. For `k > 2`, gating
coefficient tests use identifiable baseline contrasts, e.g.
`component3_vs_component1`. Reported p-values should be interpreted
together with fit convergence, label ambiguity, component proportion,
and null-optimizer diagnostics.

</div>

</div>
