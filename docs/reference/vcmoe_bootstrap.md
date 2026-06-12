<div id="main" class="col-md-9" role="main">

# Parametric bootstrap inference for a VCMoE fit

<div class="ref-description section level2">

Runs parametric bootstrap inference for a fitted Gaussian, Binomial, or
Negative-Binomial VCMoE model with `k = 2:10`.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_bootstrap(
  fit,
  data,
  u = NULL,
  B = 200L,
  coefficient_set = c("expert", "gating"),
  seed = NULL,
  control = list(),
  min_successful = max(20L, ceiling(0.5 * B)),
  keep_fits = FALSE,
  verbose = FALSE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   fit:

    A fitted `vcmoe` object. Bootstrap inference supports `k = 2:10`.

-   data:

    Original data frame used to fit `fit`. The function resamples from
    `data[fit$rows_used, ]`.

-   u:

    Optional original `u` values or column name. If `NULL`, the stored
    `u` column from `fit` is reused when available.

-   B:

    Number of parametric bootstrap replicates.

-   coefficient\_set:

    Coefficient sets to store and summarize: `"expert"`, `"gating"`, or
    both.

-   seed:

    Optional random seed.

-   control:

    Named list passed to bootstrap refits. Bandwidth is not reselected
    inside bootstrap v0.

-   min\_successful:

    Minimum number of successful replicates expected for reliable
    inference. The object is returned when at least two replicates
    succeed, but a warning is recorded below this threshold.

-   keep\_fits:

    Whether to store successful bootstrap fit objects.

-   verbose:

    Whether to print replicate progress messages.

</div>

<div class="section level2">

## Value

An object of class `vcmoe_bootstrap` with fields `fit`, `replicates`,
`replicate_summary`, `alignment_summary`, `settings`, `warnings`, and
optionally `fits`.

</div>

<div class="section level2">

## Details

For Gaussian fits, each bootstrap data set draws a latent component from
the fitted gating probabilities and then draws the response from the
selected component Normal distribution. For Binomial fits, each
bootstrap data set draws success counts from the selected component
success probability. For Negative-Binomial fits, each bootstrap data set
draws counts from the selected component mean and theta. Bernoulli and
grouped `cbind(success, failure)` response formats are preserved for
Binomial fits.

Each bootstrap replicate is refit with the same formula, family, number
of components, bandwidth, `u_grid`, and label strategy as the reference
fit. After the usual within-grid label alignment, one global component
permutation matches the bootstrap coefficient paths back to the
reference fit. Ambiguous bootstrap-to-reference matches are recorded in
`alignment_summary`. Exact permutation matching is used for small `k`;
assignment-based matching is used when exhaustive permutation is
infeasible.

Binomial expert coefficients and intervals are on the logit coefficient
scale. Negative-Binomial expert coefficients and intervals are on the
log mean count scale.

</div>

</div>
