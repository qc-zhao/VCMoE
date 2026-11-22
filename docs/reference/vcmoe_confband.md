<div id="main" class="col-md-9" role="main">

# Confidence bands for a VCMoE fit

<div class="ref-description section level2">

Pointwise intervals and finite-grid simultaneous bands for the
coefficient functions of a fitted VCMoE model, from the local sandwich
covariance.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_confband(
  fit,
  data = NULL,
  level = 0.95,
  type = c("pointwise", "simultaneous"),
  coefficient_set = c("expert", "gating", "sigma", "theta"),
  strict = TRUE,
  control = list(),
  cluster = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   fit:

    A `vcmoe` fit with `k = 2:10`.

-   data:

    Optional original data frame. The current implementation uses the
    data stored in `fit$fitted`; refit with `keep_data = TRUE` if
    needed.

-   level:

    Confidence level.

-   type:

    Interval columns to expose as `lower` and `upper`.

-   coefficient\_set:

    Coefficient blocks to return.

-   strict:

    Whether weak local fits should return blocked intervals.

-   control:

    Optional inference controls. Without `cluster`, the covariance is
    HC0 and the simultaneous band uses the analytic Epanechnikov path
    critical value;
    `simultaneous_method = "multiplier_coefficient_path"` calibrates
    each coefficient over the complete fitted grid by a Gaussian
    multiplier bootstrap instead. With `cluster`,
    `covariance_adjustment` may be `"CR0"`, `"CR1"`, `"CR2"` (default)
    or `"CR3"`, and `simultaneous_method` `"wild_cluster_path"`
    (default) or `"multiplier_coefficient_path"`. `multiplier_B`
    (default 199) sets the number of bootstrap draws and `seed` their
    random seed. `path_tests = TRUE` adds the zero and constancy path
    tests described below; `path_sets`, a named list of grid-index
    vectors, adds bands and tests over sub-paths (for example the
    interior grid points).

-   cluster:

    Optional atomic vector of independent-cluster identifiers (for
    example donors), one per retained fitted row, in fitted-row order
    (not the original unfiltered order). If named, names must match the
    retained design-matrix row names exactly. Missing IDs are rejected.

</div>

<div class="section level2">

## Value

A `vcmoe_confband` object with interval and diagnostic data frames; with
`cluster`, `cluster_support` reports per-component cluster support, and
with `path_tests = TRUE`, `tests` holds the path-test statistics and
p-values.

</div>

<div class="section level2">

## Details

**Independent observations.** HC0 sandwich intervals with normal
critical values. For `engine = "joint_path_em"`, the covariance is the
observed local-likelihood sandwich evaluated at the joint-path estimate;
it does not include shared-path, label-selection, or finite-grid
cross-grid responsibility uncertainty, so joint-path convergence and the
returned score-imbalance diagnostics should be inspected.

**Clustered observations.** With `cluster`, kernel-weighted scores are
summed within clusters. Cluster inference supports unrestricted
local-grid fits that retain their data. Coefficients and the penalized
bread are unchanged. CR0 uses normal critical values. CR1 scales by
`G / (G - 1)` with `t(G - 1)`, for `G` clusters with positive kernel
weight at the grid point. CR2 (bias-reduced linearization in score
space) and CR3 (one-step cluster jackknife) use Bell-McCaffrey degrees
of freedom under a working model of independent observations within
clusters; both use cluster Hessian blocks. The default simultaneous band
studentizes Webb-weighted wild cluster score draws with the same
adjustment at every grid point; one weight per cluster is shared by all
grid points and coefficients. When the number of distinct Webb weight
vectors, `6^G`, is at most `multiplier_B`, all of them are used.
Clusters must be independent; smoothing bias and bandwidth selection are
not corrected. With fewer than 20 clusters, or fewer than 10 effective
clusters, at a grid point a warning flags the intervals as sensitivity
analyses.

**Paths and tests.** A band covers one coefficient's complete fitted
grid, not a continuous domain or all coefficients jointly. A blocked
grid point blocks that coefficient's complete-path band, not its other
pointwise intervals; a band over a sub-path in `path_sets` remains
available. With `path_tests = TRUE` each complete path also gets two
sup-t tests from the same bootstrap draws as its band: `beta(u) = 0`
(the band-inversion test) and `beta(u)` constant, which contrasts each
grid estimate with the path mean and studentizes the contrast with the
same small-cluster adjustment. These tests are the recommended constancy
tests for clustered data; the likelihood-ratio test of `vcmoe_glrt()`
assumes independent observations.

</div>

</div>
