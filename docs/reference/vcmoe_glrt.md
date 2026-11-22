<div id="main" class="col-md-9" role="main">

# Generalized likelihood-ratio test for VCMoE coefficient variation

<div class="ref-description section level2">

Generalized likelihood-ratio test for VCMoE coefficient variation

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_glrt(
  fit,
  data,
  test = c("coefficient", "constant_block", "constant_all"),
  coefficient_set = c("expert", "gating", "sigma", "theta"),
  component = NULL,
  term = NULL,
  calibration = c("none", "bootstrap", "analytic_epanechnikov", "both",
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
    `"constant_block"` tests all expert or gating functions jointly;
    `"constant_all"` tests all fitted coefficient functions jointly.

-   coefficient\_set:

    Coefficient block for coefficient-specific or block-constant tests.

-   component:

    Component label or index for coefficient-specific tests.

-   term:

    Term name for coefficient-specific tests.

-   calibration:

    Calibration method. The default `"none"` returns the statistic
    without attaching a reference distribution.
    `"analytic_epanechnikov"` uses the Epanechnikov modified chi-square
    approximation; `"bootstrap"` uses parametric bootstrap calibration;
    `"both"` reports both.

-   B:

    Number of bootstrap calibration replicates.

-   seed:

    Optional random seed.

-   control:

    Controls for constrained null optimization and diagnostics:
    `lambda_type` (`"observation"`, the default; `"fit_objective"`, the
    scale used by VCMoE 0.2.0; or `"grid_sum"`), `strict` (default
    `FALSE`), `C_const` (default 1), `null_projection` for joint-path
    nulls (`"uniform"`, the default, or `"sample_weighted"`, the rule of
    VCMoE 0.2.0), `maxit` and `reltol` for the null optimizer, and
    `lambda_tol`.

-   refit\_control:

    Controls overriding bootstrap full-model refits.

-   verbose:

    Whether to message bootstrap progress.

</div>

<div class="section level2">

## Value

A `vcmoe_glrt` object.

</div>

<div class="section level2">

## Details

**Statistic.** `lambda` is the difference between the full and null
observation-level mixture log-likelihoods,
`sum_i log sum_k pi_k(u_i) f_k(y_i; u_i)`, with each coefficient path
evaluated at the observation's own `u_i` by the linear interpolation
used by `predict.vcmoe()`. This is the likelihood of the Wilks-type
result for local likelihood ratio tests. `lambda_type = "fit_objective"`
instead uses the objective the engine maximises (the sum of grid-point
local likelihoods for local-grid fits, the nearest-grid observation
likelihood for joint-path fits). A sum over grid points counts each
observation about `1 / (grid spacing)` times, so analytic p-values on
that scale are far too small; it is kept only to reproduce earlier
results. The fitted-objective ratio is always reported as
`lambda_fit_objective`, and a negative value blocks the test as a
null-optimizer failure. A negative observation-level ratio is a valid,
non-significant outcome.

**Analytic calibration.** The statistic `rK * max(lambda, 0)` is
compared with a chi-square distribution with
`delta = rK * p_test * C_const * |U| * (K(0) - int K^2 / 2) / h` degrees
of freedom, where `p_test` is the number of coefficient functions
constrained under the null: 1 for one coefficient function, `k - 1`
contrasts for a gating term, and the whole block for block tests. The
approximation is intended for independent observations and joint-path
fits; for local-grid fits it has not been validated and
`settings$theory_status` says so. With clustered observations (for
example cells from the same donor), use the cluster-robust path tests of
`vcmoe_confband()` instead.

**Null fits.** Local-grid fits use a constrained BFGS null optimizer.
Joint-path fits use a grid-projected null: after every M-step each
constrained coefficient path is replaced by its average over the grid
points (equal weights by default; with
`null_projection = "sample_weighted"`, weights proportional to the
number of observations assigned to each nearest grid point), and
constrained local slopes are set to zero. An optimizer error, or a null
fit that beats the full fit on the fitted objective, blocks the test and
returns an `NA` p-value. A null fit that stops at its iteration limit
gives a p-value with the warning `null_optimizer_not_fully_converged`
and `null_converged = FALSE`; `strict = TRUE` blocks it, together with
every other warning.

**Bootstrap.** Bootstrap refits preserve the engine of `fit`. Each refit
numbers its components independently, so a component-specific test
follows the matched reference component into the refit; for a gating
contrast with `k >= 3` both components of the contrast are followed.
Draws whose null optimizer failed, or whose null beats the full refit,
are excluded; draws whose null fit only reached its iteration limit are
kept. `bootstrap_counts` reports the requested, successful and excluded
draws.

</div>

</div>
