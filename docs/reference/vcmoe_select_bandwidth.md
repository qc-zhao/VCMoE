<div id="main" class="col-md-9" role="main">

# Select a VCMoE bandwidth by K-fold cross-validation

<div class="ref-description section level2">

Selects the kernel bandwidth for a VCMoE model using random K-fold
held-out predictive log-likelihood. The selected bandwidth is the
candidate with the largest held-out likelihood after ranking fully
successful candidates ahead of partial-failure candidates.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_select_bandwidth(
  formula,
  data,
  u,
  k = 2L,
  family = "gaussian",
  bandwidth_grid = NULL,
  folds = 5L,
  u_grid = NULL,
  control = list(),
  label = "align",
  parameterization = "a1_epanechnikov_scaled",
  u_scale = c("unit", "none"),
  seed = NULL,
  refit = TRUE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   formula:

    A formula of the form `y ~ expert_terms | gating_terms`.

-   data:

    A data frame.

-   u:

    Continuous index column name or numeric vector.

-   k:

    Number of mixture components. Values from `2` through `10` are
    supported.

-   family:

    Model family. `"gaussian"`, `"binomial"`, and `"negative-binomial"`
    are supported.

-   bandwidth\_grid:

    Candidate bandwidth values. If `NULL`, uses multiples of the default
    bandwidth.

-   folds:

    Number of random cross-validation folds.

-   u\_grid:

    Grid where coefficient functions are estimated.

-   control:

    Named list passed to `vcmoe_fit()`.

-   label:

    Label strategy passed to `vcmoe_fit()`.

-   u\_scale:

    `u` scaling strategy passed to `vcmoe_fit()`. The default selects
    bandwidths on the unit-scaled analysis domain.

-   parameterization:

    Estimator convention passed to `vcmoe_fit()` for each CV fold and
    the optional final refit.

-   seed:

    Optional random seed for fold assignment and, when `control$seed` is
    absent, deterministic CV refits.

-   refit:

    Whether to refit the final model on all data using the selected
    bandwidth.

</div>

<div class="section level2">

## Value

An object of class `vcmoe_bandwidth_selection` with fields
`best_bandwidth`, `cv_summary`, `cv_details`, `cv_folds`, `fit`, and
`settings`.

</div>

<div class="section level2">

## Details

The default candidate grid is the current Silverman-style default
bandwidth multiplied by `c(0.5, 0.75, 1, 1.25, 1.5, 2)`. Fold assignment
is made only among complete rows that `vcmoe_fit()` would keep.

For Gaussian models, validation scoring uses
`log sum_c pi_c Normal(y | mu_c, sigma_c)`. For Binomial models, scoring
uses `log sum_c pi_c Binomial(success | trials, p_c)`. Binomial
Bernoulli and grouped `cbind(success, failure)` responses are supported.
For Negative-Binomial models, scoring uses
`log sum_c pi_c NB(y | mu_c, theta_c)`.

Bandwidth selection supports `k = 2:10` for Gaussian, Binomial, and
Negative-Binomial models. High-k candidates use the same held-out
predictive likelihood scoring and should be interpreted together with
the returned fit diagnostics. The selected object records the fitting
parameterization and `u` scaling strategy in `settings$parameterization`
and `settings$u_scale`.

</div>

</div>
