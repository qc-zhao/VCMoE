<div id="main" class="col-md-9" role="main">

# Fit a varying-coefficient mixture-of-experts model

<div class="ref-description section level2">

Fits a Gaussian, Binomial, or Negative-Binomial VCMoE model by
local-linear EM and aligns component labels across the
coefficient-function grid.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_fit(formula, data, u, k = 2L, family = "gaussian",
  bandwidth = NULL, u_grid = NULL, control = list(), label = "align",
  parameterization = "a1_epanechnikov_scaled", u_scale = c("unit", "none"))
```

</div>

</div>

<div class="section level2">

## Arguments

-   formula:

    A formula of the form `y ~ expert_terms | gating_terms`. For grouped
    Binomial data, use
    `cbind(success, failure) ~ expert_terms | gating_terms`. For
    Negative-Binomial count data, use expert-side
    `offset(log_size_factor)` for library-size or size-factor offsets.

-   data:

    A data frame.

-   u:

    Continuous index column name or numeric vector.

-   k:

    Number of mixture components. Values 2 through 10 are accepted.
    High-k fits are candidate support and require diagnostics.

-   family:

    Model family. `"gaussian"`, `"binomial"`, and `"negative-binomial"`
    are implemented.

-   bandwidth:

    Kernel bandwidth. If `NULL`, a default is used.

-   u\_grid:

    Grid where coefficient functions are estimated.

-   control:

    Named list overriding EM and label-alignment settings.

-   label:

    Label strategy. `"align"` uses exact global alignment for `k <= 6`
    and sequential assignment for `k >= 7`. `"global"` requests exact
    global alignment when feasible and falls back to the same sequential
    assignment path for `k >= 7`. `"greedy"` keeps the older one-step
    alignment.

-   u\_scale:

    How to transform `u` before fitting. The default `"unit"` maps
    complete-row `u` values to `[0, 1]`; `"none"` leaves `u` on the
    supplied scale. `bandwidth` and `u_grid` are interpreted on the
    transformed analysis scale.

-   parameterization:

    Estimator convention. The public package uses
    `"a1_epanechnikov_scaled"`: Epanechnikov density weights
    `0.75 * (1 - t^2)_+ / h` with `t = (u - u0) / h` and stores
    local-linear slopes on the scaled `(u - u0) / h` basis.

</div>

<div class="section level2">

## Value

An object of class `vcmoe`.

</div>

<div class="section level2">

## Details

Rows with missing or non-finite response, covariates, or `u` are removed
consistently before fitting, with a warning.

By default, local EM uses Epanechnikov density kernel weights and a
scaled local-linear basis. Complete-row `u` values are mapped to
`[0, 1]` before fitting unless `u_scale = "none"`. Fitted objects store
both the analysis-scale grid and the original-scale grid metadata.

For `family = "gaussian"`, the expert mean and log-standard-deviation
blocks are both local-linear inside each fit. The variance block uses
`log(sigma_ic) = delta_c0(u0) + delta_c1(u0) * (u_i - u0) / h`;
`coef(fit, "sigma")` returns `exp(delta_c0(u0))` and
`coef(fit, "sigma_slope")` returns `delta_c1(u0)`.

For `family = "binomial"`, Bernoulli responses must be 0/1 values and
grouped responses must use non-negative finite whole-number
success/failure counts with positive row totals. Binomial expert
coefficients are on the logit success-probability scale. Binomial
`k >= 3` fits are accepted as experimental stress support and should not
be treated as stable inference-ready support. Binomial expert logistic
M-steps use `control$binomial_ridge` with default value `1`, and
`control$binomial_structured_starts = TRUE` uses structured local starts
to improve single-trial Bernoulli stability. For single-trial Bernoulli
responses only, the gating ridge default is also strengthened to
`control$ridge = 1` unless the user explicitly supplies `control$ridge`;
grouped Binomial responses keep the global default.

For `family = "negative-binomial"`, responses must be finite
non-negative whole-number counts. Expert coefficients are on the log
mean count scale. Use `offset(log_size_factor)` in the expert formula to
account for library size or size factors. Gating-side offsets are not
supported in v0. `coef(fit, "theta")` returns component-specific NB size
parameters.

The default label strategy performs local EM at each grid point, then
applies a post-processing dynamic-programming alignment over the full
grid. Transition costs use the previous local-linear slope to predict
the next coefficient value, plus gating, posterior, and, for Gaussian
fits, variance consistency terms. This improves label path tracking
without changing the local EM estimating equations.

For `k <= 6`, `label = "align"` uses exact derivative-aware global
alignment. For `k >= 7`, `label = "align"` and `label = "global"` use
sequential pairwise assignment alignment to avoid factorial permutation
growth and record the best-vs-second-best assignment margin for
ambiguity diagnostics.

</div>

</div>
