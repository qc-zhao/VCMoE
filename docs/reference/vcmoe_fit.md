<div id="main" class="col-md-9" role="main">

# Fit a varying-coefficient mixture-of-experts model

<div class="ref-description section level2">

Fit a varying-coefficient mixture-of-experts model

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_fit(
  formula,
  data,
  u,
  k = 2L,
  family = "gaussian",
  bandwidth = NULL,
  u_grid = NULL,
  control = list(),
  label = "align",
  parameterization = "a1_epanechnikov_scaled",
  u_scale = c("unit", "none"),
  engine = c("local_grid_em", "joint_path_em"),
  progress = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   formula:

    A formula of the form `y ~ expert_terms | gating_terms`. For grouped
    binomial data, use
    `cbind(success, failure) ~ expert_terms | gating_terms`. For
    Negative-Binomial count data, expert-side `offset(log_size_factor)`
    terms are supported.

-   data:

    A data frame.

-   u:

    Continuous index column name or numeric vector.

-   k:

    Number of mixture components. `k = 2` is the primary v0 target;
    `k = 3:10` are high-k candidate support and require diagnostics
    before being treated as stable.

-   family:

    Model family. `"gaussian"`, `"binomial"`, and `"negative-binomial"`
    are implemented.

-   bandwidth:

    Kernel bandwidth. If `NULL`, a Silverman-style default is used for
    `u`.

-   u\_grid:

    Grid where coefficient functions are estimated.

-   control:

    Named list overriding EM and label-alignment settings. Common
    entries: `maxit` (default 100), `tol` (default `1e-5`), `n_starts`
    (default 3) and `seed`. The local EM stops when the relative change
    of the penalized local objective, `|delta F| / (|F| + 1)`, falls
    below `tol`; with large samples this can stop well before the local
    maximum, so check `vcmoe_diagnostics()` and use `vcmoe_refine()`
    before inference. Negative-Binomial: `negbin_theta_ridge` and
    `negbin_theta_target` add a ridge on `log(theta)` towards the
    target, scaled by the kernel mass of the grid point
    (`negbin_theta_penalty_mass = "kernel"`, the default; `"component"`
    is the rule of VCMoE 0.2.0). Grid-parallel fitting (local-grid
    engine): `grid_subset` fits only the given grid indices and returns
    a `vcmoe_grid_part` to merge with `vcmoe_merge_grid_fits()`;
    `grid_seed = "per_grid"` reseeds each grid point with
    `seed + index`; `grid_warm_start = FALSE` drops the warm start from
    the previous grid point; `grid_init_locals` replaces the starts at
    given grid points (see `vcmoe_grid_neighbor_inits()`). The defaults
    reproduce the serial fit. Known labels (local-grid engine):
    `fixed_membership` gives each row of `data` a component, as indices
    in `1..k`, a factor (levels in component order) or an
    `nrow(data) x k` 0/1 matrix. Each local fit is then the
    kernel-weighted complete-data fit (component-wise experts and a
    multinomial-logit gate on the labels), with one start and no label
    alignment, and inference uses the complete-data local likelihood.
    Joint path: `joint_path_global_start` (default `TRUE`) makes the
    first start the posterior of a mixture of experts without `u`
    dependence, which keeps component labels consistent along `u`.

-   label:

    Label strategy. `"align"` uses exact global alignment for `k <= 6`
    and sequential assignment with ambiguity margins for `k >= 7`;
    `"global"` requests exact global alignment when feasible and falls
    back to the same sequential assignment path for `k >= 7`; `"greedy"`
    keeps the older one-step alignment.

-   parameterization:

    Estimator convention. The public package uses
    `"a1_epanechnikov_scaled"`: Epanechnikov density weights
    `K((u-u0)/h)/h` and the scaled local-linear basis `(u-u0)/h`.

-   u\_scale:

    How to transform `u` before fitting. The default `"unit"` maps
    complete-row `u` values to `[0, 1]`; `"none"` leaves `u` on the
    supplied scale. `bandwidth` and `u_grid` are interpreted on the
    transformed analysis scale.

-   engine:

    Fitting engine. The default `"local_grid_em"` preserves the original
    independent local-grid EM path. `"joint_path_em"` updates a shared
    observation-level responsibility path across grid points.

-   progress:

    Joint-path progress reporting. The default `NULL` is silent; use
    `TRUE` for messages or a single CSV file path for structured
    logging. Iteration frequency is controlled by
    `control$progress_every`.

</div>

<div class="section level2">

## Value

An object of class `vcmoe`, or a `vcmoe_grid_part` when
`control$grid_subset` is set.

</div>

<div class="section level2">

## Details

At each grid point, starts that met the stopping rule are preferred and
the one with the highest penalized local objective is kept; when no
start met it, the best finite start is kept with `converged = FALSE`.
`fit$diagnostics` records, per grid point, the stop reason, the last
objective change, the selected start and whether an unconverged start
reached a higher objective. With global label alignment, a grid point
relabelled against a confident match with its neighbour is flagged in
`alignment_override` with a warning; compare with `label = "greedy"`.

Rows with missing or non-finite response, covariates, or `u` are removed
consistently before fitting, with a warning. For single-trial Bernoulli
responses, the default gating ridge is strengthened to
`control$ridge = 1` unless the user explicitly supplies `control$ridge`;
grouped Binomial and other families keep the global default. Joint-path
traces record the sample-level nearest-grid log-likelihood as a
diagnostic criterion. The label-consistent updates do not guarantee that
this diagnostic is monotone; convergence is based on posterior and
parameter deltas instead.

</div>

</div>
