<div id="main" class="col-md-9" role="main">

# Cluster-robust gene screen for Negative-Binomial expert effects

<div class="ref-description section level2">

Screens candidate covariates (typically genes) for an association with
the response inside the latent components, while the fitted VCMoE (the
null model) is held fixed. For each candidate, every grid point of the
null fit is augmented by the covariate, entering each expert's local
linear predictor (the tested coefficients) and each gating contrast
(nuisance), with coefficients starting at zero. A one-step estimate of
the expert coefficients is studentized with the CR2 cluster-robust
covariance and calibrated over the grid by a studentized wild cluster
bootstrap, as in the path tests of `vcmoe_confband()`. The null model is
fitted once, so all candidates share its moments and the screen scales
to thousands of genes.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_gene_screen(
  fit,
  X,
  cluster,
  profile_genes = NULL,
  grids = NULL,
  B = 9999L,
  seed = NULL,
  adjust = c("BH", "minP", "both"),
  chunk_size = 200L,
  control = list(),
  profile_control = list()
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   fit:

    The null model: a single-response Negative-Binomial local-grid
    `vcmoe` fit that retains its data.

-   X:

    Numeric matrix of candidate covariates, one row per retained fitted
    row (in fitted-row order) and one uniquely named column per
    candidate. Columns are centred internally.

-   cluster:

    Independent-cluster identifiers (for example donors), one per
    retained fitted row. Required.

-   profile\_genes:

    Optional column names of `X` to screen after profiling their gating
    effect: the candidate first enters the gating of the null model and
    every screened grid point is refitted by local EM, giving the
    restricted local fit under the null (no expert effect, gating effect
    free). Use this for strong markers of the latent components, whose
    gating effect is far from zero; the one-step screen holds it at
    zero.

-   grids:

    Grid indices that form the tested path. Defaults to every grid
    point. A candidate is reported with an `NA` p-value when any of
    these grid points is blocked for it.

-   B:

    Number of wild bootstrap draws (Webb weights). One set of weights is
    generated per call and shared by every candidate, chunk and profiled
    refit, including when `seed = NULL`. When `6^G <= B` for `G`
    clusters, all `6^G` weight vectors are used.

-   seed:

    Optional seed for the bootstrap weights.

-   adjust:

    Multiplicity adjustment over the full candidate set: `"BH"`
    (Benjamini-Hochberg, default), `"minP"` (single-step Westfall-Young
    min-p from the shared draws) or `"both"`. Candidates with an `NA`
    p-value stay in the family size of the BH adjustment.

-   chunk\_size:

    Number of candidates processed together; smaller values use less
    memory and do not change the results.

-   control:

    Optional inference controls passed to the local moment computations
    (for example `strict`, `max_hessian_dim`).

-   profile\_control:

    List with `maxit` (default 300) and `tol` (default `1e-9`) for the
    profiled refits.

</div>

<div class="section level2">

## Value

A list of class `vcmoe_gene_screen` with `genes` (one row per candidate:
sup-t statistic, raw and adjusted p-values, maximising component and
grid, per-component statistics and directions, per-SD effects, grid
usage and status), `pointwise` (per grid and coefficient one-step
estimates, standard errors and Bell-McCaffrey degrees of freedom),
`grids` (grid-level status), `profile` (profiled refits), `draws`
(shared bootstrap maxima) and `settings`.

</div>

<div class="section level2">

## Details

The screen tests expert effects, so the null model should explain
component membership well: include covariates that distinguish the
latent components (for example known subtype markers) in its gating
formula. Otherwise strong markers can leak into the expert test. Results
describe associations on the log-mean scale of each expert;
`effect_per_sd` gives the one-step level coefficient at the maximising
grid point per standard deviation of the candidate.

</div>

</div>
