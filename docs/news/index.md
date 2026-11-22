<div id="main" class="col-md-9" role="main">

# Changelog

<div class="section level2">

## VCMoE 0.3.0

<div class="section level3">

### Fitting and convergence

-   The local-grid EM now stops on the penalized local objective it
    ascends (`loglik - penalty`), not on the unpenalized log-likelihood,
    which can stop on a flat stretch while the objective still climbs.
    Starts that met the stopping rule are preferred, and among them the
    one with the highest penalized objective is kept; if none met it,
    the best finite start is kept with `converged = FALSE`. `loglik`
    still reports the unpenalized log-likelihood. Default fits can
    therefore differ slightly from 0.2.0.
-   New `vcmoe_refine()` continues the local EM from a fitted solution
    with a tighter tolerance, keeping the labels. The default relative
    `tol = 1e-5` can stop well before the local maximum in large
    samples; refining before inference usually makes the sandwich bread
    positive definite.
-   `vcmoe_diagnostics()` reports per grid point the stop reason, the
    last objective change, the selected start, whether an unconverged
    start reached a higher objective, label-alignment overrides and
    refinement.
-   Global label alignment that relabels a grid point against a
    confident match with its neighbour is flagged (`alignment_override`)
    with a warning.
-   Grid-parallel fitting: `control$grid_subset`, `grid_seed`,
    `grid_warm_start` and `grid_init_locals`, with
    `vcmoe_merge_grid_fits()` and `vcmoe_grid_neighbor_inits()`. Merged
    parts reproduce the serial fit exactly under per-grid seeds without
    warm starts. Merges can switch to `label = "greedy"`.
-   Negative-Binomial theta ridges are scaled by the kernel mass of the
    grid point (`negbin_theta_penalty_mass = "kernel"`), so the fitted
    theta solves the estimating equation the sandwich differentiates.
    Objects from 0.2.0 are read with their original rule; strict
    inference on such fits with a theta ridge is blocked with a request
    to refit.
-   Known-label fits through `control$fixed_membership`.
-   Joint-path fits start by default from a mixture of experts without
    `u` dependence (`joint_path_global_start = TRUE`), which keeps
    component labels consistent along `u`. Joint-path fits can therefore
    differ from 0.2.0.

</div>

<div class="section level3">

### Generalized likelihood-ratio tests

-   `vcmoe_glrt()` now compares observation-level mixture
    log-likelihoods, with each coefficient path evaluated at the
    observation’s own `u`. The 0.2.0 local-grid statistic summed
    grid-point local likelihoods, which counts each observation about
    `1 / (grid spacing)` times; its analytic p-values were far too
    small, and results computed with it should be recomputed.
    `control$lambda_type = "fit_objective"` reproduces the 0.2.0 scale.
-   Analytic calibration counts only the coefficient functions
    constrained under the null: the default `C_const` is now 1 (it was
    2).
-   `strict = FALSE` is the default. Optimizer failures still block the
    test; a null fit that stopped at its iteration limit gives a p-value
    with a warning and `null_converged = FALSE`.
-   Joint-path null fits project constrained paths onto their
    equal-weight grid average by default
    (`control$null_projection = "uniform"`), the setting used in the
    validation of the analytic calibration;
    `null_projection = "sample_weighted"` restores the 0.2.0 weighting
    by nearest-grid observation counts.
-   Bootstrap refits follow the tested component, and for `k >= 3` both
    components of a gating contrast. Draws whose null optimizer failed
    or whose null beats the full refit are excluded; `bootstrap_counts`
    reports the requested, successful and excluded draws.
-   GLRT and reduced-fit metadata changed: the criterion-matching entry
    was removed from `settings`, and `theory_status`,
    `calibration_status` and `statistic_criterion` take new values that
    describe the statistic used. Joint-path confidence-band metadata
    uses `"joint_path_observed_local_likelihood_plugin"`.

</div>

<div class="section level3">

### Clustered observations

-   `vcmoe_confband(fit, cluster = id)` gives cluster-robust intervals
    for local-grid fits: CR2 with Bell-McCaffrey degrees of freedom by
    default (CR0, CR1 and CR3 are available), and finite-grid
    simultaneous bands from a studentized wild cluster bootstrap with
    Webb weights. Diagnostics report active and effective clusters,
    degrees of freedom and block reasons per grid point.
    `control$path_tests = TRUE` adds zero and constancy path tests, and
    `control$path_sets` bands over sub-paths.
-   New `vcmoe_gene_screen()` screens candidate covariates for expert
    associations in Negative-Binomial fits, with CR2 studentization, a
    shared wild cluster bootstrap, optional gating profiling of strong
    markers and Benjamini-Hochberg or min-p adjustment over the full
    candidate set.
-   `predict()` bounds Negative-Binomial `theta` to
    `[negbin_theta_min, negbin_theta_max]` before interpolating, as the
    likelihood does, so a stored unbounded value no longer gives `NaN`
    predictions.
-   Strict inference now blocks only the grid points that did not
    converge. Negative-Binomial eta-bound checks use only observations
    with positive kernel weight at the grid point.

</div>

</div>

<div class="section level2">

## VCMoE 0.2.0

CRAN release: 2026-07-15

<div class="section level3">

### Fitting engines

-   Adds stable single-response joint-path EM through
    `vcmoe_fit(..., engine = "joint_path_em")`.
-   Keeps `engine = "local_grid_em"` as the default and preserves the
    0.1.0 local-grid numerical path when `engine` is omitted.
-   Adds dense-grid cost guardrails, nearest-grid observation
    assignment, joint-path convergence traces, and optional progress
    logging.

</div>

<div class="section level3">

### Inference and selection

-   Makes analytic confidence bands available through the common
    fitted-object interface for both engines. Joint-path bands identify
    their local-likelihood sandwich plug-in scope and expose
    score-imbalance diagnostics.
-   Makes parametric bootstrap refits and bandwidth selection
    engine-aware.
-   Adds sample-weighted projected joint-path constrained-null fitting
    for GLRT and `vcmoe_fit_reduced()` while retaining the established
    local-grid BFGS null optimizer. Joint-path GLRT uses sample-level
    nearest-grid likelihood contributions and blocks failed or
    nonconverged null fits.
-   Changes the GLRT calibration default to `"none"`. The statistic is
    available without a reference distribution; analytic Epanechnikov
    calibration is an explicit approximation, and parametric bootstrap
    calibration is engine-preserving.
-   Adds criterion, projection, constrained-MLE, calibration, and
    theory-status metadata to GLRT results. Analytic confidence bands
    now state whether their covariance target matches the fitted
    estimator and which shared-path and label uncertainties are
    excluded.

</div>

<div class="section level3">

### Scope

-   Version 0.2.0 remains a single-response package for Gaussian,
    Binomial, and Negative-Binomial models with 2 through 10 components.

</div>

</div>

</div>
