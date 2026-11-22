<div id="main" class="col-md-9" role="main">

# Refine a local-grid VCMoE fit

<div class="ref-description section level2">

Continues the local EM at selected grid points from the fitted solution,
with a tighter tolerance and more iterations. The relative stopping rule
of `vcmoe_fit()` (`tol = 1e-5` by default) corresponds to a large
absolute change of the local objective when the sample is large, so a
fit can stop before the local maximum is reached; the bread used by
`vcmoe_confband()` is then often not positive definite. Refining fixes
this without refitting from scratch. Each refined grid point keeps the
labels of the original fit; a refinement that would change them is
reported, not applied silently.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_refine(fit, grids = NULL, tol = 1e-08, maxit = 1000L)
```

</div>

</div>

<div class="section level2">

## Arguments

-   fit:

    A local-grid `vcmoe` fit that retains its data
    (`control$keep_data = TRUE`).

-   grids:

    Grid indices to refine. The default refines every grid point.

-   tol:

    Relative tolerance on the change of the penalized local objective.

-   maxit:

    Maximum number of EM iterations per grid point.

</div>

<div class="section level2">

## Value

The refined `vcmoe` fit. `diagnostics$refined` marks the refined grid
points, `diagnostics$refine_objective_gain` gives the change of the
penalized local objective and `diagnostics$refine_relabelled` flags grid
points whose refinement moved to a relabelled solution (kept with the
original labels).

</div>

</div>
