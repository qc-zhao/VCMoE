<div id="main" class="col-md-9" role="main">

# Merge grid-subset VCMoE fits

<div class="ref-description section level2">

Combines `vcmoe_grid_part` objects returned by `vcmoe_fit()` with
`control$grid_subset` (for example one per cluster job) into an ordinary
`vcmoe` fit. Every grid point must be covered. When several parts cover
the same grid point (for example a cold-start pass and a neighbour
warm-start pass), a converged local fit is preferred, then the one with
the higher penalized local objective. Label alignment is then run over
the whole grid as in a serial fit.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_merge_grid_fits(parts, label = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   parts:

    A list of `vcmoe_grid_part` objects from the same data and settings.

-   label:

    Optional label-alignment override: `"align"`, `"global"` or
    `"greedy"`. The default keeps the alignment of the original call.

</div>

<div class="section level2">

## Value

A `vcmoe` object; `diagnostics$grid_source` records which part supplied
each grid point.

</div>

</div>
