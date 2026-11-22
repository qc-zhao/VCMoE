<div id="main" class="col-md-9" role="main">

# Neighbour warm starts for a second grid-parallel pass

<div class="ref-description section level2">

For each requested grid point, returns the best local fits across
`parts` (ranked as in `vcmoe_merge_grid_fits()`) at the neighbouring
grid points, ready to pass as `control$grid_init_locals` in a second
`vcmoe_fit()` call with `control$grid_subset`. This restores the serial
fit's warm start in both directions without the serial dependency.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_grid_neighbor_inits(parts, grid_ids = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   parts:

    A list of `vcmoe_grid_part` objects.

-   grid\_ids:

    Grid indices to build starts for; default all covered ones.

</div>

<div class="section level2">

## Value

A named list (by grid index) of lists of local fits.

</div>

</div>
