<div id="main" class="col-md-9" role="main">

# Simulate Negative-Binomial VCMoE count data

<div class="ref-description section level2">

Generates Negative-Binomial VCMoE simulations for gene-expression count
examples and tests.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
simulate_vcmoe_negbin(n = 300L, k = 2L, seed = NULL,
  separation = 1, u = NULL, scenario = "well_separated",
  size_factor = NULL, mean_count = 5)
```

</div>

</div>

<div class="section level2">

## Arguments

-   n:

    Number of observations.

-   k:

    Number of components. Values 2 through 10 are supported.

-   seed:

    Optional random seed.

-   separation:

    Controls expert separation.

-   u:

    Optional numeric vector of index values.

-   scenario:

    Simulation scenario: `"well_separated"`, `"moderate"`,
    `"near_overlap"`, `"crossing"`, or `"imbalanced_gating"`.

-   size\_factor:

    Optional positive size factors. If `NULL`, log-normal size factors
    are generated.

-   mean\_count:

    Baseline count scale.

</div>

<div class="section level2">

## Value

A list with `data` and `truth`. Expert truth is on the log mean count
scale. The data include `size_factor` and `log_size_factor` for use with
`offset(log_size_factor)`.

</div>

</div>
