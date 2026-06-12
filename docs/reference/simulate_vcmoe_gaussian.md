<div id="main" class="col-md-9" role="main">

# Simulate Gaussian VCMoE data

<div class="ref-description section level2">

Generates a small Gaussian no-offset VCMoE simulation for tutorials and
tests.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
simulate_vcmoe_gaussian(n = 300L, k = 2L, seed = NULL,
  separation = 1, u = NULL, scenario = "well_separated")
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

</div>

<div class="section level2">

## Value

A list with `data` and `truth`. The `truth` entry includes component
coefficients, gating logits, probabilities, means, standard deviations,
and sampled class labels.

</div>

</div>
