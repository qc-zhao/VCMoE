<div id="main" class="col-md-9" role="main">

# Simulate Binomial VCMoE data

<div class="ref-description section level2">

Generates Binomial VCMoE simulations for Bernoulli and grouped-count
examples and tests.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
simulate_vcmoe_binomial(n = 300L, k = 2L, seed = NULL,
  separation = 1, u = NULL, scenario = "well_separated", trials = 1L)
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

-   trials:

    Binomial trial counts. Use `1` for Bernoulli data, or a positive
    integer scalar/vector for grouped Binomial data.

</div>

<div class="section level2">

## Value

A list with `data` and `truth`. Expert truth is on the logit scale. The
`truth` entry includes component coefficients, gating logits, component
probabilities, component-specific success probabilities, sampled class
labels, and success/failure counts.

</div>

</div>
