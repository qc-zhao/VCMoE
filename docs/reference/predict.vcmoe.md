<div id="main" class="col-md-9" role="main">

# Predict from a VCMoE fit

<div class="ref-description section level2">

Returns fitted means, component-specific means, posterior probabilities,
or gating probabilities.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
# S3 method for class 'vcmoe'
predict(object, newdata = NULL, u = NULL,
  type = c("mean", "posterior", "component", "prior"), ...)
```

</div>

</div>

<div class="section level2">

## Arguments

-   object:

    A `vcmoe` object.

-   newdata:

    Optional data frame.

-   u:

    Optional index values for `newdata`.

-   type:

    Prediction type.

-   ...:

    Unused.

</div>

<div class="section level2">

## Value

A vector or matrix depending on `type`.

</div>

<div class="section level2">

## Details

For Gaussian fits, `type = "component"` returns component-specific means
and `type = "mean"` returns the posterior-weighted fitted mean. For
Binomial fits, `type = "component"` returns component-specific success
probabilities and `type = "mean"` returns the marginal success
probability. For Negative-Binomial fits, `type = "component"` returns
component-specific mean counts and `type = "mean"` returns the marginal
mean count.

</div>

</div>
