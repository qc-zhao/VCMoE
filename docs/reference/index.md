<div id="main" class="col-md-9" role="main">

# Package index

<div class="section level2">

## Fit and inspect models

</div>

<div class="section level2">

-   `vcmoe_fit()` : Fit a varying-coefficient mixture-of-experts model
-   `vcmoe_refine()` : Refine a local-grid VCMoE fit
-   `vcmoe_merge_grid_fits()` : Merge grid-subset VCMoE fits
-   `vcmoe_grid_neighbor_inits()` : Neighbour warm starts for a second
    grid-parallel pass
-   `vcmoe_fit_reduced()` : Fit a block-constant reduced VCMoE model
-   `coef(<vcmoe>)` : Extract VCMoE coefficients
-   `predict(<vcmoe>)` : Predict from a VCMoE fit
-   `vcmoe_diagnostics()` : Summarize VCMoE fit diagnostics
-   `vcmoe_parameterization()` : Inspect VCMoE Parameterization Metadata
-   `vcmoe_gating_contrasts()` : Report Identifiable VCMoE Gating
    Contrasts
-   `vcmoe_scaled_slopes()` : Inspect VCMoE Local-Linear Slopes On The
    Scaled Basis

</div>

<div class="section level2">

## Simulation

</div>

<div class="section level2">

-   `simulate_vcmoe_gaussian()` : Simulate Gaussian VCMoE data
-   `simulate_vcmoe_binomial()` : Simulate Binomial VCMoE data
-   `simulate_vcmoe_negbin()` : Simulate Negative-Binomial VCMoE count
    data

</div>

<div class="section level2">

## Tuning and inference

</div>

<div class="section level2">

-   `vcmoe_select_bandwidth()` : Select a VCMoE bandwidth by K-fold
    cross-validation
-   `vcmoe_bootstrap()` : Parametric bootstrap inference for a VCMoE fit
-   `confint(<vcmoe_bootstrap>)` : Bootstrap confidence intervals for
    VCMoE coefficients
-   `plot_inference()` : Plot bootstrap inference intervals
-   `vcmoe_confband()` : Confidence bands for a VCMoE fit
-   `vcmoe_glrt()` : Generalized likelihood-ratio test for VCMoE
    coefficient variation
-   `vcmoe_gene_screen()` : Cluster-robust gene screen for
    Negative-Binomial expert effects

</div>

<div class="section level2">

## Plots

</div>

<div class="section level2">

-   `plot_coefficients()` : Plot fitted coefficient functions
-   `plot_posterior()` : Plot fitted posterior summaries
-   `plot_diagnostics()` : Plot VCMoE fit diagnostics

</div>

</div>
