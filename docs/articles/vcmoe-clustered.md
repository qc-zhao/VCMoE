<div id="main" class="col-md-9" role="main">

# Clustered Observations: CR2 Bands and Gene Screening

Single-cell data are usually collected from a handful of donors, and
cells from the same donor share unmeasured effects. When `u` is a
donor-level quantity (for example the age or time point of each donor),
all cells at a given `u` come from one or a few donors, so cell-level
standard errors can be far too small. This tutorial goes from a fit to
cluster-robust confidence bands and then to a screen of candidate genes,
treating donors as independent clusters.

<div id="cb1" class="sourceCode">

``` r
library(VCMoE)
```

</div>

<div class="section level2">

## Simulate donor-clustered counts

Twelve donors, each observed at its own time point `u`, contribute 80
cells. Cells belong to one of two latent components; `marker` shifts the
component probability, and the response depends on `signal` with
opposite signs in the two components. A donor effect enters every cell
of that donor.

<div id="cb2" class="sourceCode">

``` r
set.seed(1)
n_donor <- 12L
cells <- 80L
donor <- rep(seq_len(n_donor), each = cells)
u <- ((seq_len(n_donor) - 1) / (n_donor - 1))[donor]
n <- length(donor)
marker <- rnorm(n)
component <- ifelse(rbinom(n, 1, plogis(-0.3 + 1.8 * marker + 0.8 * u)) == 1, 1L, 2L)
signal <- rnorm(n)
library_size <- exp(rnorm(n, 0, 0.3))
eta <- log(library_size) +
  ifelse(component == 1, 1 + 0.3 * u + 0.4 * signal, 0.2 - 0.5 * u - 0.2 * signal) +
  rnorm(n_donor, 0, 0.1)[donor]
dat <- data.frame(
  y = rnbinom(n, size = 6, mu = exp(eta)),
  log_library = log(library_size),
  u = u, donor = donor, marker = marker, signal = signal,
  noise1 = rnorm(n), noise2 = rnorm(n)
)
```

</div>

</div>

<div class="section level2">

## Fit the model and refine it

The null model for the screen contains the library-size offset in the
experts and the component marker in the gating. Fitting one grid point
per donor time point keeps every grid point supported by data. The theta
ridge stabilises the dispersion when few cells carry a component.

<div id="cb3" class="sourceCode">

``` r
fit <- vcmoe_fit(
  y ~ offset(log_library) | marker,
  data = dat, u = "u", k = 2L, family = "negative-binomial",
  bandwidth = 0.3, u_grid = sort(unique(dat$u)),
  control = list(seed = 3, n_starts = 2, maxit = 300,
                 negbin_theta_ridge = 0.05, negbin_theta_target = 8)
)
fit <- vcmoe_refine(fit, tol = 1e-9, maxit = 2000)
head(vcmoe_diagnostics(fit)[, c("u", "converged", "stop_reason", "objective_delta", "refined")])
#>                      u converged stop_reason objective_delta refined
#> 0           0.00000000      TRUE   tolerance    6.369944e-07    TRUE
#> 0.090909091 0.09090909      TRUE   tolerance    1.194968e-06    TRUE
#> 0.18181818  0.18181818      TRUE   tolerance    1.560205e-06    TRUE
#> 0.27272727  0.27272727      TRUE   tolerance    1.624012e-06    TRUE
#> 0.36363636  0.36363636      TRUE   tolerance    1.599456e-06    TRUE
#> 0.45454545  0.45454545      TRUE   tolerance    1.664505e-06    TRUE
```

</div>

`vcmoe_refine()` continues the local EM with a tight tolerance. The
default stopping rule is relative, so with many cells it can stop before
the local maximum; inference then often finds a sandwich bread that is
not positive definite.

</div>

<div class="section level2">

## Cluster-robust bands and path tests

Passing the donor of every fitted row gives CR2 intervals with
Bell-McCaffrey degrees of freedom and simultaneous bands from a wild
cluster bootstrap. Path tests ask whether each coefficient function is
zero or constant along `u`; `path_sets` adds a band over the interior
grid points.

<div id="cb4" class="sourceCode">

``` r
band <- vcmoe_confband(
  fit, cluster = dat$donor[fit$rows_used], type = "simultaneous",
  coefficient_set = c("expert", "gating"),
  control = list(multiplier_B = 999L, seed = 1, path_tests = TRUE,
                 path_sets = list(interior = 2:11))
)
band
#> VCMoE confidence bands
#>   family: negative-binomial
#>   components: 2
#>   engine: local_grid_em
#>   type: simultaneous
#>   covariance: CR2
#>   level: 0.95
#>   coverage theory: local_likelihood_asymptotic_with_bias_and_boundary_caveats
#>   interval rows ok: 96/96
head(band$intervals[, c("u", "parameter", "estimate", "se", "df", "lower", "upper", "status")])
#>   u                                                      parameter   estimate
#> 1 0                        expert:component1:(Intercept):intercept  0.1918718
#> 2 0                        expert:component2:(Intercept):intercept  1.2833384
#> 3 0                            expert:component1:(Intercept):slope  0.5511223
#> 4 0                            expert:component2:(Intercept):slope  0.9268440
#> 5 0 gating_contrast:component1_vs_component2:(Intercept):intercept  0.9697666
#> 6 0      gating_contrast:component1_vs_component2:marker:intercept -1.9174393
#>           se       df      lower      upper status
#> 1 0.01369189 1.135170  0.1208365  0.2629072     ok
#> 2 0.03736900 1.066773  1.1057436  1.4609333     ok
#> 3 0.03579823 1.806939  0.3294041  0.7728404     ok
#> 4 0.22534548 1.758167 -0.1679176  2.0216057     ok
#> 5 0.11294509 1.148191  0.5098455  1.4296876     ok
#> 6 0.25524554 1.182445 -3.2485912 -0.5862874     ok
band$tests[, c("parameter", "path_set", "zero_p_value", "constancy_p_value")]
#>                                                         parameter path_set
#> 1                         expert:component1:(Intercept):intercept     full
#> 2                         expert:component2:(Intercept):intercept     full
#> 3                             expert:component1:(Intercept):slope     full
#> 4                             expert:component2:(Intercept):slope     full
#> 5  gating_contrast:component1_vs_component2:(Intercept):intercept     full
#> 6       gating_contrast:component1_vs_component2:marker:intercept     full
#> 7      gating_contrast:component1_vs_component2:(Intercept):slope     full
#> 8           gating_contrast:component1_vs_component2:marker:slope     full
#> 11                        expert:component1:(Intercept):intercept interior
#> 12                        expert:component2:(Intercept):intercept interior
#> 13                            expert:component1:(Intercept):slope interior
#> 14                            expert:component2:(Intercept):slope interior
#> 15 gating_contrast:component1_vs_component2:(Intercept):intercept interior
#> 16      gating_contrast:component1_vs_component2:marker:intercept interior
#> 17     gating_contrast:component1_vs_component2:(Intercept):slope interior
#> 18          gating_contrast:component1_vs_component2:marker:slope interior
#>    zero_p_value constancy_p_value
#> 1         0.001             0.159
#> 2         0.001             0.259
#> 3         0.001             0.001
#> 4         0.213             0.285
#> 5         0.001             0.001
#> 6         0.001             0.382
#> 7         0.098             0.028
#> 8         0.212             0.026
#> 11        0.019             0.147
#> 12        0.001             0.785
#> 13        0.519             0.615
#> 14        0.260             0.188
#> 15        0.001             0.005
#> 16        0.001             0.641
#> 17        0.215             0.569
#> 18        0.212             0.062
```

</div>

With only twelve donors the function warns that the intervals are
sensitivity analyses; `band$diagnostics` reports the active and
effective number of donors at every grid point.

</div>

<div class="section level2">

## Screen candidate genes

`vcmoe_gene_screen()` asks, for each candidate, whether it is associated
with the response within the components, given the fitted model.
Candidates enter the gating as nuisance, so a gene that only marks the
components is not flagged. A candidate that strongly marks the
components can be listed in `profile_genes`, which refits its gating
effect before testing.

<div id="cb5" class="sourceCode">

``` r
X <- as.matrix(dat[fit$rows_used, c("signal", "noise1", "noise2")])
screen <- vcmoe_gene_screen(fit, X, cluster = dat$donor[fit$rows_used],
                            grids = 2:11, B = 999L, seed = 2)
screen
#> VCMoE gene screen
#>   candidates: 3 (profiled: 0)
#>   grid points: 10; clusters: 12; bootstrap draws: 999
#>   status ok: 3/3
#>   BH q <= 0.05: 1
screen$genes[, c("gene", "statistic", "p_value", "q_value", "argmax_component", "effect_per_sd", "status")]
#>     gene statistic p_value q_value argmax_component effect_per_sd status
#> 1 signal  6.339454   0.001  0.0030       component2    0.15932999     ok
#> 2 noise1  2.772646   0.271  0.4065       component2    0.10768978     ok
#> 3 noise2  2.124346   0.745  0.7450       component2    0.08123199     ok
```

</div>

`signal` is the only candidate with a within-component effect, and it is
the only one with a small adjusted p-value. The screen tests
associations on the log-mean scale of each expert; `effect_per_sd` gives
the size of the strongest one per standard deviation of the candidate.

</div>

</div>
