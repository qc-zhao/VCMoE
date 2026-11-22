grid_test_fit <- function(dat, formula, k, family, control, u_grid, bandwidth = 0.25) {
  suppressWarnings(vcmoe_fit(formula, data = dat, u = "u", k = k, family = family,
                             bandwidth = bandwidth, u_grid = u_grid, control = control))
}

grid_fit_core <- function(fit) {
  diagnostics <- fit$diagnostics[setdiff(names(fit$diagnostics), c("grid_source", "grid_candidates", "grid_seed"))]
  list(coefficients = fit$coefficients, posterior = fit$posterior, diagnostics = diagnostics,
       control = fit$control, fitted = fit$fitted, u_grid = fit$u_grid)
}

test_that("merging per-grid-seed parts reproduces the serial fit exactly", {
  sim <- simulate_vcmoe_negbin(n = 90L, k = 3L, seed = 21L)
  grid <- seq(0.1, 0.9, length.out = 5)
  control <- list(seed = 4, n_starts = 2, maxit = 3, warn_ambiguous = FALSE,
                  grid_seed = "per_grid", grid_warm_start = FALSE)
  form <- y ~ z1 + offset(log_size_factor) | x1
  serial <- grid_test_fit(sim$data, form, 3L, "negative-binomial", control, grid)
  parts <- lapply(list(c(1L, 4L), c(2L, 5L), 3L), function(ids) {
    grid_test_fit(sim$data, form, 3L, "negative-binomial", c(control, list(grid_subset = ids)), grid)
  })
  expect_s3_class(parts[[1L]], "vcmoe_grid_part")
  expect_output(print(parts[[1L]]), "grid points: 1, 4 of 5")
  merged <- vcmoe_merge_grid_fits(parts)
  expect_s3_class(merged, "vcmoe")
  expect_identical(grid_fit_core(merged), grid_fit_core(serial))
  expect_identical(merged$diagnostics$grid_source, c(1L, 2L, 3L, 1L, 2L))
  expect_null(merged$control$grid_subset)
})

test_that("one full-grid part with default settings merges to the default serial fit", {
  sim <- simulate_vcmoe_gaussian(n = 250L, k = 2L, seed = 22L)
  grid <- seq(0.1, 0.9, length.out = 4)
  control <- list(seed = 5, n_starts = 2, maxit = 3, warn_ambiguous = FALSE)
  serial <- grid_test_fit(sim$data, y ~ z1 | x1, 2L, "gaussian", control, grid)
  part <- grid_test_fit(sim$data, y ~ z1 | x1, 2L, "gaussian", c(control, list(grid_subset = 1:4)), grid)
  expect_identical(grid_fit_core(vcmoe_merge_grid_fits(part)), grid_fit_core(serial))
})

test_that("merge reproduces sequential alignment for k > 6", {
  sim <- simulate_vcmoe_gaussian(n = 100L, k = 7L, seed = 23L)
  grid <- seq(0.2, 0.8, length.out = 3)
  control <- list(seed = 6, n_starts = 1, maxit = 2, warn_ambiguous = FALSE,
                  grid_seed = "per_grid", grid_warm_start = FALSE)
  serial <- grid_test_fit(sim$data, y ~ z1 | x1, 7L, "gaussian", control, grid, bandwidth = 0.3)
  parts <- lapply(1:3, function(id) {
    grid_test_fit(sim$data, y ~ z1 | x1, 7L, "gaussian", c(control, list(grid_subset = id)), grid, bandwidth = 0.3)
  })
  merged <- vcmoe_merge_grid_fits(parts)
  expect_identical(merged$label, "sequential")
  expect_identical(grid_fit_core(merged), grid_fit_core(serial))
})

test_that("a neighbour warm-start pass never lowers the penalized local objective", {
  skip_on_cran()
  sim <- simulate_vcmoe_negbin(n = 300L, k = 2L, seed = 24L)
  grid <- seq(0.1, 0.9, length.out = 5)
  form <- y ~ z1 + offset(log_size_factor) | x1
  control <- list(seed = 7, n_starts = 1, maxit = 500, warn_ambiguous = FALSE,
                  grid_seed = "per_grid", grid_warm_start = FALSE)
  pass1 <- lapply(1:5, function(id) {
    grid_test_fit(sim$data, form, 2L, "negative-binomial", c(control, list(grid_subset = id)), grid)
  })
  inits <- vcmoe_grid_neighbor_inits(pass1)
  expect_named(inits, as.character(1:5))
  expect_length(inits[["1"]], 1L)
  expect_length(inits[["3"]], 2L)
  pass2 <- lapply(1:5, function(id) {
    grid_test_fit(sim$data, form, 2L, "negative-binomial",
                  c(control, list(grid_subset = id, grid_init_locals = inits[as.character(id)])), grid)
  })
  expect_true(pass2[[3L]]$second_pass)
  first <- vcmoe_merge_grid_fits(pass1)
  both <- vcmoe_merge_grid_fits(c(pass1, pass2))
  expect_true(all(first$diagnostics$converged))
  expect_true(all(both$diagnostics$objective >= first$diagnostics$objective - 1e-8))
  expect_identical(both$diagnostics$grid_candidates, rep(2L, 5L))
})

test_that("grid options and merges are validated", {
  sim <- simulate_vcmoe_gaussian(n = 120L, k = 2L, seed = 25L)
  grid <- c(0.3, 0.7)
  control <- list(seed = 8, n_starts = 1, maxit = 5, warn_ambiguous = FALSE)
  expect_error(grid_test_fit(sim$data, y ~ z1 | x1, 2L, "gaussian", c(control, list(grid_subset = 3L)), grid),
               "grid_subset")
  expect_error(grid_test_fit(sim$data, y ~ z1 | x1, 2L, "gaussian",
                             list(grid_seed = "per_grid", n_starts = 1, maxit = 5), grid),
               "requires `control\\$seed`")
  expect_error(grid_test_fit(sim$data, y ~ z1 | x1, 2L, "gaussian",
                             c(control, list(grid_init_locals = list(`1` = list(list())))), grid),
               "requires `grid_subset`")
  p1 <- grid_test_fit(sim$data, y ~ z1 | x1, 2L, "gaussian", c(control, list(grid_subset = 1L)), grid)
  expect_error(vcmoe_merge_grid_fits(list(p1)), "not covered by any part: 2")
  other <- sim$data
  other$y <- other$y + 1
  p2 <- grid_test_fit(other, y ~ z1 | x1, 2L, "gaussian", c(control, list(grid_subset = 2L)), grid)
  expect_error(vcmoe_merge_grid_fits(list(p1, p2)), "same data")
  expect_error(vcmoe_merge_grid_fits(list(1)), "vcmoe_grid_part")
})

test_that("merging can switch to greedy alignment", {
  skip_on_cran()
  sim <- simulate_vcmoe_negbin(n = 300L, k = 2L, seed = 30003L, separation = 1.5)
  grid <- seq(0.1, 0.9, length.out = 7)
  form <- y ~ z1 + offset(log_size_factor) | x1
  control <- list(seed = 31000, n_starts = 3, maxit = 200, warn_ambiguous = FALSE,
                  grid_seed = "per_grid", grid_warm_start = FALSE)
  parts <- lapply(list(1:4, 5:7), function(ids) {
    grid_test_fit(sim$data, form, 2L, "negative-binomial", c(control, list(grid_subset = ids)), grid, bandwidth = 0.3)
  })
  greedy_serial <- suppressWarnings(vcmoe_fit(form, data = sim$data, u = "u", family = "negative-binomial",
                                              bandwidth = 0.3, u_grid = grid, control = control, label = "greedy"))
  merged <- vcmoe_merge_grid_fits(parts, label = "greedy")
  expect_identical(merged$label, "greedy")
  expect_identical(merged$coefficients, greedy_serial$coefficients)
})
