test_that("vcmoe_diagnostics returns grid-level Gaussian diagnostics", {
  sim <- simulate_vcmoe_gaussian(n = 70, k = 2, seed = 201, separation = 1.2)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    bandwidth = 0.35,
    u_grid = seq(0.2, 0.8, length.out = 3),
    control = list(maxit = 12, n_starts = 1, seed = 202, warn_ambiguous = FALSE)
  ))

  diagnostics <- vcmoe_diagnostics(fit)
  expect_equal(nrow(diagnostics), 3)
  expect_true(all(c(
    "u",
    "converged",
    "posterior_entropy",
    "ambiguous",
    "effective_n",
    "posterior_component1",
    "posterior_component2"
  ) %in% names(diagnostics)))
  expect_true(all(is.finite(diagnostics$effective_n)))
  expect_true(all(diagnostics$effective_n > 0))
  expect_equal(length(fit$diagnostics$ambiguous), 3)

  weights <- VCMoE:::.local_weights(
    fit$fitted$u,
    fit$u_grid[[1L]],
    fit$bandwidth,
    VCMoE:::.vcmoe_estimation_spec(fit$parameterization_id)
  )
  posterior_grid <- fit$posterior[, , 1L]
  manual_component_mean <- colSums(posterior_grid * weights) / sum(weights)
  expect_equal(
    as.numeric(diagnostics[1L, c("posterior_component1", "posterior_component2")]),
    as.numeric(manual_component_mean),
    tolerance = 1e-8
  )
  expect_equal(
    diagnostics$posterior_entropy[[1L]],
    VCMoE:::.weighted_entropy(posterior_grid, weights),
    tolerance = 1e-8
  )
})

test_that("vcmoe_diagnostics includes Binomial optimizer summaries", {
  sim <- simulate_vcmoe_binomial(n = 70, k = 2, seed = 211, separation = 1.2)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    family = "binomial",
    bandwidth = 0.40,
    u_grid = seq(0.2, 0.8, length.out = 3),
    control = list(maxit = 10, n_starts = 1, seed = 212, warn_ambiguous = FALSE)
  ))

  diagnostics <- vcmoe_diagnostics(fit)
  expect_true(all(c(
    "expert_optimizer_nonzero_count",
    "expert_optimizer_na_count",
    "expert_gradient_norm_max",
    "expert_coef_norm_max"
  ) %in% names(diagnostics)))
  expect_true(any(is.finite(diagnostics$expert_gradient_norm_max)))
  expect_true(all(diagnostics$min_posterior_mean >= 0))
  expect_true(all(diagnostics$max_posterior_mean <= 1))
})

test_that("plot_diagnostics returns a ggplot object", {
  sim <- simulate_vcmoe_gaussian(n = 50, k = 2, seed = 221, separation = 1.1)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    bandwidth = 0.40,
    u_grid = seq(0.3, 0.7, length.out = 2),
    control = list(maxit = 8, n_starts = 1, seed = 222, warn_ambiguous = FALSE)
  ))

  expect_s3_class(plot_diagnostics(fit), "ggplot")
})

test_that("high-k diagnostics and plots return stable public shapes", {
  sim <- simulate_vcmoe_gaussian(n = 120, k = 10, seed = 225, separation = 1.7)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ 1 | 1,
    data = sim$data,
    u = "u",
    k = 10,
    bandwidth = 0.55,
    u_grid = c(0.5),
    control = list(maxit = 3, n_starts = 1, seed = 226, warn_ambiguous = FALSE)
  ))

  diagnostics <- vcmoe_diagnostics(fit)
  posterior_cols <- paste0("posterior_component", seq_len(10))
  expect_true(all(posterior_cols %in% names(diagnostics)))
  expect_equal(nrow(diagnostics), 1)
  expect_s3_class(plot_coefficients(fit, type = "expert"), "ggplot")
  expect_s3_class(plot_posterior(fit), "ggplot")
  expect_s3_class(plot_diagnostics(fit), "ggplot")
})

test_that("diagnostics record ambiguity and support fits without retained data", {
  sim <- simulate_vcmoe_gaussian(n = 55, k = 2, seed = 231, separation = 0.03)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    bandwidth = 0.45,
    u_grid = seq(0.25, 0.75, length.out = 2),
    control = list(
      maxit = 8,
      n_starts = 1,
      seed = 232,
      min_separation = 100,
      warn_ambiguous = FALSE,
      keep_data = FALSE
    )
  ))

  diagnostics <- vcmoe_diagnostics(fit)
  expect_true(any(diagnostics$ambiguous))
  expect_true(all(is.na(diagnostics$effective_n)))
})
