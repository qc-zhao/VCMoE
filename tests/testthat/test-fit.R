test_that("vcmoe_fit returns stable coefficient and prediction shapes", {
  sim <- simulate_vcmoe_gaussian(n = 90, k = 2, seed = 2, separation = 1.2)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    bandwidth = 0.35,
    u_grid = seq(0.1, 0.9, length.out = 5),
    control = list(maxit = 25, n_starts = 1, seed = 3, warn_ambiguous = FALSE)
  ))

  expect_s3_class(fit, "vcmoe")
  expect_equal(fit$label, "global")
  expect_equal(fit$diagnostics$alignment_method, "global")
  expect_equal(dim(coef(fit, "expert")), c(5, 2, 2))
  expect_equal(dim(coef(fit, "gating")), c(5, 2, 2))
  expect_equal(dim(coef(fit, "sigma")), c(5, 2))
  expect_equal(dim(coef(fit, "sigma_slope")), c(5, 2))
  expect_true(all(is.finite(coef(fit, "sigma_slope"))))

  fitted_mean <- predict(fit, type = "mean")
  posterior <- predict(fit, type = "posterior")
  component_mean <- predict(fit, type = "component")

  expect_length(fitted_mean, nrow(sim$data))
  expect_equal(dim(posterior), c(nrow(sim$data), 2))
  expect_equal(dim(component_mean), c(nrow(sim$data), 2))
  expect_equal(as.numeric(rowSums(posterior)), rep(1, nrow(sim$data)), tolerance = 1e-6)
  expect_true(any(fit$diagnostics$converged))
  expect_equal(dim(fit$diagnostics$global_permutations), c(5, 2))
})

test_that("greedy label strategy remains available for comparison", {
  sim <- simulate_vcmoe_gaussian(n = 70, k = 2, seed = 12, separation = 1.1)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    label = "greedy",
    bandwidth = 0.40,
    u_grid = seq(0.2, 0.8, length.out = 3),
    control = list(maxit = 15, n_starts = 1, seed = 13, warn_ambiguous = FALSE)
  ))

  expect_equal(fit$label, "greedy")
  expect_equal(fit$diagnostics$alignment_method, "greedy")
  expect_equal(dim(fit$diagnostics$global_permutations), c(3, 2))
  expect_true(all(is.na(fit$diagnostics$global_transition_cost)))
})

test_that("near-overlap components are flagged as label ambiguous", {
  sim <- simulate_vcmoe_gaussian(n = 70, k = 2, seed = 4, separation = 0.02)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    bandwidth = 0.45,
    u_grid = seq(0.2, 0.8, length.out = 3),
    control = list(
      maxit = 15,
      n_starts = 1,
      seed = 5,
      min_separation = 100,
      warn_ambiguous = FALSE
    )
  ))

  expect_gt(length(fit$diagnostics$warnings), 0)
})

test_that("missing rows are removed consistently across response, covariates, and u", {
  sim <- simulate_vcmoe_gaussian(n = 60, k = 2, seed = 6, separation = 1)
  sim$data$z1[3] <- NA
  sim$data$x1[5] <- NA
  sim$data$u[7] <- NA

  expect_warning(
    fit <- vcmoe_fit(
      y ~ z1 | x1,
      data = sim$data,
      u = "u",
      bandwidth = 0.40,
      u_grid = seq(0.2, 0.8, length.out = 3),
      control = list(maxit = 15, n_starts = 1, seed = 7, warn_ambiguous = FALSE)
    ),
    "Removed 3 row"
  )

  expect_equal(fit$rows_removed, 3)
  expect_equal(length(predict(fit, type = "mean")), 57)
  expect_equal(dim(predict(fit, type = "posterior")), c(57, 2))
})

test_that("k = 3 and k = 5 smoke fits return valid diagnostics", {
  for (k in c(3L, 5L)) {
    sim <- simulate_vcmoe_gaussian(
      n = if (k == 3L) 90 else 120,
      k = k,
      seed = 20 + k,
      separation = 1.4,
      scenario = "well_separated"
    )
    fit <- suppressWarnings(vcmoe_fit(
      y ~ z1 | x1,
      data = sim$data,
      u = "u",
      k = k,
      bandwidth = 0.45,
      u_grid = seq(0.25, 0.75, length.out = 2),
      control = list(maxit = 12, n_starts = 1, seed = 30 + k, warn_ambiguous = FALSE)
    ))

    expect_equal(dim(coef(fit, "expert")), c(2, k, 2))
    expect_equal(dim(coef(fit, "gating")), c(2, k, 2))
    expect_equal(dim(coef(fit, "sigma")), c(2, k))
    expect_equal(dim(coef(fit, "sigma_slope")), c(2, k))
    expect_equal(dim(fit$diagnostics$permutations), c(2, k))
    expect_length(fit$diagnostics$loglik, 2)
    expect_length(fit$diagnostics$posterior_entropy, 2)

    posterior <- predict(fit, type = "posterior")
    expect_equal(dim(posterior), c(nrow(sim$data), k))
    expect_equal(as.numeric(rowSums(posterior)), rep(1, nrow(sim$data)), tolerance = 1e-6)
  }
})

test_that("k = 10 public fit path returns stable shapes and sequential alignment", {
  sim <- simulate_vcmoe_gaussian(
    n = 130,
    k = 10,
    seed = 36,
    separation = 1.7,
    scenario = "well_separated"
  )
  fit <- suppressWarnings(vcmoe_fit(
    y ~ 1 | 1,
    data = sim$data,
    u = "u",
    k = 10,
    label = "global",
    bandwidth = 0.55,
    u_grid = c(0.5),
    control = list(maxit = 3, n_starts = 1, seed = 37, warn_ambiguous = FALSE)
  ))

  expect_s3_class(fit, "vcmoe")
  expect_equal(fit$k, 10)
  expect_equal(fit$diagnostics$alignment_method, "sequential")
  expect_equal(dim(coef(fit, "expert")), c(1, 10, 1))
  expect_equal(dim(coef(fit, "gating")), c(1, 10, 1))
  expect_equal(dim(coef(fit, "sigma")), c(1, 10))
  expect_equal(dim(fit$diagnostics$permutations), c(1, 10))
  expect_equal(dim(predict(fit, type = "posterior")), c(nrow(sim$data), 10))
  expect_equal(dim(predict(fit, type = "component")), c(nrow(sim$data), 10))
  expect_length(predict(fit, type = "mean"), nrow(sim$data))
  expect_equal(as.numeric(rowSums(predict(fit, type = "posterior"))),
               rep(1, nrow(sim$data)), tolerance = 1e-6)
})

test_that("Binomial and Negative-Binomial k = 10 smoke fits return public shapes", {
  grouped <- simulate_vcmoe_binomial(n = 120, k = 10, seed = 38, separation = 1.7, trials = 6)
  grouped_fit <- suppressWarnings(vcmoe_fit(
    cbind(success, failure) ~ 1 | 1,
    data = grouped$data,
    u = "u",
    k = 10,
    family = "binomial",
    bandwidth = 0.55,
    u_grid = c(0.5),
    control = list(maxit = 3, n_starts = 1, seed = 39, warn_ambiguous = FALSE)
  ))
  expect_equal(grouped_fit$k, 10)
  expect_equal(grouped_fit$family, "binomial")
  expect_equal(dim(coef(grouped_fit, "expert")), c(1, 10, 1))
  expect_equal(dim(predict(grouped_fit, type = "posterior")), c(nrow(grouped$data), 10))
  expect_true(all(predict(grouped_fit, type = "component") >= 0 & predict(grouped_fit, type = "component") <= 1))

  nb <- simulate_vcmoe_negbin(n = 120, k = 10, seed = 40, separation = 1.7)
  nb_fit <- suppressWarnings(vcmoe_fit(
    y ~ 1 + offset(log_size_factor) | 1,
    data = nb$data,
    u = "u",
    k = 10,
    family = "negative-binomial",
    bandwidth = 0.55,
    u_grid = c(0.5),
    control = list(
      maxit = 3,
      n_starts = 1,
      seed = 41,
      warn_ambiguous = FALSE,
      ridge = 3,
      negbin_theta_ridge = 0.05,
      negbin_theta_target = 8
    )
  ))
  expect_equal(nb_fit$k, 10)
  expect_equal(nb_fit$family, "negative-binomial")
  expect_equal(dim(coef(nb_fit, "expert")), c(1, 10, 1))
  expect_equal(dim(coef(nb_fit, "theta")), c(1, 10))
  expect_equal(dim(predict(nb_fit, type = "posterior")), c(nrow(nb$data), 10))
  expect_true(all(is.finite(predict(nb_fit, type = "mean"))))
  expect_true(all(predict(nb_fit, type = "mean") >= 0))
})

test_that("binomial Bernoulli fit returns probabilities and stable shapes", {
  sim <- simulate_vcmoe_binomial(n = 90, k = 2, seed = 41, separation = 1.3)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    family = "binomial",
    bandwidth = 0.40,
    u_grid = seq(0.2, 0.8, length.out = 3),
    control = list(maxit = 15, n_starts = 1, seed = 42, warn_ambiguous = FALSE)
  ))

  expect_s3_class(fit, "vcmoe")
  expect_equal(fit$family, "binomial")
  expect_null(coef(fit, "sigma"))
  expect_null(coef(fit, "sigma_slope"))
  expect_equal(dim(coef(fit, "expert")), c(3, 2, 2))
  expect_equal(dim(coef(fit, "gating")), c(3, 2, 2))
  expect_equal(dim(fit$diagnostics$expert_optimizer_convergence), c(3, 2))
  expect_equal(dim(fit$diagnostics$expert_gradient_norm), c(3, 2))
  expect_equal(dim(fit$diagnostics$expert_coef_norm), c(3, 2))
  expect_true(any(is.finite(fit$diagnostics$expert_gradient_norm)))

  marginal <- predict(fit, type = "mean")
  component <- predict(fit, type = "component")
  posterior <- predict(fit, type = "posterior")
  expect_equal(dim(component), c(nrow(sim$data), 2))
  expect_true(all(marginal >= 0 & marginal <= 1))
  expect_true(all(component >= 0 & component <= 1))
  expect_equal(as.numeric(rowSums(posterior)), rep(1, nrow(sim$data)), tolerance = 1e-6)
})

test_that("Bernoulli-only default strengthens gating ridge", {
  bernoulli <- simulate_vcmoe_binomial(n = 45, k = 2, seed = 410, separation = 1.2)
  bernoulli_fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = bernoulli$data,
    u = "u",
    family = "binomial",
    bandwidth = 0.45,
    u_grid = c(0.35, 0.65),
    control = list(maxit = 2, n_starts = 1, seed = 411, warn_ambiguous = FALSE)
  ))
  expect_equal(bernoulli_fit$response_info$type, "bernoulli")
  expect_equal(bernoulli_fit$control$ridge, 1)
  expect_equal(vcmoe_parameterization(bernoulli_fit)$optimization$ridge, 1)

  explicit_fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = bernoulli$data,
    u = "u",
    family = "binomial",
    bandwidth = 0.45,
    u_grid = c(0.35, 0.65),
    control = list(maxit = 2, n_starts = 1, seed = 412, warn_ambiguous = FALSE, ridge = 0.2)
  ))
  expect_equal(explicit_fit$control$ridge, 0.2)

  grouped <- simulate_vcmoe_binomial(n = 45, k = 2, seed = 413, separation = 1.2, trials = 5)
  grouped_fit <- suppressWarnings(vcmoe_fit(
    cbind(success, failure) ~ z1 | x1,
    data = grouped$data,
    u = "u",
    family = "binomial",
    bandwidth = 0.45,
    u_grid = c(0.35, 0.65),
    control = list(maxit = 2, n_starts = 1, seed = 414, warn_ambiguous = FALSE)
  ))
  expect_equal(grouped_fit$response_info$type, "grouped")
  expect_equal(grouped_fit$control$ridge, VCMoE:::.vcmoe_default_control(list())$ridge)

  grouped_one_trial <- bernoulli$data
  grouped_one_trial$success <- grouped_one_trial$y
  grouped_one_trial$failure <- 1 - grouped_one_trial$y
  grouped_one_trial_fit <- suppressWarnings(vcmoe_fit(
    cbind(success, failure) ~ z1 | x1,
    data = grouped_one_trial,
    u = "u",
    family = "binomial",
    bandwidth = 0.45,
    u_grid = c(0.35, 0.65),
    control = list(maxit = 2, n_starts = 1, seed = 415, warn_ambiguous = FALSE)
  ))
  expect_equal(grouped_one_trial_fit$response_info$type, "grouped")
  expect_equal(grouped_one_trial_fit$control$ridge, VCMoE:::.vcmoe_default_control(list())$ridge)
})

test_that("binomial mean prediction is marginal and does not depend on observed response", {
  sim <- simulate_vcmoe_binomial(n = 80, k = 2, seed = 49, separation = 1.2)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    family = "binomial",
    bandwidth = 0.40,
    u_grid = seq(0.2, 0.8, length.out = 3),
    control = list(maxit = 12, n_starts = 1, seed = 50, warn_ambiguous = FALSE)
  ))

  with_response <- predict(fit, type = "mean")
  without_response <- predict(fit, newdata = sim$data[, setdiff(names(sim$data), "y")], type = "mean")
  prior <- predict(fit, type = "prior")
  component <- predict(fit, type = "component")

  expect_equal(with_response, without_response, tolerance = 1e-8)
  expect_equal(with_response, rowSums(prior * component), tolerance = 1e-8)
})

test_that("binomial grouped fit supports cbind response", {
  sim <- simulate_vcmoe_binomial(n = 80, k = 2, seed = 43, separation = 1.2, trials = 6)
  fit <- suppressWarnings(vcmoe_fit(
    cbind(success, failure) ~ z1 | x1,
    data = sim$data,
    u = "u",
    family = "binomial",
    bandwidth = 0.40,
    u_grid = seq(0.2, 0.8, length.out = 3),
    control = list(maxit = 12, n_starts = 1, seed = 44, warn_ambiguous = FALSE)
  ))

  expect_equal(fit$response_info$type, "grouped")
  expect_equal(fit$response_info$success, "success")
  expect_equal(fit$response_info$failure, "failure")
  expect_equal(dim(predict(fit, type = "posterior")), c(nrow(sim$data), 2))
  expect_true(all(predict(fit, type = "mean") >= 0 & predict(fit, type = "mean") <= 1))
})

test_that("invalid binomial responses fail clearly", {
  sim <- simulate_vcmoe_binomial(n = 30, k = 2, seed = 45)
  bad <- sim$data
  bad$y[1] <- 0.5
  expect_error(
    vcmoe_fit(y ~ z1 | x1, bad, u = "u", family = "binomial"),
    "0/1"
  )

  grouped <- simulate_vcmoe_binomial(n = 30, k = 2, seed = 46, trials = 4)$data
  grouped$success[1] <- -1
  expect_error(
    vcmoe_fit(cbind(success, failure) ~ z1 | x1, grouped, u = "u", family = "binomial"),
    "non-negative"
  )

  grouped <- simulate_vcmoe_binomial(n = 30, k = 2, seed = 47, trials = 4)$data
  grouped$success[1] <- Inf
  expect_error(
    vcmoe_fit(cbind(success, failure) ~ z1 | x1, grouped, u = "u", family = "binomial"),
    "finite whole"
  )

  grouped <- simulate_vcmoe_binomial(n = 30, k = 2, seed = 48, trials = 4)$data
  grouped$success[1] <- 0
  grouped$failure[1] <- 0
  expect_error(
    vcmoe_fit(cbind(success, failure) ~ z1 | x1, grouped, u = "u", family = "binomial"),
    "positive trial"
  )
})

test_that("binomial k = 3 and k = 5 smoke fits return valid diagnostics", {
  expect_warning(
    vcmoe_fit(
      y ~ z1 | x1,
      data = simulate_vcmoe_binomial(n = 30, k = 3, seed = 51)$data,
      u = "u",
      k = 3,
      family = "binomial",
      bandwidth = 0.50,
      u_grid = seq(0.3, 0.7, length.out = 2),
      control = list(maxit = 2, n_starts = 1, seed = 52, warn_ambiguous = FALSE)
    ),
    "experimental"
  )

  for (k in c(3L, 5L)) {
    sim <- simulate_vcmoe_binomial(
      n = if (k == 3L) 90 else 110,
      k = k,
      seed = 50 + k,
      separation = 1.4,
      scenario = "well_separated"
    )
    fit <- suppressWarnings(vcmoe_fit(
      y ~ z1 | x1,
      data = sim$data,
      u = "u",
      k = k,
      family = "binomial",
      bandwidth = 0.45,
      u_grid = seq(0.25, 0.75, length.out = 2),
      control = list(maxit = 8, n_starts = 1, seed = 60 + k, warn_ambiguous = FALSE)
    ))

    expect_equal(dim(coef(fit, "expert")), c(2, k, 2))
    expect_equal(dim(coef(fit, "gating")), c(2, k, 2))
    expect_null(coef(fit, "sigma"))
    expect_equal(dim(fit$diagnostics$permutations), c(2, k))
    posterior <- predict(fit, type = "posterior")
    component <- predict(fit, type = "component")
    expect_equal(dim(posterior), c(nrow(sim$data), k))
    expect_equal(dim(component), c(nrow(sim$data), k))
    expect_equal(as.numeric(rowSums(posterior)), rep(1, nrow(sim$data)), tolerance = 1e-6)
    expect_true(all(component >= 0 & component <= 1))
  }
})
