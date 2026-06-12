test_that("simulate_vcmoe_negbin returns count data and truth", {
  sim <- simulate_vcmoe_negbin(n = 40, k = 3, seed = 401, separation = 1.2)
  expect_true(all(c("y", "u", "z1", "x1", "size_factor", "log_size_factor", "component") %in%
    names(sim$data)))
  expect_true(all(sim$data$y >= 0))
  expect_true(all(abs(sim$data$y - round(sim$data$y)) < 1e-8))
  expect_equal(dim(sim$truth$expert), c(40, 3, 2))
  expect_equal(dim(sim$truth$theta), c(40, 3))
  expect_equal(dim(sim$truth$mean), c(40, 3))
})

test_that("invalid Negative-Binomial responses fail clearly", {
  sim <- simulate_vcmoe_negbin(n = 25, k = 2, seed = 402)$data
  bad <- sim
  bad$y[1] <- -1
  expect_error(
    vcmoe_fit(y ~ z1 + offset(log_size_factor) | x1, bad, u = "u", family = "negative-binomial"),
    "non-negative"
  )

  bad <- sim
  bad$y[1] <- 1.5
  expect_error(
    vcmoe_fit(y ~ z1 + offset(log_size_factor) | x1, bad, u = "u", family = "negative-binomial"),
    "whole numbers"
  )

  bad <- sim
  bad$y[1] <- Inf
  expect_error(
    vcmoe_fit(y ~ z1 + offset(log_size_factor) | x1, bad, u = "u", family = "negative-binomial"),
    "whole numbers"
  )
})

test_that("Negative-Binomial offset affects predictions and is not a coefficient term", {
  sim <- simulate_vcmoe_negbin(n = 55, k = 2, seed = 403, separation = 1.2)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 + offset(log_size_factor) | x1,
    data = sim$data,
    u = "u",
    family = "negative-binomial",
    bandwidth = 0.42,
    u_grid = seq(0.3, 0.7, length.out = 2),
    control = list(maxit = 6, n_starts = 1, seed = 404, warn_ambiguous = FALSE)
  ))

  expect_false(any(grepl("offset", dimnames(coef(fit, "expert"))[[3L]], fixed = TRUE)))
  expect_equal(dim(coef(fit, "theta")), c(2, 2))
  expect_true(all(is.finite(coef(fit, "theta"))))
  expect_true(all(coef(fit, "theta") > 0))

  low <- sim$data[1:6, ]
  high <- low
  high$log_size_factor <- low$log_size_factor + log(2)
  high$size_factor <- low$size_factor * 2
  pred_low <- predict(fit, newdata = low, type = "mean")
  pred_high <- predict(fit, newdata = high, type = "mean")
  expect_true(mean(pred_high / pred_low) > 1.4)
})

test_that("Negative-Binomial fit and predictions have stable shapes", {
  sim <- simulate_vcmoe_negbin(n = 60, k = 2, seed = 405, separation = 1.25)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 + offset(log_size_factor) | x1,
    data = sim$data,
    u = "u",
    family = "negative-binomial",
    bandwidth = 0.40,
    u_grid = seq(0.25, 0.75, length.out = 3),
    control = list(maxit = 8, n_starts = 1, seed = 406, warn_ambiguous = FALSE)
  ))

  expect_equal(fit$family, "negative-binomial")
  expect_equal(dim(coef(fit, "expert")), c(3, 2, 2))
  expect_equal(dim(coef(fit, "gating")), c(3, 2, 2))
  expect_equal(dim(coef(fit, "theta")), c(3, 2))
  expect_true(all(is.finite(predict(fit, type = "mean"))))
  expect_true(all(predict(fit, type = "mean") >= 0))
  expect_equal(dim(predict(fit, type = "component")), c(nrow(sim$data), 2))
  expect_true(all(predict(fit, type = "component") >= 0))
  posterior <- predict(fit, type = "posterior")
  expect_equal(dim(posterior), c(nrow(sim$data), 2))
  expect_equal(rowSums(posterior), rep(1, nrow(sim$data)), tolerance = 1e-8)

  diagnostics <- vcmoe_diagnostics(fit)
  expect_true(all(c(
    "theta_min", "theta_max", "theta_boundary_count", "theta_boundary",
    "negbin_eta_min_observed", "negbin_eta_max_observed",
    "negbin_eta_clipping_count", "min_component_effective_n"
  ) %in% names(diagnostics)))
  expect_true(all(is.finite(diagnostics$min_component_effective_n)))
  expect_true(all(diagnostics$negbin_eta_clipping_count >= 0))
})

test_that("Negative-Binomial theta ridge can prevent Poisson-boundary estimates", {
  y <- rep(5, 25)
  mu <- rep(5, 25)
  weights <- rep(1, 25)
  unpenalized <- VCMoE:::.fit_negbin_theta(
    y,
    mu,
    weights,
    previous_theta = 8,
    control = VCMoE:::.vcmoe_default_control(list(negbin_theta_max = 1e4))
  )
  penalized <- VCMoE:::.fit_negbin_theta(
    y,
    mu,
    weights,
    previous_theta = 8,
    control = VCMoE:::.vcmoe_default_control(list(
      negbin_theta_max = 1e4,
      negbin_theta_ridge = 5,
      negbin_theta_target = 8
    ))
  )

  expect_gt(unpenalized, penalized)
  expect_lt(penalized, 100)
  expect_gt(penalized, 0)
})

test_that("Negative-Binomial bandwidth selection and bootstrap support k = 2", {
  sim <- simulate_vcmoe_negbin(n = 45, k = 2, seed = 407, separation = 1.2)
  selection <- suppressWarnings(vcmoe_select_bandwidth(
    y ~ z1 + offset(log_size_factor) | x1,
    data = sim$data,
    u = "u",
    family = "negative-binomial",
    bandwidth_grid = c(0.35),
    folds = 2,
    u_grid = seq(0.35, 0.65, length.out = 2),
    control = list(maxit = 4, n_starts = 1, seed = 408, warn_ambiguous = FALSE),
    seed = 409
  ))
  expect_s3_class(selection, "vcmoe_bandwidth_selection")
  expect_equal(selection$best_bandwidth, 0.35)
  expect_true(is.finite(selection$cv_summary$total_loglik))

  boot <- suppressWarnings(vcmoe_bootstrap(
    selection$fit,
    sim$data,
    B = 2,
    seed = 410,
    min_successful = 2,
    coefficient_set = "expert",
    control = list(maxit = 3, n_starts = 1, warn_ambiguous = FALSE)
  ))
  expect_s3_class(boot, "vcmoe_bootstrap")
  expect_equal(dim(boot$replicates$expert)[1:3], dim(coef(selection$fit, "expert")))
  expect_equal(sum(boot$replicate_summary$status == "ok"), 2)
})

test_that("high-k assignment alignment supports k = 8 and global falls back to sequential", {
  costs <- matrix(10, nrow = 8, ncol = 8)
  target <- c(3, 1, 8, 2, 7, 6, 4, 5)
  costs[cbind(target, seq_len(8))] <- seq(0.01, 0.08, length.out = 8)
  assignment <- VCMoE:::.assignment_permutation(costs)
  expect_equal(assignment$permutation, target)
  expect_true(is.finite(assignment$second_score))
  expect_gt(assignment$second_score, assignment$score)

  tied_assignment <- VCMoE:::.assignment_permutation(matrix(0, nrow = 8, ncol = 8))
  expect_equal(tied_assignment$second_score - tied_assignment$score, 0)

  sim <- simulate_vcmoe_gaussian(n = 80, k = 8, seed = 411, separation = 1.3)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    k = 8,
    bandwidth = 0.50,
    u_grid = seq(0.35, 0.65, length.out = 2),
    control = list(maxit = 3, n_starts = 1, seed = 412, warn_ambiguous = FALSE)
  ))
  expect_equal(fit$diagnostics$alignment_method, "sequential")
  expect_equal(dim(coef(fit, "expert")), c(2, 8, 2))
  global_fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    k = 8,
    bandwidth = 0.50,
    u_grid = seq(0.35, 0.65, length.out = 2),
    label = "global",
    control = list(maxit = 2, n_starts = 1, seed = 413, warn_ambiguous = FALSE)
  ))
  expect_equal(global_fit$diagnostics$alignment_method, "sequential")
  expect_equal(dim(coef(global_fit, "expert")), c(2, 8, 2))
})
