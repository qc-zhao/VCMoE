glrt_gaussian_fit <- function(n = 90, k = 2, seed = 3361, u_grid = c(0.2, 0.4, 0.6, 0.8),
                              bandwidth = 0.35, maxit = 15) {
  sim <- simulate_vcmoe_gaussian(n = n, k = k, seed = seed, separation = 1.5)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1, data = sim$data, u = "u", k = k,
    bandwidth = bandwidth, u_grid = u_grid,
    control = list(maxit = maxit, n_starts = 1, seed = seed + 1, warn_ambiguous = FALSE)
  ))
  list(sim = sim, fit = fit)
}

test_that("lambda uses the observation-level likelihood", {
  case <- glrt_gaussian_fit(u_grid = c(0.3, 0.7))
  fit <- case$fit
  args <- list(
    fit = fit, data = case$sim$data, test = "coefficient", coefficient_set = "expert",
    component = 1, term = "z1", calibration = "analytic_epanechnikov",
    control = list(maxit = 3)
  )
  glrt <- suppressWarnings(do.call(vcmoe_glrt, args))
  expect_identical(glrt$lambda_type, "observation")
  expect_gt(glrt$lambda, 0)
  expect_false(isTRUE(all.equal(glrt$lambda, glrt$lambda_fit_objective)))
  expect_equal(
    VCMoE:::.vcmoe_glrt_observation_loglik(fit),
    VCMoE:::.vcmoe_predictive_loglik(fit, newdata = case$sim$data),
    tolerance = 1e-8
  )
  expect_equal(
    glrt$lambda,
    VCMoE:::.vcmoe_glrt_observation_loglik(fit) - VCMoE:::.vcmoe_glrt_observation_loglik(glrt$null_fit),
    tolerance = 1e-8
  )
  expect_equal(glrt$lambda_observation, glrt$lambda)
  expect_equal(
    glrt$lambda_fit_objective,
    VCMoE:::.vcmoe_glrt_global_pseudologlik(fit) - VCMoE:::.vcmoe_glrt_global_pseudologlik(glrt$null_fit),
    tolerance = 1e-8
  )
  expect_equal(glrt$lambda_grid_sum, glrt$lambda_fit_objective, tolerance = 1e-8)
  expect_equal(glrt$analytic_statistic, glrt$rK * max(0, glrt$lambda), tolerance = 1e-8)
  expect_equal(glrt$C_const, 1)
  expect_equal(glrt$delta, glrt$rK * glrt$p_test * (0.75 - 0.5 * 0.6) / 0.35, tolerance = 1e-8)

  args$control$lambda_type <- "fit_objective"
  legacy <- suppressWarnings(do.call(vcmoe_glrt, args))
  expect_identical(legacy$lambda_type, "fit_objective")
  expect_equal(legacy$lambda, glrt$lambda_fit_objective, tolerance = 1e-8)
  expect_equal(legacy$lambda_observation, glrt$lambda, tolerance = 1e-8)
  expect_identical(legacy$settings$theory_status, "analytic_calibration_invalid_for_grid_summed_statistic")

  args$control$lambda_type <- "bogus"
  expect_error(suppressWarnings(do.call(vcmoe_glrt, args)), "should be one of")
})

test_that("the observation-level likelihood matches a hand computation", {
  case <- glrt_gaussian_fit(n = 40, u_grid = c(0.25, 0.75), bandwidth = 0.6, maxit = 5)
  fit <- case$fit
  d <- case$sim$data
  # Interpolate every coefficient path to each observation's own u (sigma on
  # the log scale).
  interp <- function(path) stats::approx(fit$u_grid, path, xout = fit$fitted$u, rule = 2)$y
  k <- fit$k
  loglik_terms <- matrix(NA_real_, nrow(d), k)
  gating_eta <- sapply(seq_len(k), function(c) {
    interp(fit$coefficients$gating[, c, "(Intercept)"]) + d$x1 * interp(fit$coefficients$gating[, c, "x1"])
  })
  prior <- exp(gating_eta - apply(gating_eta, 1L, max))
  prior <- prior / rowSums(prior)
  for (c in seq_len(k)) {
    mean <- interp(fit$coefficients$expert[, c, "(Intercept)"]) + d$z1 * interp(fit$coefficients$expert[, c, "z1"])
    sd <- exp(interp(log(fit$coefficients$sigma[, c])))
    loglik_terms[, c] <- log(prior[, c]) + stats::dnorm(d$y, mean, sd, log = TRUE)
  }
  by_hand <- sum(apply(loglik_terms, 1L, function(x) max(x) + log(sum(exp(x - max(x))))))
  expect_equal(VCMoE:::.vcmoe_glrt_observation_loglik(fit), by_hand, tolerance = 1e-6)
})

test_that("a denser grid does not count observations more than once", {
  skip_on_cran()
  # Constant coefficient paths: the observation-level likelihood does not
  # depend on the grid, while the grid sum grows with the number of grid points.
  case <- glrt_gaussian_fit(n = 60, u_grid = c(0.3, 0.7), bandwidth = 0.5, maxit = 5)
  coarse <- case$fit
  dense <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1, data = case$sim$data, u = "u", bandwidth = 0.5, u_grid = seq(0.3, 0.7, length.out = 9),
    control = list(maxit = 5, n_starts = 1, seed = 3362, warn_ambiguous = FALSE)
  ))
  freeze <- function(fit, template) {
    for (set in c("expert", "gating")) {
      value <- template$coefficients[[set]][1L, , , drop = FALSE]
      fit$coefficients[[set]][] <- value[rep(1L, length(fit$u_grid)), , , drop = FALSE]
      fit$coefficients[[paste0(set, "_slope")]][] <- 0
    }
    fit$coefficients$sigma[] <- matrix(template$coefficients$sigma[1L, ], nrow = length(fit$u_grid),
                                       ncol = fit$k, byrow = TRUE)
    fit$coefficients$sigma_slope[] <- 0
    fit
  }
  a <- freeze(coarse, coarse)
  b <- freeze(dense, coarse)
  expect_equal(VCMoE:::.vcmoe_glrt_observation_loglik(a), VCMoE:::.vcmoe_glrt_observation_loglik(b),
               tolerance = 1e-8)
  ratio <- VCMoE:::.vcmoe_dev_global_pseudologlik(b) / VCMoE:::.vcmoe_dev_global_pseudologlik(a)
  expect_gt(ratio, 3)
})

test_that("analytic degrees of freedom count the functions constrained under the null", {
  default <- VCMoE:::.vcmoe_glrt_rk_delta_epanechnikov(bandwidth = 0.22, p_test = 1)
  expect_equal(default$C_const, 1)
  expect_equal(default$delta, 8.653392 / 2, tolerance = 1e-6)
  expect_equal(VCMoE:::.vcmoe_glrt_control(list())$C_const, 1)
  expect_false(VCMoE:::.vcmoe_glrt_control(list())$strict)
  counts <- lapply(2:3, function(k) {
    fit <- glrt_gaussian_fit(n = 60 * k, k = k, seed = 3370 + k, u_grid = c(0.35, 0.65),
                             bandwidth = 0.45, maxit = 3)$fit
    p_test <- function(...) VCMoE:::.vcmoe_glrt_constraint(fit, ...)$p_test
    c(
      expert = p_test("coefficient", "expert", 1, "z1"),
      gating = p_test("coefficient", "gating", if (k == 2) 1 else 2, "x1"),
      gating_block = p_test("constant_block", "gating", NULL, NULL),
      expert_block = p_test("constant_block", "expert", NULL, NULL),
      all = p_test("constant_all", "expert", NULL, NULL)
    )
  })
  expect_equal(counts[[1L]], c(expert = 1, gating = 1, gating_block = 2, expert_block = 6, all = 8))
  expect_equal(counts[[2L]], c(expert = 1, gating = 1, gating_block = 4, expert_block = 9, all = 13))
})

swap_fit <- function(fit, permutation) {
  out <- fit
  for (set in c("expert", "expert_slope", "gating", "gating_slope", "sigma", "sigma_slope")) {
    value <- fit$coefficients[[set]]
    out$coefficients[[set]] <- if (length(dim(value)) == 3L) {
      value[, permutation, , drop = FALSE]
    } else {
      value[, permutation, drop = FALSE]
    }
    dimnames(out$coefficients[[set]]) <- dimnames(value)
  }
  out$posterior <- fit$posterior[, permutation, , drop = FALSE]
  dimnames(out$posterior) <- dimnames(fit$posterior)
  out
}

test_that("bootstrap tests follow the reference expert component into a relabelled refit", {
  fit <- glrt_gaussian_fit(n = 50, seed = 331, u_grid = c(0.3, 0.7), bandwidth = 0.4, maxit = 8)$fit
  swapped <- swap_fit(fit, c(2L, 1L))
  expert <- VCMoE:::.vcmoe_glrt_constraint(fit, "coefficient", "expert", 1, "z1")
  same <- VCMoE:::.vcmoe_glrt_bootstrap_constraint(fit, fit, expert)
  expect_equal(same$constraint$shared_parameters, "expert:component1:z1:intercept")
  moved <- VCMoE:::.vcmoe_glrt_bootstrap_constraint(fit, swapped, expert)
  expect_equal(moved$constraint$shared_parameters, "expert:component2:z1:intercept")
  expect_equal(moved$constraint$fixed_zero_parameters, "expert:component2:z1:slope")
  expect_equal(moved$component_map, "component1->component2")
  gating <- VCMoE:::.vcmoe_glrt_constraint(fit, "constant_block", "gating", NULL, NULL)
  expect_identical(VCMoE:::.vcmoe_glrt_bootstrap_constraint(fit, swapped, gating)$constraint, gating)
})

test_that("bootstrap gating-contrast tests follow both components when k = 3", {
  skip_on_cran()
  fit <- glrt_gaussian_fit(n = 240, k = 3, seed = 3391, u_grid = c(0.3, 0.7), bandwidth = 0.45, maxit = 10)$fit
  # Refit labels a permutation of the reference components: refit component j
  # is reference component perm[j].
  perm <- c(3L, 1L, 2L)
  relabelled <- swap_fit(fit, perm)
  contrast <- VCMoE:::.vcmoe_glrt_constraint(fit, "coefficient", "gating", 2, "x1")
  mapped <- VCMoE:::.vcmoe_glrt_bootstrap_constraint(fit, relabelled, contrast)
  expect_identical(mapped$constraint, contrast)
  # The mapped refit carries the reference labels again, so the null tests the
  # same pair of components as in the reference fit.
  for (set in c("expert", "gating", "gating_slope")) {
    expect_equal(mapped$fit$coefficients[[set]], fit$coefficients[[set]])
  }
  expect_equal(mapped$fit$posterior, fit$posterior)
  null_reference <- VCMoE:::.vcmoe_glrt_null_fit(fit, contrast, control = list(maxit = 20))
  null_mapped <- VCMoE:::.vcmoe_glrt_null_fit(mapped$fit, mapped$constraint, control = list(maxit = 20))
  expect_equal(VCMoE:::.vcmoe_glrt_statistic(fit, null_reference)$lambda,
               VCMoE:::.vcmoe_glrt_statistic(mapped$fit, null_mapped)$lambda, tolerance = 1e-6)
  # Without the mapping the refit's contrast 2 would compare other components.
  unmapped <- VCMoE:::.vcmoe_glrt_null_fit(relabelled, contrast, control = list(maxit = 20))
  expect_false(isTRUE(all.equal(VCMoE:::.vcmoe_glrt_statistic(fit, null_reference)$lambda,
                                VCMoE:::.vcmoe_glrt_statistic(relabelled, unmapped)$lambda,
                                tolerance = 1e-6)))
})

test_that("bootstrap draws report counts and keep nulls that only hit the iteration limit", {
  skip_on_cran()
  case <- glrt_gaussian_fit(n = 70, seed = 3401, u_grid = c(0.3, 0.7), bandwidth = 0.45, maxit = 20)
  glrt <- suppressWarnings(vcmoe_glrt(
    case$fit, case$sim$data, test = "coefficient", coefficient_set = "expert", component = 1, term = "z1",
    calibration = "bootstrap", B = 3, seed = 11, control = list(maxit = 2),
    refit_control = list(maxit = 10, n_starts = 1)
  ))
  counts <- glrt$bootstrap_counts
  expect_identical(counts$requested, 3L)
  expect_identical(counts$successful + counts$failed + counts$excluded, 3L)
  expect_true(all(c("null_converged", "component_map", "lambda_fit_objective") %in% names(glrt$replicate_summary)))
  expect_gt(counts$null_not_converged_kept, 0L)
  expect_true(is.finite(glrt$bootstrap_p_value))
})

test_that("joint-path fits start from the u-free mixture posterior by default", {
  sim <- simulate_vcmoe_gaussian(n = 80, k = 2, seed = 9011, separation = 1.5)
  fit_with <- function(global) {
    control <- list(maxit = 1, n_starts = 1, seed = 9012, warn_ambiguous = FALSE)
    if (!is.null(global)) control$joint_path_global_start <- global
    suppressWarnings(vcmoe_fit(
      y ~ z1 | x1, data = sim$data, u = "u", bandwidth = 0.45,
      u_grid = seq(0.2, 0.8, length.out = 3), engine = "joint_path_em", control = control
    ))
  }
  default_fit <- fit_with(NULL)
  global_fit <- fit_with(TRUE)
  random_fit <- fit_with(FALSE)
  expect_true(isTRUE(default_fit$control$joint_path_global_start))
  expect_equal(coef(default_fit, "expert"), coef(global_fit, "expert"))
  expect_false(isTRUE(all.equal(coef(global_fit, "expert"), coef(random_fit, "expert"))))
  spec <- VCMoE:::.vcmoe_dev_parameter_spec(global_fit, 1L)
  tau <- VCMoE:::.joint_path_global_tau(spec$y, spec$trials, global_fit$family, spec$y, spec$z_design,
                                        spec$x_design, spec$offset, 2L, spec$control, spec$estimation_spec)
  expect_equal(dim(tau), c(nrow(sim$data), 2L))
  expect_equal(as.numeric(rowSums(tau)), rep(1, nrow(sim$data)), tolerance = 1e-7)
})

test_that("observation-level likelihoods bound theta as the fit does", {
  sim <- simulate_vcmoe_negbin(n = 150, k = 2, seed = 3501, separation = 1.5)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 + offset(log_size_factor) | x1, data = sim$data, u = "u", family = "negative-binomial",
    bandwidth = 0.4, u_grid = c(0.3, 0.5, 0.7),
    control = list(maxit = 50, n_starts = 1, seed = 3502, warn_ambiguous = FALSE)
  ))
  capped <- fit
  capped$coefficients$theta[1L, 1L] <- fit$control$negbin_theta_max
  unbounded <- fit
  unbounded$coefficients$theta[1L, 1L] <- Inf
  expect_true(is.finite(VCMoE:::.vcmoe_glrt_observation_loglik(unbounded)))
  expect_equal(VCMoE:::.vcmoe_glrt_observation_loglik(unbounded),
               VCMoE:::.vcmoe_glrt_observation_loglik(capped))
})
