stop_rule_fit <- function(k, theta_target, tol = 1e-8, seed = 4401, maxit = 2000) {
  sim <- simulate_vcmoe_negbin(n = 300, k = k, seed = seed, separation = 1.5)
  suppressWarnings(vcmoe_fit(
    y ~ z1 + offset(log_size_factor) | x1,
    data = sim$data, u = "u", family = "negative-binomial", k = k,
    bandwidth = 0.4, u_grid = c(0.3, 0.5, 0.7),
    control = list(maxit = maxit, tol = tol, n_starts = 1, seed = seed + 1, warn_ambiguous = FALSE,
                   negbin_theta_ridge = 0.05, negbin_theta_target = theta_target)
  ))
}

converged_stop_rule_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) fit <<- stop_rule_fit(2L, 8, tol = 1e-6)
    fit
  }
})

stop_rule_params <- function(fit, g) {
  list(
    expert_coef = cbind(matrix(fit$coefficients$expert[g, , ], nrow = fit$k),
                        matrix(fit$coefficients$expert_slope[g, , ], nrow = fit$k)),
    gating_coef = cbind(matrix(fit$coefficients$gating[g, , ], nrow = fit$k),
                        matrix(fit$coefficients$gating_slope[g, , ], nrow = fit$k)),
    sigma = NULL, sigma_slope = NULL,
    theta = as.numeric(fit$coefficients$theta[g, ])
  )
}

test_that("the local EM penalty is the penalty of the inference objective", {
  skip_on_cran()
  for (case in list(list(k = 2L, target = 8), list(k = 3L, target = 8), list(k = 3L, target = NULL))) {
    fit <- stop_rule_fit(case$k, case$target)
    control <- VCMoE:::.vcmoe_default_control(VCMoE:::.vcmoe_legacy_fit_control(fit$control))
    for (g in seq_along(fit$u_grid)) {
      par <- VCMoE:::.vcmoe_dev_parameter_vector(fit, g)
      bread_penalty <- VCMoE:::.vcmoe_dev_local_loglik(fit, grid_id = g, par = par, penalized = FALSE) -
        VCMoE:::.vcmoe_dev_local_loglik(fit, grid_id = g, par = par, penalized = TRUE)
      spec <- VCMoE:::.vcmoe_dev_parameter_spec(fit, g)
      em_penalty <- VCMoE:::.local_ridge_penalty("negative-binomial", stop_rule_params(fit, g), spec$weights,
                                                 posterior = NULL, control = control)
      expect_gt(em_penalty, 0)
      expect_equal(em_penalty, bread_penalty, tolerance = 1e-8)
    }
  }
})

test_that("the local EM stops at the first step whose penalized change is within tol", {
  skip_on_cran()
  # The EM is deterministic, so rerunning it with maxit = n - 1 and n - 2
  # recovers the iterates before the stop at iteration n.
  tol <- 1e-6
  for (case in list(list(k = 2L, seed = 4701), list(k = 3L, seed = 4501), list(k = 3L, seed = 4401))) {
    fit <- stop_rule_fit(case$k, 8, tol = 1e-4, seed = case$seed)
    for (g in seq_along(fit$u_grid)) {
      spec <- VCMoE:::.vcmoe_dev_parameter_spec(fit, g)
      init <- stop_rule_params(fit, g)
      run <- function(maxit) {
        control <- utils::modifyList(spec$control, list(maxit = as.integer(maxit), tol = tol))
        VCMoE:::.local_em(spec$y, spec$trials, fit$family, NULL, spec$z_design, spec$x_design, spec$offset,
                          spec$u, spec$u0, fit$bandwidth, fit$k, control, spec$estimation_spec,
                          init_params = init)
      }
      objective <- function(local) {
        params <- VCMoE:::.local_to_params(local)
        VCMoE:::.local_loglik(spec$y, spec$trials, fit$family, spec$z_design, spec$x_design, spec$du,
                              spec$offset, spec$weights, params, spec$control) -
          VCMoE:::.local_ridge_penalty(fit$family, params, spec$weights, NULL, spec$control)
      }
      stopped <- run(5000L)
      expect_true(stopped$converged)
      expect_identical(stopped$stop_reason, "tolerance")
      expect_equal(stopped$objective, objective(stopped), tolerance = 1e-10)
      n <- stopped$iterations
      if (n < 3L) next
      before <- run(n - 1L)
      earlier <- run(n - 2L)
      expect_false(before$converged)
      expect_identical(before$stop_reason, "maxit")
      threshold <- function(value) tol * (abs(value) + 1)
      expect_lte(abs(objective(stopped) - objective(before)), threshold(objective(before)))
      expect_gt(abs(objective(before) - objective(earlier)), threshold(objective(earlier)))
    }
  }
})

test_that("start selection prefers converged starts and then the penalized objective", {
  make <- function(objective, converged) list(objective = objective, converged = converged, loglik = objective + 1)
  candidates <- list(make(-10, TRUE), make(-5, FALSE), make(-8, TRUE))
  picked <- VCMoE:::.select_local_start(candidates, loglik = c(-9, -4, -7), objective = c(-10, -5, -8),
                                        start_converged = c(TRUE, FALSE, TRUE), u0 = 0.5)
  expect_identical(picked$selected_start, 3L)
  expect_true(picked$better_unconverged_start)
  none <- VCMoE:::.select_local_start(candidates, loglik = c(-9, -4, -7), objective = c(-10, -5, -8),
                                      start_converged = c(FALSE, FALSE, FALSE), u0 = 0.5)
  expect_identical(none$selected_start, 2L)
  expect_false(none$better_unconverged_start)
  legacy <- VCMoE:::.select_local_start(candidates, loglik = c(-9, -4, -7), objective = c(-10, -5, -8),
                                        start_converged = c(TRUE, FALSE, TRUE), u0 = 0.5, legacy = TRUE)
  expect_identical(legacy$selected_start, 2L)
  expect_error(VCMoE:::.select_local_start(candidates, loglik = rep(-Inf, 3), objective = rep(-Inf, 3),
                                           start_converged = rep(FALSE, 3), u0 = 0.5),
               "All EM starts failed")
})

test_that("fits report the stop reason, objective and start choice per grid point", {
  fit <- stop_rule_fit(2L, 8, tol = 1e-5, maxit = 3)
  d <- fit$diagnostics
  expect_identical(d$stop_reason, rep("maxit", 3L))
  expect_false(any(d$converged))
  expect_true(all(is.finite(d$objective)))
  expect_true(all(is.finite(d$objective_delta)))
  expect_length(d$start_objective, 3L)
  diag_table <- vcmoe_diagnostics(fit)
  expect_true(all(c("stop_reason", "objective_delta", "selected_start") %in% names(diag_table)))
  converged <- converged_stop_rule_fit()
  expect_identical(converged$diagnostics$stop_reason, rep("tolerance", 3L))
})

test_that("strict inference blocks only the grid points that did not converge", {
  fit <- converged_stop_rule_fit()
  fit$diagnostics$converged[2L] <- FALSE
  control <- VCMoE:::.vcmoe_dev_control(list())
  reasons <- vapply(1:3, function(g) VCMoE:::.vcmoe_dev_grid_block_reason(fit, g, control), character(1L))
  expect_match(reasons[[2L]], "grid_not_converged")
  expect_false(grepl("grid_not_converged", reasons[[1L]]))
  expect_false(grepl("grid_not_converged", reasons[[3L]]))
  expect_false(any(grepl("low_converged_fraction", reasons)))
  gated <- VCMoE:::.vcmoe_dev_control(list(min_converged_fraction = 0.9))
  expect_match(VCMoE:::.vcmoe_dev_grid_block_reason(fit, 1L, gated), "low_converged_fraction")
})

test_that("theta ridge mass is recorded and 0.2.0 objects are read with the old rule", {
  fit <- converged_stop_rule_fit()
  expect_identical(fit$control$negbin_theta_penalty_mass, "kernel")
  control <- VCMoE:::.vcmoe_dev_control(list())
  expect_false(grepl("component_mass", VCMoE:::.vcmoe_dev_grid_block_reason(fit, 2L, control)))
  old <- fit
  old$control$negbin_theta_penalty_mass <- NULL
  expect_identical(VCMoE:::.vcmoe_legacy_fit_control(old$control)$negbin_theta_penalty_mass, "component")
  expect_match(VCMoE:::.vcmoe_dev_grid_block_reason(old, 2L, control), "component_mass_theta_penalty_refit_required")
  no_ridge <- old
  no_ridge$control$negbin_theta_ridge <- 0
  expect_false(grepl("component_mass", VCMoE:::.vcmoe_dev_grid_block_reason(no_ridge, 2L, control)))
  expect_error(VCMoE:::.vcmoe_theta_penalty_mass(list(negbin_theta_penalty_mass = "other")),
               "negbin_theta_penalty_mass")
})
