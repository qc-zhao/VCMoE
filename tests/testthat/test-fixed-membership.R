fixed_test_fit <- function(dat, membership, u_grid = c(0.3, 0.7), ...) {
  vcmoe_fit(y ~ z1 + offset(log_size_factor) | x1, data = dat, u = "u", k = 2L,
            family = "negative-binomial", bandwidth = 0.3, u_grid = u_grid,
            control = list(fixed_membership = membership, maxit = 500, tol = 1e-10,
                           gating_maxit = 200, negbin_theta_ridge = 0, ...))
}

test_that("known-label local fits are kernel-weighted NB and logistic regressions", {
  skip_if_not_installed("MASS")
  dat <- simulate_vcmoe_negbin(n = 600L, k = 2L, seed = 31L, separation = 1.2)$data
  fit <- fixed_test_fit(dat, dat$component)
  expect_identical(vcmoe_parameterization(fit)$membership, "fixed")
  expect_identical(fit$diagnostics$alignment_method, "none")
  expect_true(all(fit$diagnostics$converged))
  expect_identical(fit$diagnostics$selected_start, c(1L, 1L))
  expect_equal(unname(fit$posterior[, , 2L]), fit$control$fixed_membership)
  spec <- .vcmoe_estimation_spec("a1_epanechnikov_scaled")
  for (g in seq_along(fit$u_grid)) {
    w <- .local_weights(fit$fitted$u, fit$u_grid[g], fit$bandwidth, spec)
    du <- .local_du(fit$fitted$u, fit$u_grid[g], fit$bandwidth, spec)
    keep <- w > 0
    local <- data.frame(y = dat$y, z1 = dat$z1, x1 = dat$x1, du = du, off = dat$log_size_factor, w = w)[keep, ]
    for (component in 1:2) {
      rows <- dat$component[keep] == component
      ref <- suppressWarnings(MASS::glm.nb(y ~ z1 + du + I(du * z1) + offset(off), data = local[rows, ],
                                           weights = w, control = stats::glm.control(epsilon = 1e-12, maxit = 200)))
      expect_equal(unname(c(fit$coefficients$expert[g, component, ], fit$coefficients$expert_slope[g, component, ])),
                   unname(stats::coef(ref)), tolerance = 1e-3)
      expect_equal(unname(fit$coefficients$theta[g, component]), ref$theta, tolerance = 1e-3)
    }
    local$h <- as.numeric(dat$component[keep] == 1L)
    gate <- suppressWarnings(stats::glm(h ~ x1 + du + I(du * x1), family = stats::binomial(), data = local,
                                        weights = w, control = stats::glm.control(epsilon = 1e-12)))
    contrast <- c(fit$coefficients$gating[g, 1L, ] - fit$coefficients$gating[g, 2L, ],
                  fit$coefficients$gating_slope[g, 1L, ] - fit$coefficients$gating_slope[g, 2L, ])
    expect_equal(unname(contrast), unname(stats::coef(gate)), tolerance = 1e-5)
  }
})

test_that("known-label inference uses the complete-data local likelihood", {
  dat <- simulate_vcmoe_negbin(n = 400L, k = 2L, seed = 32L, separation = 1.2)$data
  fit <- fixed_test_fit(dat, dat$component, u_grid = 0.5)
  spec <- .vcmoe_dev_parameter_spec(fit, 1L)
  par <- .vcmoe_dev_parameter_vector(fit, 1L)
  gating_design <- .make_expert_design(spec$x_design, spec$du)
  prior <- .gating_prob(gating_design, cbind(fit$coefficients$gating[1L, , ], fit$coefficients$gating_slope[1L, , ]))
  expert_design <- .make_expert_design(spec$z_design, spec$du)
  loglik <- 0
  score <- list()
  for (component in 1:2) {
    beta <- c(fit$coefficients$expert[1L, component, ], fit$coefficients$expert_slope[1L, component, ])
    mu <- exp(spec$offset + as.numeric(expert_design %*% beta))
    theta <- fit$coefficients$theta[1L, component]
    rows <- dat$component == component
    loglik <- loglik + sum((spec$weights * (log(prior[, component]) +
                                              stats::dnbinom(spec$y, size = theta, mu = mu, log = TRUE)))[rows])
    score[[component]] <- spec$weights * rows * (spec$y - mu) / (1 + mu / theta) * expert_design
  }
  expect_equal(.vcmoe_dev_local_loglik(fit, grid_id = 1L, par = par), loglik, tolerance = 1e-10)
  # Parameter order: intercept block for each component, then slope block.
  expected <- cbind(score[[1L]][, 1:2], score[[2L]][, 1:2], score[[1L]][, 3:4], score[[2L]][, 3:4])
  expert <- spec$parameter_table$coefficient_set == "expert"
  expect_equal(unname(.vcmoe_dev_score_matrix(par, spec)[, expert]), unname(expected), tolerance = 1e-6)
})

test_that("membership given as indices, a factor or a matrix gives the same fit", {
  dat <- simulate_vcmoe_negbin(n = 300L, k = 2L, seed = 33L, separation = 1.2)$data
  by_index <- fixed_test_fit(dat, dat$component, u_grid = 0.5)
  by_factor <- fixed_test_fit(dat, factor(c("a", "b")[dat$component]), u_grid = 0.5)
  by_matrix <- fixed_test_fit(dat, cbind(dat$component == 1L, dat$component == 2L) + 0, u_grid = 0.5)
  expect_identical(by_factor$coefficients, by_index$coefficients)
  expect_identical(by_matrix$coefficients, by_index$coefficients)
})

test_that("invalid or unsupported fixed membership fails clearly", {
  dat <- simulate_vcmoe_negbin(n = 120L, k = 2L, seed = 34L)$data
  expect_error(fixed_test_fit(dat, dat$component[-1L]), "component in 1..k")
  expect_error(fixed_test_fit(dat, replace(dat$component, 1L, 3L)), "component in 1..k")
  expect_error(fixed_test_fit(dat, rep(1L, nrow(dat))), "at least one row")
  expect_error(fixed_test_fit(dat, matrix(1, nrow(dat), 2L)), "one 1 per row")
  expect_error(
    vcmoe_fit(y ~ z1 + offset(log_size_factor) | x1, data = dat, u = "u", k = 2L,
              family = "negative-binomial", engine = "joint_path_em",
              control = list(fixed_membership = dat$component)),
    "local_grid_em"
  )
})

test_that("latent fits report latent membership", {
  dat <- simulate_vcmoe_gaussian(n = 150L, k = 2L, seed = 35L)$data
  fit <- suppressWarnings(vcmoe_fit(y ~ z1 | x1, data = dat, u = "u", k = 2L, u_grid = 0.5,
                                    control = list(seed = 1, maxit = 20)))
  expect_identical(vcmoe_parameterization(fit)$membership, "latent")
})

test_that("known-label fits give cluster-robust bands", {
  dat <- simulate_vcmoe_negbin(n = 400L, k = 2L, seed = 36L, separation = 1.2)$data
  fit <- fixed_test_fit(dat, dat$component, u_grid = c(0.3, 0.5, 0.7))
  ids <- rep(seq_len(20L), length.out = nrow(dat))
  band <- suppressWarnings(vcmoe_confband(fit, cluster = ids, type = "simultaneous",
                                          control = list(multiplier_B = 99L, seed = 1)))
  expect_true(all(band$intervals$status[band$intervals$coefficient_set == "expert"] == "ok"))
  expect_true(all(band$settings$covariance_adjustment == "CR2"))
})
