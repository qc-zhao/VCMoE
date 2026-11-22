cluster_test_fit <- function(drop_row = FALSE, k = 2L) {
  sim <- simulate_vcmoe_gaussian(
    n = 64L, k = k, seed = 9101, separation = 2,
    u = seq(0, 1, length.out = 64L)
  )
  rownames(sim$data) <- paste0("cell_", seq_len(nrow(sim$data)))
  if (drop_row) sim$data$z1[7L] <- NA_real_
  suppressWarnings(vcmoe_fit(
    y ~ z1 | x1, data = sim$data, u = "u", k = k,
    bandwidth = .75, u_grid = c(.3, .7),
    control = list(maxit = 3L, n_starts = 1L, seed = 9102,
                   warn_ambiguous = FALSE, keep_data = TRUE)
  ))
}

cluster_test_moments <- function(fit, data = NULL, grid_id, control = list()) {
  spec <- VCMoE:::.vcmoe_dev_parameter_spec(fit, grid_id)
  par <- VCMoE:::.vcmoe_dev_parameter_vector(fit, grid_id)
  n <- length(spec$u)
  d <- length(par)
  score <- (outer(seq_len(n), seq_len(d), function(i, j) {
    .4 + sin(i * (j + .7) / 11) + cos(i / 3 + j)
  })) * spec$weights
  colnames(score) <- names(par)
  bread <- diag(seq(2, 3, length.out = d)) + matrix(.015, d, d)
  inv <- solve(bread)
  list(status = "ok", block_reason = NA_character_, parameter = par,
       spec = spec, score = score, bread = bread, meat = crossprod(score),
       covariance = inv %*% crossprod(score) %*% t(inv),
       condition = kappa(bread), min_eigenvalue = min(eigen(bread)$values),
       covariance_adjustment = "HC0", covariance_factor = 1,
       local_effective_n = VCMoE:::.vcmoe_dev_effective_n(spec$weights),
       leverage_mean = NA_real_, leverage_max = NA_real_)
}

# Pins the original CR0 + Gaussian multiplier algebra; the cluster default is CR2 + wild.
cluster_test_control <- function() list(multiplier_B = 127L, seed = 9103, covariance_adjustment = "CR0",
                                        simultaneous_method = "multiplier_coefficient_path")

test_that("cluster IDs obey the retained-row contract without guessing", {
  fit <- cluster_test_fit(drop_row = TRUE)
  n <- length(fit$fitted$u)
  ids <- rep(c("donor_b", "donor_a", "donor_c"), length.out = n)
  named <- stats::setNames(ids, rownames(fit$fitted$z_design))
  expect_equal(VCMoE:::.vcmoe_cluster_index(fit, ids),
               VCMoE:::.vcmoe_cluster_index(fit, named))
  expect_equal(VCMoE:::.vcmoe_cluster_index(fit, factor(ids)),
               match(ids, unique(ids)))
  expect_equal(VCMoE:::.vcmoe_cluster_index(fit, NULL), seq_len(n))
  expect_false(7L %in% fit$rows_used)
  expect_error(VCMoE:::.vcmoe_cluster_index(fit, rep(1:2, 32L)), "retained fitted row")
  expect_error(VCMoE:::.vcmoe_cluster_index(fit, named[rev(seq_len(n))]), "exactly the same order")
  for (bad in list(rep(NA_character_, n), rep(" ", n), rep(Inf, n),
                   matrix(ids, ncol = 1L), as.list(ids), ids[-1L])) {
    expect_error(VCMoE:::.vcmoe_cluster_index(fit, bad), "retained fitted row")
  }
  expect_error(VCMoE:::.vcmoe_cluster_index(fit, rep("one", n)), "at least two")
  expect_identical(tail(names(formals(vcmoe_confband)), 1L), "cluster")
})

test_that("cluster design diagnostics use weight sums, not cluster averages", {
  design <- VCMoE:::.vcmoe_cluster_design(c(1, 2, 0, 3, 4, 0), c(1, 1, 2, 3, 3, 4))
  expect_equal(design$n_clusters, 4L)
  expect_equal(design$active_clusters, 2L)
  expect_equal(design$effective_clusters, 100 / 58)
  expect_equal(design$max_cluster_weight_share, .7)
})

test_that("CR0 matches uncentered weighted cluster-score algebra", {
  fit <- cluster_test_fit()
  ids <- rep(seq_len(16L), each = 4L)
  testthat::local_mocked_bindings(.package = "VCMoE",
    .vcmoe_dev_score_hessian = cluster_test_moments)
  band <- suppressWarnings(vcmoe_confband(fit, cluster = ids, strict = FALSE,
                                         control = cluster_test_control()))
  for (j in seq_along(fit$u_grid)) {
    sh <- cluster_test_moments(fit, grid_id = j)
    sum_score <- rowsum(sh$score, ids, reorder = FALSE)
    inv <- solve(sh$bread)
    expected <- inv %*% crossprod(sum_score) %*% t(inv)
    rows <- band$intervals[band$intervals$grid_id == j, ]
    expect_equal(rows$se^2, unname(diag(expected)), tolerance = 1e-10)
    expect_equal(rows$estimate, unname(sh$parameter), tolerance = 1e-12)
    expect_equal(rows$covariance_factor, rep(1, nrow(rows)))
    expect_true(all(rows$pointwise_status == "ok"))
    centered <- scale(sum_score, center = TRUE, scale = FALSE)
    expect_gt(max(abs(expected - inv %*% crossprod(centered) %*% t(inv))), .01)
  }
  named <- stats::setNames(ids, rownames(fit$fitted$z_design))
  named_band <- suppressWarnings(vcmoe_confband(fit, cluster = named, strict = FALSE,
                                               control = cluster_test_control()))
  expect_equal(band$intervals, named_band$intervals)
  relabeled <- suppressWarnings(vcmoe_confband(fit, cluster = paste0("e", 100 - ids),
                                               strict = FALSE, control = cluster_test_control()))
  expect_equal(band$intervals$se, relabeled$intervals$se)
  expect_equal(band$intervals$simultaneous_critical, relabeled$intervals$simultaneous_critical)
  expect_identical(band$settings$covariance_adjustment, "CR0")
  expect_identical(band$settings$simultaneous_method, "multiplier_coefficient_path")
  expect_equal(band$settings$multiplier_B, 127L)
  expect_true(all(band$intervals$few_clusters))
  expect_true(all(band$cluster_support$effective_clusters > 0))
})

test_that("singleton clusters reproduce HC0 covariance and path bands", {
  fit <- cluster_test_fit()
  testthat::local_mocked_bindings(.package = "VCMoE",
    .vcmoe_dev_score_hessian = cluster_test_moments)
  ctl <- utils::modifyList(cluster_test_control(), list(covariance_adjustment = "HC0"))
  hc <- vcmoe_confband(fit, strict = FALSE, control = ctl)
  cr <- suppressWarnings(vcmoe_confband(fit, cluster = seq_along(fit$fitted$u),
                                       strict = FALSE, control = cluster_test_control()))
  cols <- c("estimate", "se", "pointwise_lower", "pointwise_upper",
            "simultaneous_lower", "simultaneous_upper", "simultaneous_critical")
  expect_equal(hc$intervals[cols], cr$intervals[cols], tolerance = 1e-10)
  sh <- cluster_test_moments(fit, grid_id = 1L)
  expect_equal(hc$intervals$se[hc$intervals$grid_id == 1L]^2,
               unname(diag(sh$covariance)), tolerance = 1e-10)
})

test_that("paired row and donor permutations leave cluster variances unchanged", {
  fit <- cluster_test_fit()
  ids <- rep(1:16, each = 4L)
  perm <- rev(seq_along(ids))
  permuted <- fit
  for (field in c("y", "u", "u_original", "expert_offset")) {
    permuted$fitted[[field]] <- fit$fitted[[field]][perm]
  }
  permuted$fitted$response$y <- fit$fitted$response$y[perm]
  permuted$fitted$z_design <- fit$fitted$z_design[perm, , drop = FALSE]
  permuted$fitted$x_design <- fit$fitted$x_design[perm, , drop = FALSE]
  permuted$test_permutation <- perm
  mocked <- function(fit, data = NULL, grid_id, control = list()) {
    original <- fit
    if (!is.null(fit$test_permutation)) {
      for (field in c("y", "u", "u_original", "expert_offset")) {
        original$fitted[[field]] <- fit$fitted[[field]][perm]
      }
      original$fitted$response$y <- fit$fitted$response$y[perm]
      original$fitted$z_design <- fit$fitted$z_design[perm, , drop = FALSE]
      original$fitted$x_design <- fit$fitted$x_design[perm, , drop = FALSE]
    }
    sh <- cluster_test_moments(original, grid_id = grid_id)
    if (!is.null(fit$test_permutation)) sh$score <- sh$score[perm, , drop = FALSE]
    sh
  }
  testthat::local_mocked_bindings(.package = "VCMoE", .vcmoe_dev_score_hessian = mocked)
  original_band <- suppressWarnings(vcmoe_confband(fit, cluster = ids, strict = FALSE,
                                                   control = cluster_test_control()))
  permuted_band <- suppressWarnings(vcmoe_confband(permuted, cluster = ids[perm], strict = FALSE,
                                                   control = cluster_test_control()))
  expect_equal(original_band$intervals$se, permuted_band$intervals$se, tolerance = 1e-10)
  expect_equal(original_band$intervals$effective_clusters, permuted_band$intervals$effective_clusters)
})

test_that("singular cluster meat preserves positive scalar variances", {
  fit <- cluster_test_fit()
  testthat::local_mocked_bindings(.package = "VCMoE",
    .vcmoe_dev_score_hessian = cluster_test_moments)
  band <- suppressWarnings(vcmoe_confband(fit, cluster = rep(1:2, each = 32L),
                                         strict = FALSE, control = cluster_test_control()))
  expect_true(all(band$diagnostics$meat_rank <= 2L))
  expect_true(all(band$diagnostics$covariance_rank <= 2L))
  expect_true(all(band$diagnostics$parameter_count > 2L))
  expect_true(all(is.finite(band$intervals$se) & band$intervals$se > 0))
  expect_true(all(band$intervals$path_complete))
})

test_that("one blocked grid blocks full paths without erasing valid pointwise intervals", {
  fit <- cluster_test_fit()
  testthat::local_mocked_bindings(.package = "VCMoE",
    .vcmoe_dev_score_hessian = cluster_test_moments,
    .vcmoe_dev_grid_block_reason = function(fit, grid_id, control) {
      if (grid_id == 2L) "ambiguous_label" else ""
    })
  band <- suppressWarnings(vcmoe_confband(fit, cluster = rep(1:16, each = 4L),
                                         type = "simultaneous", control = cluster_test_control()))
  expect_true(all(band$intervals$status == "blocked"))
  expect_true(all(!band$intervals$path_complete))
  expect_true(all(is.na(band$intervals$lower)))
  first <- band$intervals$grid_id == 1L
  expect_true(all(band$intervals$pointwise_status[first] == "ok"))
  expect_true(all(is.finite(band$intervals$pointwise_lower[first])))
  expect_true(all(band$intervals$pointwise_block_reason[!first] == "ambiguous_label"))
  expect_true(all(band$intervals$simultaneous_block_reason == "incomplete_grid_path"))
  point <- suppressWarnings(vcmoe_confband(fit, cluster = rep(1:16, each = 4L),
                                          type = "pointwise", control = cluster_test_control()))
  expect_true(all(point$intervals$status[first] == "ok"))
  expect_true(all(point$intervals$status[!first] == "blocked"))
})

test_that("local cluster support and bread failures remain explicitly blocked", {
  fit <- cluster_test_fit()
  fit$bandwidth <- .03
  fit$u_grid <- fit$fitted$u[c(15L, 45L)]
  calls <- 0L
  mocked <- function(fit, data = NULL, grid_id, control = list()) {
    calls <<- calls + 1L
    cluster_test_moments(fit, grid_id = grid_id)
  }
  testthat::local_mocked_bindings(.package = "VCMoE", .vcmoe_dev_score_hessian = mocked)
  band <- suppressWarnings(vcmoe_confband(fit, cluster = rep(1:2, each = 32L),
    strict = FALSE, control = cluster_test_control()))
  expect_equal(calls, 0L)
  expect_true(all(band$intervals$pointwise_block_reason == "fewer_than_two_active_clusters"))
  expect_true(all(is.na(band$intervals$lower)))
})

test_that("zero variance affects its own scalar path and not other coefficients", {
  fit <- cluster_test_fit()
  mocked <- function(fit, data = NULL, grid_id, control = list()) {
    sh <- cluster_test_moments(fit, grid_id = grid_id)
    sh$bread <- diag(length(sh$parameter))
    if (grid_id == 2L) sh$score[, 1L] <- rep(c(-1, 1), nrow(sh$score) / 2L)
    sh
  }
  testthat::local_mocked_bindings(.package = "VCMoE", .vcmoe_dev_score_hessian = mocked)
  band <- suppressWarnings(vcmoe_confband(fit, cluster = rep(1:32, each = 2L),
    type = "simultaneous", strict = FALSE, control = cluster_test_control()))
  key <- band$intervals$parameter[1L]
  selected <- band$intervals$parameter == key
  expect_true(all(band$intervals$simultaneous_status[selected] == "blocked"))
  expect_true(all(band$intervals$simultaneous_status[!selected] == "ok"))
  expect_equal(band$intervals$pointwise_status[selected], c("ok", "blocked"))
  expect_equal(band$intervals$pointwise_block_reason[selected][2L], "invalid_scalar_variance")
})

test_that("shared Gaussian multipliers calibrate each curve across grid points", {
  fn <- VCMoE:::.vcmoe_path_multiplier
  identical_path <- cbind(c(1, 0), c(1, 0))
  dependent <- fn(identical_path, c("a", "a"), .95, 12000L, seed = 9104)
  independent <- fn(diag(2), c("a", "a"), .95, 12000L, seed = 9104)
  separate <- fn(diag(2), c("a", "b"), .95, 12000L, seed = 9104)
  expect_equal(unname(dependent), qnorm(.975), tolerance = .05)
  expect_equal(unname(independent), qnorm((1 + sqrt(.95)) / 2), tolerance = .05)
  expect_gt(independent[[1L]] - dependent[[1L]], .15)
  expect_equal(unname(separate), rep(qnorm(.975), 2L), tolerance = .05)
  gram_path <- rbind(identical_path, matrix(0, 4L, 2L))
  gram <- fn(gram_path, c("a", "a"), .95, 12000L, seed = 9104)
  expect_equal(unname(gram), qnorm(.975), tolerance = .05)
  expect_identical(fn(diag(2), c("a", "a"), .95, 129L, 9),
                   fn(diag(2), c("a", "a"), .95, 129L, 9))
})

test_that("seeded multiplier preserves both existing and absent RNG state", {
  set.seed(9105)
  original <- .Random.seed
  VCMoE:::.vcmoe_path_multiplier(diag(2), c("a", "a"), .95, 5L, 3)
  expect_identical(.Random.seed, original)
  on.exit(assign(".Random.seed", original, envir = .GlobalEnv), add = TRUE)
  rm(".Random.seed", envir = .GlobalEnv)
  VCMoE:::.vcmoe_path_multiplier(diag(2), c("a", "a"), .95, 5L, 3)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  for (bad in list(0, 1, 2.5, Inf, NA_real_, "10")) {
    expect_error(VCMoE:::.vcmoe_path_multiplier(diag(2), c("a", "a"), .95, bad, 3), "multiplier_B")
  }
  expect_error(VCMoE:::.vcmoe_path_multiplier(diag(2), c("a", "a"), .95, 5L, -1), "seed")
})

test_that("public cluster mode rejects incompatible fit types and legacy critical values", {
  fit <- cluster_test_fit()
  ids <- rep(1:16, each = 4L)
  for (change in c("joint", "reduced")) {
    altered <- fit
    if (change == "joint") altered$engine_id <- "experimental_joint_path_em"
    if (change == "reduced") altered$diagnostics$reduced_fit <- TRUE
    expect_error(vcmoe_confband(altered, cluster = ids), "unrestricted local-grid")
  }
  for (ctl in list(list(simultaneous_method = "analytic"),
                   list(simultaneous_method = "multiplier_grid"),
                   list(covariance_adjustment = "HC0"))) {
    expect_error(vcmoe_confband(fit, cluster = ids, control = ctl), "requires CR0")
  }
  expect_error(vcmoe_confband(fit, cluster = ids, level = 1), "Confidence|confidence|level")
})

test_that("paired comparison reuses moments and preserves old analytic numbers", {
  fit <- cluster_test_fit()
  calls <- 0L
  counted <- function(fit, data = NULL, grid_id, control = list()) {
    calls <<- calls + 1L
    cluster_test_moments(fit, grid_id = grid_id)
  }
  testthat::local_mocked_bindings(.package = "VCMoE", .vcmoe_dev_score_hessian = counted)
  old <- vcmoe_confband(fit, type = "simultaneous", strict = FALSE)
  expect_identical(old$settings$covariance_adjustment, "HC0")
  expect_identical(old$settings$simultaneous_method, "analytic_epanechnikov_path")
  calls <- 0L
  stable_fields <- c("coefficients", "fitted", "diagnostics", "posterior", "u_grid", "bandwidth")
  snapshot <- fit[stable_fields]
  all <- suppressWarnings(VCMoE:::.vcmoe_confband_comparison(
    fit, rep(1:16, each = 4L), strict = FALSE, control = cluster_test_control()))
  expect_equal(calls, length(fit$u_grid))
  expect_identical(fit[stable_fields], snapshot)
  numeric_cols <- c("estimate", "se", "pointwise_lower", "pointwise_upper",
                    "simultaneous_lower", "simultaneous_upper", "simultaneous_critical")
  expect_equal(all$hc0_analytic$intervals[numeric_cols], old$intervals[numeric_cols], tolerance = 1e-12)
  expect_identical(all$cluster_multiplier$settings$covariance_adjustment, "CR0")
  expect_identical(all$hc0_multiplier$settings$covariance_adjustment, "HC0")
})

test_that("K2 and K3 gating intervals use identifiable contrasts", {
  testthat::local_mocked_bindings(.package = "VCMoE", .vcmoe_dev_score_hessian = cluster_test_moments)
  for (k in 2:3) {
    fit <- cluster_test_fit(k = k)
    band <- suppressWarnings(vcmoe_confband(fit, cluster = rep(1:16, each = 4L),
      coefficient_set = "gating", strict = FALSE, control = cluster_test_control()))
    for (j in seq_along(fit$u_grid)) {
      coefficients <- matrix(fit$coefficients$gating[j, , ], nrow = k)
      slopes <- matrix(fit$coefficients$gating_slope[j, , ], nrow = k)
      expected <- VCMoE:::.vcmoe_dev_gating_contrasts(coefficients, slopes)
      rows <- band$intervals[band$intervals$grid_id == j & band$intervals$block == "intercept", ]
      expect_equal(rows$estimate, as.vector(t(expected$intercept)), tolerance = 1e-12)
      expect_equal(unique(rows$component), expected$labels)
    }
    expect_true(all(band$intervals$coefficient_set == "gating_contrast"))
  }
})

test_that("real Gaussian, binomial and NB offset fits exercise cluster inference", {
  for (case in c("gaussian", "gaussian_k3", "binomial", "negative-binomial")) {
    k <- if (case == "gaussian_k3") 3L else 2L
    family <- if (case == "gaussian_k3") "gaussian" else case
    n <- 120L
    u <- rep(seq(0, 1, length.out = 20L), each = 6L)
    sim <- switch(family,
      gaussian = simulate_vcmoe_gaussian(n, k, seed = 9111, separation = 2, u = u),
      binomial = simulate_vcmoe_binomial(n, k, seed = 9112, separation = 1.4, u = u, trials = 8L),
      `negative-binomial` = simulate_vcmoe_negbin(n, k, seed = 9113, separation = 1.2, u = u, mean_count = 8))
    formula <- switch(family,
      gaussian = y ~ z1 | x1,
      binomial = cbind(success, failure) ~ z1 | x1,
      `negative-binomial` = y ~ z1 + offset(log_size_factor) | x1)
    fit <- suppressWarnings(vcmoe_fit(formula, data = sim$data, u = "u", k = k,
      family = family, bandwidth = 1.1, u_grid = .5,
      control = list(maxit = 20L, n_starts = 1L, seed = 9114, warn_ambiguous = FALSE)))
    stable_fields <- c("coefficients", "fitted", "diagnostics", "posterior", "u_grid", "bandwidth")
    snapshot <- fit[stable_fields]
    band <- suppressWarnings(vcmoe_confband(fit, cluster = rep(1:20, each = 6L),
      strict = FALSE, control = list(multiplier_B = 31L, seed = 9115, max_hessian_dim = 60L)))
    expect_s3_class(band, "vcmoe_confband")
    expect_identical(fit[stable_fields], snapshot)
    expect_true(nrow(band$intervals) > 0L)
    ok <- band$intervals$pointwise_status == "ok"
    expect_true(all(is.finite(band$intervals$estimate)))
    expect_true(all(is.finite(band$intervals$se[ok]) & band$intervals$se[ok] > 0))
    expect_true(all(is.finite(band$intervals$pointwise_lower[ok])))
    expect_true(all(nzchar(band$intervals$pointwise_block_reason[!ok])))
    expect_identical(band$settings$covariance_adjustment, "CR2")
    expect_identical(band$settings$simultaneous_method, "wild_cluster_path")
    expect_true(all(is.finite(band$intervals$df[ok]) & band$intervals$df[ok] > 0))
    if (family == "gaussian" && k == 2L) expect_true(any(ok))
    if (family == "negative-binomial") {
      expect_true(any(band$intervals$term == "log_theta"))
      expect_true(any(abs(fit$fitted$expert_offset) > 0))
    }
  }
})
