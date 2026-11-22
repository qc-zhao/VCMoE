small_ols <- function(seed = 9201, G = 9L, weights = NULL) {
  set.seed(seed)
  sizes <- rep(c(3L, 5L, 8L), length.out = G)
  group <- rep(seq_len(G), sizes)
  n <- length(group)
  X <- cbind(1, stats::rnorm(n), stats::rnorm(G)[group])
  y <- as.numeric(X %*% c(1, .5, -.3) + stats::rnorm(G)[group] + stats::rnorm(n))
  w <- if (is.null(weights)) rep(1, n) else weights
  A <- crossprod(X, X * w)
  beta <- solve(A, crossprod(X, w * y))
  e <- as.numeric(y - X %*% beta)
  score <- rowsum(X * (w * e), group, reorder = FALSE)
  hessian <- array(0, c(G, 3L, 3L))
  omega <- vector("list", G)
  for (g in seq_len(G)) {
    rows <- group == g
    hessian[g, , ] <- crossprod(X[rows, , drop = FALSE], X[rows, , drop = FALSE] * w[rows])
    omega[[g]] <- crossprod(X[rows, , drop = FALSE] * w[rows])
  }
  list(X = X, y = y, w = w, e = e, group = group, A = A, score = score,
       hessian = hessian, omega = omega, data = data.frame(y = y, x1 = X[, 2], x2 = X[, 3], g = group))
}

test_that("score-space CR1/CR2/CR3 equal classical cluster sandwiches for OLS", {
  skip_if_not_installed("sandwich")
  x <- small_ols()
  model <- stats::lm(y ~ x1 + x2, data = x$data)
  # vcovCL's clustered HC2/HC3 meats carry (G-1)/G; cadjust = TRUE cancels it,
  # giving the unscaled Bell-McCaffrey CR2 and jackknife-type CR3.
  for (case in list(c("CR1", "HC0"), c("CR2", "HC2"), c("CR3", "HC3"))) {
    adj <- VCMoE:::.vcmoe_cluster_adjust(x$score, x$hessian, x$omega, x$A, case[[1L]])
    reference <- sandwich::vcovCL(model, cluster = ~g, type = case[[2L]], cadjust = TRUE)
    expect_equal(unname(crossprod(adj$contributions)), unname(reference), tolerance = 1e-10)
  }
  H <- x$X %*% solve(x$A) %*% t(x$X)
  for (case in list(c("CR2", -0.5), c("CR3", -1))) {
    meat <- Reduce(`+`, lapply(seq_len(max(x$group)), function(g) {
      rows <- which(x$group == g)
      adjust <- VCMoE:::.vcmoe_psd_power(diag(length(rows)) - H[rows, rows, drop = FALSE], as.numeric(case[[2L]]))
      tcrossprod(t(x$X[rows, , drop = FALSE]) %*% adjust %*% x$e[rows])
    }))
    adj <- VCMoE:::.vcmoe_cluster_adjust(x$score, x$hessian, x$omega, x$A, case[[1L]])
    expect_equal(crossprod(adj$contributions), solve(x$A) %*% meat %*% solve(x$A), tolerance = 1e-10,
                 ignore_attr = TRUE)
  }
  cr0 <- VCMoE:::.vcmoe_cluster_adjust(x$score, x$hessian, x$omega, x$A, "CR0")
  expect_equal(unname(crossprod(cr0$contributions)),
               unname(sandwich::vcovCL(model, cluster = ~g, type = "HC0", cadjust = FALSE)), tolerance = 1e-10)
  expect_equal(cr0$df, rep(Inf, 3))
  expect_equal(VCMoE:::.vcmoe_cluster_adjust(x$score, x$hessian, x$omega, x$A, "CR1")$df, rep(8, 3))
})

test_that("Bell-McCaffrey df matches an independent residual-space OLS computation", {
  x <- small_ols(G = 7L)
  n <- nrow(x$X)
  inv <- solve(x$A)
  H <- x$X %*% inv %*% t(x$X)
  residual <- diag(n) - H
  adj <- VCMoE:::.vcmoe_cluster_adjust(x$score, x$hessian, x$omega, x$A, "CR2")
  for (k in 1:3) {
    columns <- vapply(seq_len(max(x$group)), function(g) {
      rows <- which(x$group == g)
      root <- VCMoE:::.vcmoe_psd_power(diag(length(rows)) - H[rows, rows, drop = FALSE], -0.5)
      as.numeric(residual[, rows, drop = FALSE] %*% root %*% x$X[rows, , drop = FALSE] %*% inv[, k])
    }, numeric(n))
    Q <- crossprod(columns)
    expect_equal(adj$df[[k]], sum(diag(Q))^2 / sum(Q * Q), tolerance = 1e-8)
    expect_lte(adj$df[[k]], max(x$group))
  }
})

test_that("CR2 is exactly unbiased under its working model with kernel weights and ridge", {
  set.seed(9202)
  x <- small_ols(G = 6L, weights = stats::runif(length(rep(seq_len(6L), rep(c(3L, 5L, 8L), 2L))), .2, 2))
  ridge <- diag(c(.4, .1, .2))
  bread <- x$A + ridge
  adj <- VCMoE:::.vcmoe_cluster_adjust(x$score, x$hessian, x$omega, bread, "CR2")
  inv <- solve(bread)
  total <- Reduce(`+`, x$omega)
  expected_meat <- matrix(0, 3, 3)
  for (g in seq_len(6L)) {
    P <- x$hessian[g, , ] %*% inv
    E <- x$omega[[g]] - P %*% x$omega[[g]] - x$omega[[g]] %*% t(P) + P %*% total %*% t(P)
    M <- bread %*% adj$R[[g]]
    expect_equal(M %*% E %*% t(M), x$omega[[g]], tolerance = 1e-8)
    expected_meat <- expected_meat + adj$R[[g]] %*% E %*% t(adj$R[[g]])
  }
  expect_equal(expected_meat, inv %*% total %*% inv, tolerance = 1e-8)
  expect_true(adj$cr2_min_eigenvalue > 0 && adj$cr2_min_eigenvalue <= 1)
  expect_equal(adj$cr2_dropped_directions, 0L)
})

test_that("cluster Hessian blocks sum to the penalized numerical bread", {
  sim <- simulate_vcmoe_negbin(n = 120, k = 2, seed = 9203, separation = 1.2,
                               u = rep(seq(0, 1, length.out = 20L), each = 6L), mean_count = 8)
  fit <- suppressWarnings(vcmoe_fit(y ~ z1 + offset(log_size_factor) | x1, data = sim$data, u = "u",
    k = 2, family = "negative-binomial", bandwidth = .6, u_grid = .5,
    control = list(maxit = 20L, n_starts = 1L, seed = 9204, warn_ambiguous = FALSE,
                   negbin_theta_ridge = .05, negbin_theta_target = 8)))
  group <- rep(1:20, each = 6L)
  ctl <- VCMoE:::.vcmoe_dev_control(list(max_hessian_dim = 60L))
  plain <- VCMoE:::.vcmoe_dev_score_hessian(fit, grid_id = 1L, control = ctl)
  blocks <- VCMoE:::.vcmoe_dev_score_hessian(fit, grid_id = 1L, control = ctl, group = group)
  expect_null(plain$cluster_hessian)
  expect_equal(dim(blocks$cluster_hessian), c(20L, length(plain$parameter), length(plain$parameter)))
  expect_equal(blocks$bread, plain$bread, tolerance = 1e-6)
  expect_equal(blocks$score, plain$score)
  spec <- VCMoE:::.vcmoe_dev_parameter_spec(fit, 1L)
  inactive <- rowsum(spec$weights, group)[, 1L] == 0
  expect_true(all(blocks$cluster_hessian[inactive, , ] == 0))
})

test_that("kernel-mass theta ridge makes the fitted theta score stationary", {
  sim <- simulate_vcmoe_negbin(n = 300, k = 2, seed = 9205, separation = 1.6, mean_count = 10)
  base <- list(maxit = 400L, tol = 1e-10, n_starts = 1L, seed = 9206, warn_ambiguous = FALSE,
               negbin_theta_ridge = 2, negbin_theta_target = 50)
  gradients <- list()
  for (mass in c("kernel", "component")) {
    fit <- suppressWarnings(vcmoe_fit(y ~ z1 + offset(log_size_factor) | x1, data = sim$data, u = "u",
      k = 2, family = "negative-binomial", bandwidth = .5, u_grid = .5,
      control = utils::modifyList(base, list(negbin_theta_penalty_mass = mass))))
    expect_identical(fit$control$negbin_theta_penalty_mass, mass)
    if (mass == "kernel") {
      spec <- VCMoE:::.vcmoe_dev_parameter_spec(fit, 1L)
      par <- VCMoE:::.vcmoe_dev_parameter_vector(fit, 1L)
      rows <- spec$parameter_table$term == "log_theta"
      score <- colSums(VCMoE:::.vcmoe_dev_score_matrix(par, spec))[rows]
      penalized <- score - spec$ridge[rows] * (par[rows] - spec$ridge_target[rows])
      gradients[[mass]] <- max(abs(penalized)) / max(spec$ridge[rows])
    } else {
      # Legacy fits keep their rule; the kernel-mass bread then misstates the equation.
      spec <- VCMoE:::.vcmoe_dev_parameter_spec(fit, 1L)
      par <- VCMoE:::.vcmoe_dev_parameter_vector(fit, 1L)
      rows <- spec$parameter_table$term == "log_theta"
      score <- colSums(VCMoE:::.vcmoe_dev_score_matrix(par, spec))[rows]
      expect_identical(spec$theta_penalty_mass, "component")
      gradients[[mass]] <- max(abs(score - spec$ridge[rows] * (par[rows] - spec$ridge_target[rows]))) /
        max(spec$ridge[rows])
      expect_match(VCMoE:::.vcmoe_dev_grid_block_reason(fit, 1L, VCMoE:::.vcmoe_dev_control(list())),
                   "component_mass_theta_penalty")
    }
  }
  expect_lt(gradients$kernel, 1e-3)
  expect_gt(gradients$component, 100 * gradients$kernel)
})

test_that("stored controls without the theta-mass key are read as legacy", {
  expect_identical(VCMoE:::.vcmoe_legacy_fit_control(list(negbin_theta_ridge = .05))$negbin_theta_penalty_mass,
                   "component")
  expect_identical(VCMoE:::.vcmoe_legacy_fit_control(list())$negbin_theta_penalty_mass, NULL)
  expect_identical(VCMoE:::.vcmoe_theta_penalty_mass(VCMoE:::.vcmoe_default_control(list())), "kernel")
  expect_error(VCMoE:::.vcmoe_theta_penalty_mass(list(negbin_theta_penalty_mass = "total")), "kernel")
  expect_identical(VCMoE:::.bootstrap_refit_control(list(negbin_theta_ridge = .05), list(), 1, 1)$negbin_theta_penalty_mass,
                   "component")
  expect_equal(VCMoE:::.fit_negbin_theta(c(1, 4, 9, 0), rep(3, 4), rep(.5, 4), 8,
                 VCMoE:::.vcmoe_default_control(list())),
               VCMoE:::.fit_negbin_theta(c(1, 4, 9, 0), rep(3, 4), rep(.5, 4), 8,
                 VCMoE:::.vcmoe_default_control(list()), penalty_mass = 100))
})

test_that("wild path draws are seeded, nested across path sets, and near normal for many clusters", {
  set.seed(9207)
  G <- 300L
  score <- matrix(stats::rnorm(G * 2L), G, 2L, dimnames = list(NULL, c("a", "b")))
  hessian <- array(0, c(G, 2L, 2L))
  for (g in seq_len(G)) hessian[g, , ] <- diag(2) / G
  omega <- lapply(seq_len(G), function(g) tcrossprod(score[g, ]) + diag(1e-6, 2))
  adj <- VCMoE:::.vcmoe_cluster_adjust(score, hessian, omega, diag(2) * 1.01, "CR2")
  adjustments <- list(adj, adj, adj)
  active <- rep(list(seq_len(G)), 3L)
  sets <- list(full = 1:3, inner = 2L)
  keys <- list(`1,2,3` = c("a", "b"), `2` = c("a", "b"))
  run <- function() VCMoE:::.vcmoe_wild_path_critical(adjustments, active, keys, sets, .95, 999L, 11L, G)
  set.seed(9208)
  state <- .Random.seed
  first <- run()
  expect_identical(.Random.seed, state)
  expect_identical(first, run())
  expect_true(all(first$full >= first$inner))
  expect_equal(unname(first$inner), rep(qnorm(.975), 2), tolerance = .1)
  stat <- VCMoE:::.vcmoe_wild_tstat(adj, VCMoE:::.vcmoe_webb_weights(5L, G))
  expect_equal(dim(stat), c(5L, 2L))
  expect_identical(colnames(stat), c("a", "b"))
})

test_that("public CR2 wild bands report finite df, t critical values and interior paths", {
  sim <- simulate_vcmoe_gaussian(n = 120, k = 2, seed = 9209, separation = 2,
                                 u = rep(seq(0, 1, length.out = 20L), each = 6L))
  fit <- suppressWarnings(vcmoe_fit(y ~ 1 | x1, data = sim$data, u = "u", k = 2,
    bandwidth = .5, u_grid = c(.3, .5, .7),
    control = list(maxit = 20L, n_starts = 1L, seed = 9210, warn_ambiguous = FALSE)))
  ids <- rep(1:20, each = 6L)
  ctl <- list(multiplier_B = 199L, seed = 9211, max_hessian_dim = 60L,
              covariance_adjustment = "CR2", simultaneous_method = "wild_cluster_path")
  band <- suppressWarnings(vcmoe_confband(fit, cluster = ids, type = "simultaneous", strict = FALSE, control = ctl))
  ok <- band$intervals$pointwise_status == "ok"
  expect_true(any(ok))
  expect_identical(band$settings$covariance_adjustment, "CR2")
  expect_identical(band$settings$simultaneous_method, "wild_cluster_path")
  default <- suppressWarnings(vcmoe_confband(fit, cluster = ids, type = "simultaneous", strict = FALSE,
    control = ctl[c("multiplier_B", "seed", "max_hessian_dim")]))
  expect_identical(default$intervals, band$intervals)
  expect_true(all(band$intervals$df[ok] > 0 & band$intervals$df[ok] <= band$intervals$active_clusters[ok]))
  expect_equal(band$intervals$pointwise_critical[ok], qt(.975, band$intervals$df[ok]))
  expect_true(all(is.finite(band$diagnostics$max_cluster_leverage)))
  expect_error(vcmoe_confband(fit, control = list(covariance_adjustment = "CR2")), "explicit `cluster`")
  all <- suppressWarnings(VCMoE:::.vcmoe_confband_comparison(fit, ids, strict = FALSE,
    control = list(multiplier_B = 99L, seed = 9212, max_hessian_dim = 60L, path_sets = list(interior = 2L)),
    small_sample = c("CR1", "CR2", "CR3")))
  expect_named(all, c("hc0_analytic", "hc0_multiplier", "cluster_multiplier", "cr1_wild", "cr2_wild", "cr3_wild"))
  cr2 <- all$cr2_wild$intervals
  inner <- cr2$simultaneous_status_interior == "ok"
  expect_true(any(inner))
  expect_true(all(cr2$grid_id[inner] == 2L))
  both <- intersect(cr2$parameter[inner], cr2$parameter[cr2$simultaneous_status == "ok"])
  for (key in both) {
    rows <- cr2$parameter == key
    expect_lte(max(cr2$simultaneous_critical_interior[rows], na.rm = TRUE),
               max(cr2$simultaneous_critical[rows]) + 1e-12)
  }
  expect_true("simultaneous_critical_interior" %in% names(all$cluster_multiplier$intervals))
  expect_equal(all$cr1_wild$intervals$se^2,
               all$cluster_multiplier$intervals$se^2 * all$cr1_wild$intervals$covariance_factor,
               tolerance = 1e-6)
  expect_error(VCMoE:::.vcmoe_confband_comparison(fit, ids, small_sample = "CR4"), "small_sample")
})

test_that("path tests leave bands unchanged and the zero test inverts the band", {
  sim <- simulate_vcmoe_gaussian(n = 120, k = 2, seed = 9209, separation = 2,
                                 u = rep(seq(0, 1, length.out = 20L), each = 6L))
  fit <- suppressWarnings(vcmoe_fit(y ~ z1 | 1, data = sim$data, u = "u", k = 2,
    bandwidth = .5, u_grid = c(.3, .4, .5, .6, .7),
    control = list(maxit = 20L, n_starts = 1L, seed = 9210, warn_ambiguous = FALSE)))
  ids <- rep(1:20, each = 6L)
  ctl <- list(multiplier_B = 199L, seed = 9213, max_hessian_dim = 60L, path_sets = list(interior = 2:4))
  plain <- suppressWarnings(VCMoE:::.vcmoe_confband_comparison(fit, ids, strict = FALSE, control = ctl,
                                                               small_sample = "CR2"))
  tested <- suppressWarnings(VCMoE:::.vcmoe_confband_comparison(fit, ids, strict = FALSE,
    control = c(ctl, path_tests = TRUE), small_sample = "CR2"))
  expect_null(plain$cr2_wild$tests)
  for (method in c("hc0_multiplier", "cluster_multiplier", "cr2_wild")) {
    expect_identical(tested[[method]]$intervals, plain[[method]]$intervals)
    tests <- tested[[method]]$tests
    expect_setequal(unique(tests$path_set), c("full", "interior"))
    expect_true(all(tests$zero_p_value > 0 & tests$zero_p_value <= 1))
    expect_true(all(tests$constancy_p_value > 0 & tests$constancy_p_value <= 1))
    intervals <- tested[[method]]$intervals
    for (set in c("full", "interior")) {
      suffix <- if (set == "full") "" else "_interior"
      ok <- intervals[[paste0("simultaneous_status", suffix)]] == "ok"
      excludes <- tapply(intervals[[paste0("simultaneous_lower", suffix)]][ok] > 0 |
                           intervals[[paste0("simultaneous_upper", suffix)]][ok] < 0,
                         intervals$parameter[ok], any)
      x <- tests[tests$path_set == set, ]
      far <- abs(x$zero_p_value - .05) > 2 / 200
      expect_identical((x$zero_p_value <= .05)[far], as.vector(excludes[x$parameter])[far])
    }
  }
  band <- suppressWarnings(vcmoe_confband(fit, cluster = ids, type = "simultaneous", strict = FALSE,
    coefficient_set = "expert", control = c(ctl, path_tests = TRUE)))
  expect_true(all(band$tests$coefficient_set == "expert"))
  expect_identical(band$tests, tested$cr2_wild$tests[tested$cr2_wild$tests$coefficient_set == "expert", ],
                   ignore_attr = TRUE)
})

test_that("constancy statistics vanish for constant paths and match their definition", {
  influence <- matrix(c(1, 2, -1, 0, 1, 1, 2, -2, 1), 3, 3)
  flat <- VCMoE:::.vcmoe_path_observed(influence, rep(.7, 3))
  expect_equal(flat$constant, 0)
  est <- c(.2, -.1, .5)
  out <- VCMoE:::.vcmoe_path_observed(influence, est)
  centered <- influence - rowMeans(influence)
  expect_equal(out$constant, max(abs(est - mean(est)) / sqrt(colSums(centered^2))))
  expect_equal(out$zero, max(abs(est) / sqrt(colSums(influence^2))))
  expect_true(is.na(VCMoE:::.vcmoe_path_observed(influence[, 1, drop = FALSE], est[1])$constant))
  expect_equal(VCMoE:::.vcmoe_path_p_value(c(1, 2, 3), 2), 3 / 4)
})

test_that("wild constancy draws match an explicit residual-score recomputation", {
  set.seed(9214)
  G <- 7L
  d <- 2L
  B <- 6L
  grids <- 3L
  active <- list(1:5, 2:7, c(1L, 3:7))
  adjustments <- lapply(seq_len(grids), function(j) {
    m <- length(active[[j]])
    score <- matrix(stats::rnorm(m * d), m, d, dimnames = list(NULL, c("a", "b")))
    hessian <- array(0, c(m, d, d))
    for (g in seq_len(m)) hessian[g, , ] <- crossprod(matrix(stats::rnorm(4 * d), 4, d)) / 10
    omega <- lapply(seq_len(m), function(g) tcrossprod(score[g, ]) + diag(.05, d))
    VCMoE:::.vcmoe_cluster_adjust(score, hessian, omega, Reduce(`+`, lapply(seq_len(m), function(g) hessian[g, , ])) + diag(.1, d), "CR2")
  })
  weights <- VCMoE:::.vcmoe_webb_weights(B, G)
  sets <- list(full = list(grids = 1:3, keys = c("a", "b")), pair = list(grids = 2:3, keys = "b"))
  got <- VCMoE:::.vcmoe_wild_constancy_draws(adjustments, active, weights, sets)
  # Rebuild from residual bootstrap scores U*_g - P_g sum_h U*_h, mapped back by R_g.
  expected <- lapply(sets, function(set) {
    t(vapply(seq_len(B), function(b) {
      numerator <- rows <- list()
      for (j in set$grids) {
        a <- adjustments[[j]]
        m <- length(active[[j]])
        P <- lapply(seq_len(m), function(i) solve(a$R[[i]], a$K[[i]]))
        star <- weights[b, active[[j]]] * a$score
        total <- colSums(star)
        numerator[[as.character(j)]] <- as.numeric(a$inv %*% total)
        r <- matrix(0, G, d, dimnames = list(NULL, colnames(a$score)))
        for (i in seq_len(m)) r[active[[j]][i], ] <- a$R[[i]] %*% (star[i, ] - P[[i]] %*% total)
        rows[[as.character(j)]] <- r
      }
      center_n <- Reduce(`+`, numerator) / length(numerator)
      center_r <- Reduce(`+`, rows) / length(rows)
      stat <- vapply(names(rows), function(j) {
        abs(numerator[[j]] - center_n) / sqrt(colSums((rows[[j]] - center_r)^2))
      }, numeric(d))
      apply(stat, 1L, max)[match(set$keys, c("a", "b"))]
    }, numeric(length(set$keys))))
  })
  expect_equal(unname(as.matrix(got[[1]])), unname(expected$full), tolerance = 1e-10)
  expect_equal(unname(as.matrix(got[[2]])), unname(t(expected$pair)), tolerance = 1e-10)
})

test_that("cell Hessian columns sum to the cluster Hessian slices", {
  sim <- simulate_vcmoe_negbin(n = 120, k = 2, seed = 9203, separation = 1.2,
                               u = rep(seq(0, 1, length.out = 20L), each = 6L), mean_count = 8)
  fit <- suppressWarnings(vcmoe_fit(y ~ z1 + offset(log_size_factor) | x1, data = sim$data, u = "u",
    k = 2, family = "negative-binomial", bandwidth = .6, u_grid = .5,
    control = list(maxit = 20L, n_starts = 1L, seed = 9204, warn_ambiguous = FALSE,
                   negbin_theta_ridge = .05, negbin_theta_target = 8)))
  group <- rep(1:20, each = 6L)
  spec <- VCMoE:::.vcmoe_dev_parameter_spec(fit, 1L)
  par <- VCMoE:::.vcmoe_dev_parameter_vector(fit, 1L)
  table <- spec$parameter_table
  effect <- which(table$coefficient_set %in% c("expert", "gating_contrast") & table$block == "intercept" &
                    table$term == "(Intercept)")
  expect_identical(names(par)[effect], c("expert:component1:(Intercept):intercept",
    "expert:component2:(Intercept):intercept", "gating_contrast:component1_vs_component2:(Intercept):intercept"))
  cells <- VCMoE:::.vcmoe_dev_cell_hessian_columns(par, spec, effect)
  blocks <- VCMoE:::.vcmoe_dev_cluster_hessian(par, spec, group)
  for (a in seq_along(effect)) {
    expect_equal(unname(rowsum(cells[, , a], group)), unname(blocks[, , effect[a]]), tolerance = 1e-6)
  }
})

test_that("Webb weights are enumerated when 6^G does not exceed B", {
  w <- VCMoE:::.vcmoe_webb_weights(999L, 3L)
  expect_true(attr(w, "enumerated"))
  expect_identical(dim(w), c(216L, 3L))
  expect_identical(nrow(unique(w)), 216L)
  expect_setequal(unique(as.numeric(w)), VCMoE:::.vcmoe_webb_values)
  s <- VCMoE:::.vcmoe_webb_weights(999L, 5L)
  expect_false(attr(s, "enumerated"))
  expect_identical(dim(s), c(999L, 5L))
})

test_that("few clusters give an enumerated wild band and a sensitivity warning", {
  sim <- simulate_vcmoe_negbin(n = 160, k = 2, seed = 9301, separation = 1.3,
                               u = rep(seq(0, 1, length.out = 4L), each = 40L))
  fit <- suppressWarnings(vcmoe_fit(y ~ z1 + offset(log_size_factor) | x1, data = sim$data, u = "u",
    k = 2, family = "negative-binomial", bandwidth = .8, u_grid = c(.4, .6),
    control = list(maxit = 300L, n_starts = 2L, seed = 9302, warn_ambiguous = FALSE)))
  ids <- rep(1:4, each = 40L)
  expect_warning(band <- vcmoe_confband(fit, cluster = ids, type = "simultaneous", strict = FALSE,
                                        control = list(multiplier_B = 9999L, seed = 1)),
                 "Few locally supported")
  again <- suppressWarnings(vcmoe_confband(fit, cluster = ids, type = "simultaneous", strict = FALSE,
                                           control = list(multiplier_B = 9999L, seed = 2)))
  # Enumeration makes the band independent of the seed.
  ok <- band$intervals$status == "ok"
  expect_true(any(ok))
  expect_equal(band$intervals$simultaneous_critical[ok], again$intervals$simultaneous_critical[ok])
  expect_true(all(band$diagnostics$few_clusters))
})
