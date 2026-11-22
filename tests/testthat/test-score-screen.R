make_screen_fixture <- function(grid_indices, n_starts) {
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    set.seed(1)
    G <- 12L
    m <- 80L
    donor <- rep(seq_len(G), each = m)
    u <- ((seq_len(G) - 1) / (G - 1))[donor]
    n <- length(donor)
    marker <- stats::rnorm(n)
    comp <- ifelse(stats::rbinom(n, 1, stats::plogis(-0.3 + 1.8 * marker + 0.8 * u)) == 1, 1L, 2L)
    sig <- stats::rnorm(n)
    lib <- exp(stats::rnorm(n, 0, 0.3))
    eta <- log(lib) + ifelse(comp == 1, 1 + 0.3 * u + 0.4 * sig, 0.2 - 0.5 * u - 0.2 * sig) +
      stats::rnorm(G, 0, 0.1)[donor]
    dat <- data.frame(y = stats::rnbinom(n, size = 6, mu = exp(eta)), lo = log(lib), u = u, donor = donor,
                      marker = marker, sig = sig, mk2 = stats::rnorm(n) + 1.2 * (comp == 1),
                      noise1 = stats::rnorm(n), noise2 = stats::rnorm(n))
    fit <- suppressWarnings(vcmoe_fit(y ~ offset(lo) | marker, data = dat, u = "u", k = 2L,
      family = "negative-binomial", bandwidth = 0.3, u_grid = sort(unique(dat$u))[grid_indices],
      control = list(seed = 3, n_starts = n_starts, maxit = 300, tol = 1e-7, keep_data = TRUE, warn_ambiguous = FALSE,
                     negbin_theta_ridge = 0.05, negbin_theta_target = 8, negbin_ridge = 1e-4,
                     negbin_theta_penalty_mass = "kernel"),
      label = "align", parameterization = "a1_epanechnikov_scaled", u_scale = "unit"))
    X <- as.matrix(dat[c("sig", "mk2", "noise1", "noise2")])
    cache <<- list(fit = fit, dat = dat, X = sweep(X, 2L, colMeans(X)))
    cache
  }
}

screen_toy <- make_screen_fixture(seq_len(12L), 2L)
screen_two_grids <- make_screen_fixture(c(4L, 8L), 1L)

test_that("closed-form screen blocks match finite differences of the augmented local likelihood", {
  skip_on_cran()
  toy <- screen_toy()
  group <- VCMoE:::.vcmoe_cluster_index(toy$fit, toy$dat$donor)
  ctl <- VCMoE:::.vcmoe_dev_control(list(covariance_adjustment = "HC0"))
  rel <- function(a, b) max(abs(a - b)) / max(abs(b))
  # Expert test: the gene enters experts and gating; gating diagnostic: gating only.
  for (case in list(list(families = c("expert", "gating"), tested = "expert"),
                    list(families = "gating", tested = "gating"))) {
    args <- list(toy$fit, 4L, toy$X, group, case$families, case$tested, ctl)
    mom <- do.call(VCMoE:::.vcmoe_screen_grid_moments, c(args, information = "observed"))
    complete <- do.call(VCMoE:::.vcmoe_screen_grid_moments, args)
    expect_identical(mom$status, "ok")
    expect_equal(complete$U, mom$U)
    expect_equal(complete$Hx, mom$Hx)
    for (j in 1:2) {
      blocks <- VCMoE:::.vcmoe_screen_gene_blocks(mom, j)
      aug <- VCMoE:::.vcmoe_screen_augment_fit(toy$fit, toy$X[, j], colnames(toy$X)[j], case$families)
      spec <- VCMoE:::.vcmoe_dev_parameter_spec(aug, 4L)
      par <- VCMoE:::.vcmoe_dev_parameter_vector(aug, 4L)
      score <- VCMoE:::.vcmoe_dev_score_matrix(par, spec, ctl$finite_diff_eps)
      hessian <- VCMoE:::.vcmoe_dev_cluster_hessian(par, spec, group, ctl$hessian_eps)
      ord <- match(c(names(mom$par), VCMoE:::.vcmoe_screen_new_names(mom, colnames(toy$X)[j])), names(par))
      expect_false(anyNA(ord))
      expect_identical(length(ord), length(par))
      expect_lt(rel(blocks$score, rowsum(score, group, reorder = FALSE)[mom$clusters, ord]), 1e-8)
      expect_lt(rel(blocks$hessian, hessian[mom$clusters, ord, ord]), 1e-5)
      expect_lt(rel(blocks$bread, (colSums(hessian) + diag(spec$ridge))[ord, ord]), 1e-5)
      omega <- lapply(mom$clusters, function(g) crossprod(score[group == g, ord]))
      expect_lt(max(mapply(rel, blocks$omega, omega)), 1e-8)
      # Complete-data information: the efficient information of the tested
      # coefficients given everything else is J^c + ridge; all other blocks
      # are unchanged.
      working <- VCMoE:::.vcmoe_screen_gene_blocks(complete, j)
      tested <- complete$d + which(complete$tested)
      rest <- setdiff(seq_len(ncol(working$bread)), tested)
      expect_equal(working$bread[rest, ], blocks$bread[rest, ])
      expect_equal(working$score, blocks$score)
      schur <- working$bread[tested, tested] -
        working$bread[tested, rest] %*% solve(working$bread[rest, rest], working$bread[rest, tested])
      jc <- matrix(colSums(matrix(complete$Hc[, j, ], nrow = length(complete$clusters))), complete$e, complete$e,
                   byrow = TRUE)[which(complete$tested), which(complete$tested)]
      expect_equal(schur, (jc + t(jc)) / 2 + diag(complete$ridge_new[complete$tested]), tolerance = 1e-8,
                   ignore_attr = TRUE)
    }
  }
})

test_that("one-step scores solve the estimating equation; restricted t is capped by the active clusters", {
  skip_on_cran()
  toy <- screen_toy()
  group <- VCMoE:::.vcmoe_cluster_index(toy$fit, toy$dat$donor)
  ctl <- VCMoE:::.vcmoe_dev_control(list(covariance_adjustment = "HC0"))
  mom <- VCMoE:::.vcmoe_screen_grid_moments(toy$fit, 4L, toy$X, group, c("expert", "gating"), "expert", ctl)
  blocks <- VCMoE:::.vcmoe_screen_gene_blocks(mom, 1L)
  total <- colSums(blocks$score) + c(mom$penalty, numeric(mom$e))
  delta <- solve(blocks$bread, total)
  fitted <- blocks$score - t(vapply(seq_len(nrow(blocks$score)), function(i) {
    as.numeric(blocks$hessian[i, , ] %*% delta)
  }, numeric(ncol(blocks$score))))
  ridge <- diag(blocks$bread - Reduce(`+`, lapply(seq_len(dim(blocks$hessian)[1L]), function(i) blocks$hessian[i, , ])))
  # Linearized scores solve the penalized estimating equation at the one-step point.
  expect_equal(unname(colSums(fitted) + c(mom$penalty, numeric(mom$e)) - ridge * delta), numeric(length(delta)),
               tolerance = 1e-8)
  # CR2 rescales the cluster rows, so allow its adjustment on top of the CR0 bound sqrt(G).
  bound <- 1.5 * sqrt(length(mom$clusters))
  for (j in seq_len(ncol(toy$X))) {
    restricted <- VCMoE:::.vcmoe_screen_gene_grid(mom, j, studentize = "restricted")
    expect_true(all(abs(restricted$estimate / restricted$se) <= bound))
  }
})

test_that("the screen is chunk-invariant, reproducible and ranks the within-component signal first", {
  toy <- screen_two_grids()
  grids <- seq_along(toy$fit$u_grid)
  full <- vcmoe_gene_screen(toy$fit, toy$X, toy$dat$donor, grids = grids, B = 199L, seed = 7)
  part <- VCMoE:::.vcmoe_score_screen(toy$fit, toy$X[, c("noise2", "sig")], toy$dat$donor, grids = grids,
                                      B = 199L, seed = 7, center = FALSE)
  expect_equal(part$draws[, "sig"], full$draws[, "sig"])
  expect_equal(part$genes$statistic[part$genes$gene == "sig"], full$genes$statistic[full$genes$gene == "sig"])
  expect_identical(full$genes$gene[which.min(full$genes$p_value)], "sig")
  expect_lte(full$genes$p_value[full$genes$gene == "sig"], 0.05)
  expect_true(all(full$genes$p_value > 0 & full$genes$p_value <= 1))
  expect_identical(nrow(full$pointwise), ncol(toy$X) * length(grids) * 2L * 2L)
  adjusted <- VCMoE:::.vcmoe_screen_minp(full$genes$statistic, full$draws)
  expect_true(all(adjusted >= full$genes$p_value))
  # A gene already in the null gating is tested without a second gating column.
  with_gate <- cbind(toy$X[, c("sig", "noise1")], marker = toy$dat$marker)
  mixed <- VCMoE:::.vcmoe_score_screen(toy$fit, with_gate, toy$dat$donor, grids = grids, B = 199L, seed = 7)
  expect_identical(mixed$genes$gene, colnames(with_gate))
  expect_identical(mixed$genes$adjust_gating, c(TRUE, TRUE, FALSE))
  expect_true(all(is.finite(mixed$genes$statistic)))
  expect_equal(mixed$draws[, "sig"], full$draws[, "sig"])
  gating <- VCMoE:::.vcmoe_score_screen(toy$fit, toy$X[, 1:2], toy$dat$donor, grids = grids, side = "gating", B = 0L)
  expect_null(gating$draws)
  expect_identical(gating$components, "component1_vs_component2")
})

test_that("screen blocks match finite differences when the null experts carry covariates (covariate null)", {
  skip_on_cran()
  set.seed(11)
  G <- 10L
  m <- 70L
  donor <- rep(seq_len(G), each = m)
  u <- ((seq_len(G) - 1) / (G - 1))[donor]
  n <- length(donor)
  marker <- stats::rnorm(n)
  comp <- ifelse(stats::rbinom(n, 1, stats::plogis(-0.2 + 1.8 * marker)) == 1, 1L, 2L)
  state <- stats::rnorm(n)
  lib <- exp(stats::rnorm(n, 0, 0.3))
  eta <- log(lib) + ifelse(comp == 1, 1 + 0.3 * u + 0.4 * state, 0.2 - 0.4 * u - 0.3 * state)
  dat <- data.frame(y = stats::rnbinom(n, size = 6, mu = exp(eta)), lo = log(lib), u = u, donor = donor,
                    marker = marker, statepc1 = state, g1 = state + stats::rnorm(n), g2 = stats::rnorm(n))
  fit <- suppressWarnings(vcmoe_fit(y ~ offset(lo) + statepc1 | marker, data = dat, u = "u", k = 2L,
    family = "negative-binomial", bandwidth = 0.35, u_grid = sort(unique(dat$u)),
    control = list(seed = 3, n_starts = 2, maxit = 300, tol = 1e-7, keep_data = TRUE, warn_ambiguous = FALSE,
                   negbin_theta_ridge = 0.05, negbin_theta_target = 8, negbin_ridge = 1e-4,
                   negbin_theta_penalty_mass = "kernel"),
    label = "align", parameterization = "a1_epanechnikov_scaled", u_scale = "unit"))
  expect_identical(dimnames(fit$coefficients$expert)[[3L]], c("(Intercept)", "statepc1"))
  X <- as.matrix(dat[c("g1", "g2")])
  X <- sweep(X, 2L, colMeans(X))
  group <- VCMoE:::.vcmoe_cluster_index(fit, dat$donor)
  ctl <- VCMoE:::.vcmoe_dev_control(list(covariance_adjustment = "HC0"))
  rel <- function(a, b) max(abs(a - b)) / max(abs(b))
  mom <- VCMoE:::.vcmoe_screen_grid_moments(fit, 4L, X, group, c("expert", "gating"), "expert", ctl,
                                            information = "observed")
  expect_identical(mom$status, "ok")
  for (j in 1:2) {
    blocks <- VCMoE:::.vcmoe_screen_gene_blocks(mom, j)
    aug <- VCMoE:::.vcmoe_screen_augment_fit(fit, X[, j], colnames(X)[j], c("expert", "gating"))
    spec <- VCMoE:::.vcmoe_dev_parameter_spec(aug, 4L)
    par <- VCMoE:::.vcmoe_dev_parameter_vector(aug, 4L)
    score <- VCMoE:::.vcmoe_dev_score_matrix(par, spec, ctl$finite_diff_eps)
    hessian <- VCMoE:::.vcmoe_dev_cluster_hessian(par, spec, group, ctl$hessian_eps)
    ord <- match(c(names(mom$par), VCMoE:::.vcmoe_screen_new_names(mom, colnames(X)[j])), names(par))
    expect_false(anyNA(ord))
    expect_identical(length(ord), length(par))
    expect_lt(rel(blocks$score, rowsum(score, group, reorder = FALSE)[mom$clusters, ord]), 1e-8)
    expect_lt(rel(blocks$hessian, hessian[mom$clusters, ord, ord]), 1e-5)
    expect_lt(rel(blocks$bread, (colSums(hessian) + diag(spec$ridge))[ord, ord]), 1e-5)
  }
  res <- VCMoE:::.vcmoe_score_screen(fit, X, dat$donor, grids = 2:9, B = 199L, seed = 5)
  expect_true(all(is.finite(res$genes$statistic)))
  expect_true(all(res$genes$p_value > 0 & res$genes$p_value <= 1))
})

test_that("the gating-profiled screen refits the gene's gating under H0 and screens on that fit", {
  skip_on_cran()
  toy <- screen_toy()
  grids <- 3:9
  # No EM iteration: the restricted fit is the null with a zero gating column.
  still <- VCMoE:::.vcmoe_screen_profile_fit(toy$fit, toy$X[, "mk2"], "mk2", grids, maxit = 0L)
  aug <- VCMoE:::.vcmoe_screen_augment_fit(toy$fit, toy$X[, "mk2"], "mk2", sides = "gating")
  for (g in grids) {
    expect_equal(VCMoE:::.vcmoe_dev_parameter_vector(still$fit, g), VCMoE:::.vcmoe_dev_parameter_vector(aug, g))
  }
  expect_equal(still$grids$gating_contrast, rep(0, length(grids)))
  # EM ascends the penalized local objective; mk2 marks the same generating
  # component as the gating covariate `marker`, so its fitted gating contrast
  # takes the sign of marker's (the fitted labels may be either way round).
  prof <- VCMoE:::.vcmoe_screen_profile_fit(toy$fit, toy$X[, "mk2"], "mk2", grids)
  expect_true(all(prof$grids$objective_gain > -1e-8))
  expect_true(all(prof$grids$converged))
  marker_sign <- sign(stats::median(toy$fit$coefficients$gating[grids, 1L, "marker"] -
                                      toy$fit$coefficients$gating[grids, 2L, "marker"]))
  expect_gt(marker_sign * stats::median(prof$grids$gating_contrast), 0.3)
  expect_gt(min(prof$grids$map_agreement), 0.8)
  # Grids outside the profiled set keep the null parameters.
  expect_equal(VCMoE:::.vcmoe_dev_parameter_vector(prof$fit, 10L), VCMoE:::.vcmoe_dev_parameter_vector(aug, 10L))
  with_gate <- cbind(toy$X[, c("sig", "mk2")], marker = toy$dat$marker)
  res <- VCMoE:::.vcmoe_score_screen_profiled(toy$fit, with_gate, toy$dat$donor, grids = grids, B = 199L, seed = 7)
  expect_identical(res$genes$gene, colnames(with_gate))
  expect_identical(res$genes$profiled, c(TRUE, TRUE, FALSE))
  expect_true(all(is.finite(res$genes$statistic)))
  expect_identical(dim(res$draws), c(199L, 3L))
  expect_identical(sort(unique(res$profile$gene)), c("mk2", "sig"))
  # The gene already in the null gating is screened on the null fit unchanged.
  plain <- VCMoE:::.vcmoe_score_screen(toy$fit, with_gate[, "marker", drop = FALSE], toy$dat$donor, grids = grids,
                                       B = 199L, seed = 7)
  expect_equal(res$draws[, "marker"], plain$draws[, "marker"])
})

test_that("finite differences on the kernel support equal those on every row", {
  skip_on_cran()
  toy <- screen_toy()
  group <- VCMoE:::.vcmoe_cluster_index(toy$fit, toy$dat$donor)
  ctl <- list(covariance_adjustment = "HC0")
  # Grid 3 of 12 with h = 0.3: part of the cells have zero kernel weight.
  spec <- VCMoE:::.vcmoe_dev_parameter_spec(toy$fit, 3L)
  expect_true(any(spec$weights == 0) && any(spec$weights > 0))
  fast <- VCMoE:::.vcmoe_dev_score_hessian(toy$fit, grid_id = 3L, control = ctl, group = group)
  all_rows <- function(expr) {
    old <- options(vcmoe.dev_all_rows = TRUE)
    on.exit(options(old))
    expr
  }
  slow <- all_rows(VCMoE:::.vcmoe_dev_score_hessian(toy$fit, grid_id = 3L, control = ctl, group = group))
  expect_identical(fast$status, "ok")
  expect_equal(fast$score, slow$score, tolerance = 1e-12)
  expect_equal(fast$cluster_hessian, slow$cluster_hessian, tolerance = 1e-12)
  expect_equal(fast$bread, slow$bread, tolerance = 1e-12)
  expect_true(all(fast$score[spec$weights == 0, ] == 0))
  a <- VCMoE:::.vcmoe_score_screen(toy$fit, toy$X, toy$dat$donor, grids = 2:5, B = 199L, seed = 3)
  b <- all_rows(VCMoE:::.vcmoe_score_screen(toy$fit, toy$X, toy$dat$donor, grids = 2:5, B = 199L, seed = 3))
  expect_equal(a$genes, b$genes, tolerance = 1e-10)
  expect_equal(a$draws, b$draws, tolerance = 1e-10)
})

test_that("vcmoe_gene_screen() shares one set of weights across chunks, profiles and seed = NULL", {
  skip_on_cran()
  toy <- screen_toy()
  set.seed(42)
  a <- vcmoe_gene_screen(toy$fit, toy$X, toy$dat$donor, grids = 2:11, B = 199L, seed = NULL, chunk_size = 1L,
                         adjust = "both")
  set.seed(42)
  b <- vcmoe_gene_screen(toy$fit, toy$X, toy$dat$donor, grids = 2:11, B = 199L, seed = NULL, chunk_size = 10L,
                         adjust = "both")
  expect_s3_class(a, "vcmoe_gene_screen")
  expect_equal(a$genes, b$genes)
  expect_equal(a$draws, b$draws)
  expect_identical(a$genes$gene, colnames(toy$X))
  expect_identical(a$genes$gene[which.max(a$genes$statistic)], "sig")
  expect_true(all(c("p_value", "q_value", "p_minp", "effect_per_sd", "status", "argmax_component") %in% names(a$genes)))
  expect_equal(a$genes$q_value, stats::p.adjust(a$genes$p_value, "BH"))
  expect_output(print(a), "candidates: 4")
  # Same seed: the plain screen and the internal screen agree.
  c1 <- vcmoe_gene_screen(toy$fit, toy$X, toy$dat$donor, grids = 2:11, B = 199L, seed = 7)
  internal <- VCMoE:::.vcmoe_score_screen(toy$fit, toy$X, toy$dat$donor, grids = 2:11, B = 199L, seed = 7)
  expect_equal(c1$genes$p_value, internal$genes$p_value)
  # A profiled candidate uses the same weights as the others.
  with_gate <- cbind(toy$X, marker = toy$dat$marker - mean(toy$dat$marker), mk = toy$X[, "mk2"])
  colnames(with_gate)[ncol(with_gate)] <- "mk_copy"
  prof <- vcmoe_gene_screen(toy$fit, with_gate[, c("sig", "mk_copy")], toy$dat$donor, profile_genes = "mk_copy",
                            grids = 2:11, B = 199L, seed = 7)
  expect_identical(prof$genes$profiled, c(FALSE, TRUE))
  expect_true(is.data.frame(prof$profile))
  expect_equal(prof$genes$p_value[1L], c1$genes$p_value[c1$genes$gene == "sig"])
})

test_that("vcmoe_gene_screen() keeps incomplete candidates in the BH family with NA p-values", {
  skip_on_cran()
  toy <- screen_toy()
  full <- vcmoe_gene_screen(toy$fit, toy$X, toy$dat$donor, grids = 2:11, B = 99L, seed = 3)
  original <- VCMoE:::.vcmoe_screen_gene_grid
  testthat::local_mocked_bindings(.package = "VCMoE",
    .vcmoe_screen_gene_grid = function(mom, j, weights = NULL, studentize = "onestep") {
      out <- original(mom, j, weights, studentize)
      # Block the first candidate (sig) at grid 5.
      if (j == 1L && mom$grid_id == 5L) {
        out$blocked <- TRUE
        out$estimate[] <- NA_real_
        out$se[] <- NA_real_
      }
      out
    })
  res <- vcmoe_gene_screen(toy$fit, toy$X, toy$dat$donor, grids = 2:11, B = 99L, seed = 3)
  g <- res$genes
  expect_identical(g$status, c("incomplete_path", "ok", "ok", "ok"))
  expect_true(is.na(g$p_value[1L]))
  expect_true(is.finite(g$p_value_raw_path[1L]))
  expect_identical(g$n_blocked[1L], 1L)
  # The other candidates are unchanged, and BH still counts four candidates.
  expect_equal(g$p_value[-1L], full$genes$p_value[-1L])
  expect_equal(g$q_value, stats::p.adjust(g$p_value, "BH", n = 4L))
  expect_true(all(g$q_value[-1L] >= stats::p.adjust(g$p_value[-1L], "BH")))
})

test_that("vcmoe_gene_screen() validates its inputs", {
  toy <- screen_two_grids()
  expect_error(vcmoe_gene_screen(toy$fit, toy$X, NULL), "`cluster` is required")
  expect_error(vcmoe_gene_screen(toy$fit, toy$X), "`cluster` is required")
  expect_error(vcmoe_gene_screen(toy$fit, toy$X[-1L, ], toy$dat$donor), "one row per retained")
  expect_error(vcmoe_gene_screen(toy$fit, toy$X, toy$dat$donor, profile_genes = "nope"), "not in `X`")
  expect_error(vcmoe_gene_screen(toy$fit, toy$X, toy$dat$donor, grids = 99L), "grid indices")
  expect_error(vcmoe_gene_screen(toy$fit, toy$X, toy$dat$donor, chunk_size = 0), "chunk_size")
  ids <- seq_len(nrow(toy$dat))
  expect_warning(expect_error(vcmoe_gene_screen(toy$fit, toy$X, ids, B = 60000L), "exceed"), "own cluster")
  gaussian <- toy$fit
  gaussian$family <- "gaussian"
  expect_error(vcmoe_gene_screen(gaussian, toy$X, toy$dat$donor), "Negative-Binomial")
})
