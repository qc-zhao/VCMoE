.vcmoe_cluster_require_fit <- function(fit) {
  .vcmoe_dev_require_fit(fit)
  if (!identical(fit$engine_id %||% "local_grid_em", "local_grid_em") ||
      isTRUE(fit$diagnostics$reduced_fit)) {
    stop("Cluster/path inference requires an unrestricted local-grid fit; joint-path and reduced fits are unsupported.", call. = FALSE)
  }
  invisible(TRUE)
}

.vcmoe_cluster_index <- function(fit, cluster) {
  n <- length(fit$fitted$u)
  if (is.null(cluster)) return(seq_len(n))
  if (!is.atomic(cluster) || !is.null(dim(cluster)) || length(cluster) != n ||
      anyNA(cluster) || any(!nzchar(trimws(as.character(cluster)))) ||
      (is.numeric(cluster) && any(!is.finite(cluster)))) {
    stop("`cluster` must contain one nonmissing ID per retained fitted row, in fitted-row order.", call. = FALSE)
  }
  if (!is.null(names(cluster))) {
    expected <- rownames(fit$fitted$z_design) %||% as.character(fit$rows_used)
    if (length(expected) != n || !identical(names(cluster), expected)) {
      stop("Named `cluster` IDs must match retained design-matrix row names in exactly the same order.", call. = FALSE)
    }
  }
  values <- as.character(cluster)
  group <- match(values, unique(values))
  if (max(group) < 2L) stop("Cluster inference needs at least two independent clusters.", call. = FALSE)
  group
}

.vcmoe_cluster_rank <- function(x) {
  if (is.null(x) || any(!is.finite(x))) return(NA_integer_)
  values <- eigen((x + t(x)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (max(abs(values)) == 0) return(0L)
  sum(values > max(abs(values)) * 1e-8)
}

.vcmoe_cluster_design <- function(weights, group) {
  totals <- as.numeric(rowsum(matrix(weights, ncol = 1L), group, reorder = FALSE))
  positive <- totals > 0
  list(n_clusters = max(group), active_clusters = sum(positive),
       effective_clusters = .vcmoe_dev_effective_n(totals),
       max_cluster_weight_share = if (sum(totals) > 0) max(totals) / sum(totals) else NA_real_)
}

.vcmoe_inference_moments <- function(fit, strict, control, group = NULL) {
  control <- .vcmoe_dev_control(utils::modifyList(control, list(covariance_adjustment = "HC0")))
  lapply(seq_along(fit$u_grid), function(j) {
    if (strict && nzchar(.vcmoe_dev_grid_block_reason(fit, j, control))) return(NULL)
    if (is.null(group)) return(.vcmoe_dev_score_hessian(fit, grid_id = j, control = control))
    .vcmoe_dev_score_hessian(fit, grid_id = j, control = control, group = group)
  })
}

.vcmoe_seeded <- function(seed, expr) {
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed < 0 || seed > .Machine$integer.max) stop("Invalid multiplier seed.", call. = FALSE)
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  expr
}

.vcmoe_check_multiplier_B <- function(B) {
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 2 || B != as.integer(B)) {
    stop("`multiplier_B` must be an integer of at least two.", call. = FALSE)
  }
  invisible(TRUE)
}

# Observed sup-t statistics of beta(u) = 0 and beta(u) = const over one path.
# The constancy contrast is beta_hat(u_j) minus its equal-weight path mean
# (an average-of-grid-estimates constant fit); `influence` is units x grids.
.vcmoe_path_observed <- function(influence, estimate) {
  se <- sqrt(colSums(influence^2))
  centered <- influence - rowMeans(influence)
  contrast_se <- sqrt(colSums(centered^2))
  constant_ok <- length(estimate) > 1L && all(is.finite(contrast_se) & contrast_se > 0)
  list(zero = max(abs(estimate) / se), contrast_se = contrast_se,
       constant = if (constant_ok) max(abs(estimate - mean(estimate)) / contrast_se) else NA_real_)
}

.vcmoe_path_p_value <- function(draws, observed) {
  if (!is.finite(observed)) return(NA_real_)
  (1 + sum(draws >= observed)) / (length(draws) + 1)
}

.vcmoe_path_test_frame <- function(keys, zero, zero_draws, constant, constant_draws) {
  data.frame(parameter = keys, zero_statistic = zero,
             zero_p_value = vapply(seq_along(keys), function(k) .vcmoe_path_p_value(zero_draws[, k], zero[k]), numeric(1L)),
             constancy_statistic = constant,
             constancy_p_value = vapply(seq_along(keys), function(k) .vcmoe_path_p_value(constant_draws[, k], constant[k]), numeric(1L)),
             stringsAsFactors = FALSE)
}

# All columns share one Gaussian multiplier vector. For many cells the joint
# Gram factor is distributionally identical and avoids a B-by-n allocation.
# With `estimate` (one value per influence column), the same draws also give
# the zero and constancy path tests, attached as attribute "tests".
.vcmoe_path_multiplier <- function(influence, keys, level, B, seed = NULL, estimate = NULL) {
  .vcmoe_check_multiplier_B(B)
  .vcmoe_seeded(seed, {
    if (!ncol(influence)) return(numeric())
    norm <- sqrt(colSums(influence^2))
    if (any(!is.finite(norm)) || any(norm <= 0)) stop("Invalid path influence scale.", call. = FALSE)
    scaled <- sweep(influence, 2L, norm, "/")
    if (nrow(scaled) > ncol(scaled)) {
      gram <- crossprod(scaled)
      eig <- eigen((gram + t(gram)) / 2, symmetric = TRUE)
      if (min(eig$values) < -1e-8 * max(abs(eig$values))) stop("Invalid path covariance.", call. = FALSE)
      # Remove numerical roundoff only; do not inflate low-rank covariance.
      keep <- eig$values > 0
      factor <- t(sweep(eig$vectors[, keep, drop = FALSE], 2L, sqrt(eig$values[keep]), "*"))
    } else factor <- scaled
    curve_keys <- unique(keys)
    columns <- lapply(curve_keys, function(key) which(keys == key))
    maxima <- matrix(NA_real_, nrow = B, ncol = length(columns))
    testing <- !is.null(estimate)
    if (testing) {
      observed <- lapply(columns, function(cols) .vcmoe_path_observed(influence[, cols, drop = FALSE], estimate[cols]))
      constant_maxima <- matrix(NA_real_, nrow = B, ncol = length(columns))
    }
    for (start in seq.int(1L, B, by = 128L)) {
      rows <- start:min(B, start + 127L)
      draws <- matrix(stats::rnorm(length(rows) * nrow(factor)), nrow = length(rows)) %*% factor
      for (k in seq_along(columns)) {
        maxima[rows, k] <- apply(abs(draws[, columns[[k]], drop = FALSE]), 1L, max)
        if (testing && is.finite(observed[[k]]$constant)) {
          raw <- sweep(draws[, columns[[k]], drop = FALSE], 2L, norm[columns[[k]]], "*")
          contrast <- sweep(raw - rowMeans(raw), 2L, observed[[k]]$contrast_se, "/")
          constant_maxima[rows, k] <- apply(abs(contrast), 1L, max)
        }
      }
    }
    out <- stats::setNames(vapply(seq_along(columns), function(k) {
      as.numeric(stats::quantile(maxima[, k], level, names = FALSE))
    }, numeric(1L)), curve_keys)
    if (testing) {
      attr(out, "tests") <- .vcmoe_path_test_frame(curve_keys,
        vapply(observed, `[[`, numeric(1L), "zero"), maxima,
        vapply(observed, `[[`, numeric(1L), "constant"), constant_maxima)
    }
    out
  })
}

# Pseudo-power of a symmetric PSD matrix. Directions below `tol` relative to
# the largest eigenvalue are dropped, as in generalized-inverse CR2.
.vcmoe_psd_power <- function(x, power, tol = 1e-8) {
  eig <- eigen((x + t(x)) / 2, symmetric = TRUE)
  keep <- eig$values > tol * max(abs(eig$values), 0)
  vectors <- eig$vectors[, keep, drop = FALSE]
  structure(vectors %*% (eig$values[keep]^power * t(vectors)),
            min_value = if (any(keep)) min(eig$values[keep]) else NA_real_,
            dropped = sum(!keep))
}

# Score-space small-cluster adjustment for an M-estimator with penalized
# bread A = sum_g J_g + ridge. Fitted cluster scores satisfy, to first order,
# U_hat_g = U_g - P_g sum_h U_h with P_g = J_g A^{-1}. Adjusted scores
# M_g U_hat_g give the variance rows A^{-1} M_g U_hat_g.
# CR1: M_g = sqrt(G/(G-1)) I, t(G-1).
# CR2: M_g makes E[M_g U_hat_g U_hat_g' M_g'] = Omega_g under the working
#      independent-cell covariance Omega_g (Bell-McCaffrey/Pustejovsky-Tipton;
#      equal to classical CR2 for unweighted linear least squares).
# CR3: M_g = (I - P_g)^{-1}, the one-step leave-one-cluster-out jackknife.
# CR2/CR3 use Bell-McCaffrey Satterthwaite df under the same working model.
.vcmoe_cluster_adjust <- function(score, hessian, omega, bread, method) {
  G <- nrow(score)
  d <- ncol(score)
  method <- match.arg(method, c("CR0", "CR1", "CR2", "CR3"))
  if (G < 2L) stop("Small-cluster adjustment needs at least two active clusters.", call. = FALSE)
  inv <- solve(bread)
  identity <- diag(d)
  P <- lapply(seq_len(G), function(g) hessian[g, , ] %*% inv)
  omega_total <- Reduce(`+`, omega)
  min_value <- NA_real_
  dropped <- 0L
  M <- switch(method,
    CR0 = rep(list(identity), G),
    CR1 = rep(list(sqrt(G / (G - 1)) * identity), G),
    CR2 = lapply(seq_len(G), function(g) {
      root <- .vcmoe_psd_power(omega[[g]], 0.5)
      inverse_root <- .vcmoe_psd_power(omega[[g]], -0.5)
      expected <- omega[[g]] - P[[g]] %*% omega[[g]] - omega[[g]] %*% t(P[[g]]) +
        P[[g]] %*% omega_total %*% t(P[[g]])
      standardized <- .vcmoe_psd_power(inverse_root %*% expected %*% inverse_root, -0.5)
      min_value <<- min(min_value, attr(standardized, "min_value"), na.rm = TRUE)
      # Only count leverage-one directions, not the rank Omega_g lacks by design
      # (e.g. local slopes of cluster-constant u duplicate their intercepts).
      dropped <<- dropped + attr(standardized, "dropped") - attr(inverse_root, "dropped")
      root %*% standardized %*% inverse_root
    }),
    CR3 = lapply(seq_len(G), function(g) solve(identity - P[[g]])))
  R <- lapply(seq_len(G), function(g) inv %*% M[[g]])
  K <- lapply(seq_len(G), function(g) R[[g]] %*% P[[g]])
  adjusted <- t(vapply(seq_len(G), function(g) as.numeric(M[[g]] %*% score[g, ]), numeric(d)))
  contributions <- t(vapply(seq_len(G), function(g) as.numeric(R[[g]] %*% score[g, ]), numeric(d)))
  df <- switch(method, CR0 = rep(Inf, d), CR1 = rep(G - 1, d), {
    # V_hat_k = sum_g (sum_h b_gh U_h)^2 with b_gh = R_g[k, ] (delta_gh I - P_g);
    # Q = B Omega B' gives df = tr(Q)^2 / tr(Q^2) under Gaussian U.
    vapply(seq_len(d), function(k) {
      Rk <- t(vapply(R, function(x) x[k, ], numeric(d)))
      Kk <- t(vapply(K, function(x) x[k, ], numeric(d)))
      V <- t(vapply(seq_len(G), function(g) as.numeric(omega[[g]] %*% Rk[g, ]), numeric(d)))
      cross <- V %*% t(Kk)
      Q <- diag(rowSums(V * Rk), G) - cross - t(cross) + Kk %*% omega_total %*% t(Kk)
      sum(diag(Q))^2 / sum(Q * Q)
    }, numeric(1L))
  })
  root <- .vcmoe_psd_power(bread, -0.5)
  leverage <- max(vapply(seq_len(G), function(g) {
    max(eigen(root %*% hessian[g, , ] %*% root, symmetric = TRUE, only.values = TRUE)$values)
  }, numeric(1L)))
  colnames(adjusted) <- colnames(contributions) <- colnames(score)
  list(method = method, score = adjusted, contributions = contributions, df = df,
       inv = inv, R = R, K = K, max_cluster_leverage = leverage,
       cr2_min_eigenvalue = min_value, cr2_dropped_directions = dropped)
}

# Webb six-point weights, one column per cluster. When the 6^G distinct weight
# vectors number no more than B, all of them are used (exact enumeration)
# instead of B random draws.
.vcmoe_webb_values <- c(-sqrt(1.5), -1, -sqrt(.5), sqrt(.5), 1, sqrt(1.5))

.vcmoe_webb_weights <- function(B, G) {
  if (6^G <= B) {
    grid <- as.matrix(expand.grid(rep(list(.vcmoe_webb_values), G), KEEP.OUT.ATTRS = FALSE))
    dimnames(grid) <- NULL
    return(structure(grid, enumerated = TRUE))
  }
  structure(matrix(sample(.vcmoe_webb_values, B * G, replace = TRUE), nrow = B), enumerated = FALSE)
}

# Studentized wild cluster score bootstrap (Webb weights) for one grid:
# numerator A^{-1} sum_g w_g U_tilde_g; denominator re-applies the same
# adjustment to bootstrap residual scores w_g U_tilde_g - P_g sum_h w_h U_tilde_h.
.vcmoe_wild_tstat <- function(adjustment, weights) {
  S <- weights %*% adjustment$score
  numerator <- S %*% t(adjustment$inv)
  variance <- matrix(0, nrow(S), ncol(S))
  for (g in seq_len(nrow(adjustment$score))) {
    rows <- outer(weights[, g], as.numeric(adjustment$R[[g]] %*% adjustment$score[g, ])) -
      S %*% t(adjustment$K[[g]])
    variance <- variance + rows^2
  }
  out <- numerator / sqrt(variance)
  colnames(out) <- colnames(adjustment$score)
  out
}

# Studentized wild draws of the constancy contrast for each path set. Every
# draw contrasts each grid's numerator with its path mean, and studentizes it
# by the same contrast of the per-grid bootstrap residual rows (the
# `.vcmoe_wild_tstat` denominator), cluster by cluster. Clusters inactive at a
# grid contribute zero rows there. Returns, per set, a B x keys matrix of
# maxima over the set's grids.
.vcmoe_wild_constancy_draws <- function(adjustments, active, weights, sets) {
  grids <- sort(unique(unlist(lapply(sets, `[[`, "grids"))))
  keys <- unique(unlist(lapply(sets, `[[`, "keys")))
  B <- nrow(weights)
  col <- lapply(grids, function(j) match(keys, colnames(adjustments[[j]]$score)))
  S <- lapply(grids, function(j) weights[, active[[j]], drop = FALSE] %*% adjustments[[j]]$score)
  numerator <- lapply(seq_along(grids), function(i) {
    S[[i]] %*% t(adjustments[[grids[i]]]$inv[col[[i]], , drop = FALSE])
  })
  position <- lapply(sets, function(set) match(set$grids, grids))
  variance <- lapply(position, function(pos) rep(list(matrix(0, B, length(keys))), length(pos)))
  zero <- matrix(0, B, length(keys))
  for (g in seq_len(ncol(weights))) {
    rows <- lapply(seq_along(grids), function(i) {
      a <- adjustments[[grids[i]]]
      idx <- match(g, active[[grids[i]]])
      if (is.na(idx)) return(zero)
      outer(weights[, g], as.numeric(a$R[[idx]] %*% a$score[idx, ])[col[[i]]]) -
        S[[i]] %*% t(a$K[[idx]][col[[i]], , drop = FALSE])
    })
    for (s in seq_along(sets)) {
      pos <- position[[s]]
      center <- Reduce(`+`, rows[pos]) / length(pos)
      for (p in seq_along(pos)) variance[[s]][[p]] <- variance[[s]][[p]] + (rows[[pos[p]]] - center)^2
    }
  }
  lapply(seq_along(sets), function(s) {
    pos <- position[[s]]
    center <- Reduce(`+`, numerator[pos]) / length(pos)
    t <- lapply(seq_along(pos), function(p) abs(numerator[[pos[p]]] - center) / sqrt(variance[[s]][[p]]))
    stats::setNames(as.data.frame(do.call(pmax, t)), keys)[, sets[[s]]$keys, drop = FALSE]
  })
}

.vcmoe_wild_path_critical <- function(adjustments, active, path_keys, path_sets, level, B, seed,
                                      n_clusters, observed = NULL) {
  .vcmoe_check_multiplier_B(B)
  # One weight per cluster, shared by every grid point and coefficient.
  weights <- .vcmoe_seeded(seed, .vcmoe_webb_weights(B, n_clusters))
  stats <- lapply(seq_len(max(unlist(path_sets))), function(j) {
    if (j > length(adjustments) || is.null(adjustments[[j]])) return(NULL)
    .vcmoe_wild_tstat(adjustments[[j]], weights[, active[[j]], drop = FALSE])
  })
  maxima <- list()
  out <- lapply(path_sets, function(grids) {
    keys <- path_keys[[paste(grids, collapse = ",")]] %||% character()
    stats::setNames(vapply(keys, function(key) {
      maxima[[paste(c(grids, key), collapse = ",")]] <<- do.call(pmax, lapply(grids, function(j) abs(stats[[j]][, key])))
      as.numeric(stats::quantile(maxima[[paste(c(grids, key), collapse = ",")]], level, names = FALSE))
    }, numeric(1L)), keys)
  })
  if (is.null(observed)) return(out)
  # Zero tests reuse the band maxima; constancy needs its own studentization.
  sets <- Filter(function(set) length(set$keys) > 0L, lapply(path_sets, function(grids) {
    list(grids = grids, keys = path_keys[[paste(grids, collapse = ",")]] %||% character())
  }))
  constant_draws <- if (length(sets)) .vcmoe_wild_constancy_draws(adjustments, active, weights, sets) else list()
  names(constant_draws) <- names(sets)
  for (set in names(sets)) {
    grids <- sets[[set]]$grids
    keys <- sets[[set]]$keys
    stat <- lapply(keys, function(key) {
      influence <- do.call(cbind, lapply(grids, function(j) observed$influence[[as.character(j)]][, key]))
      .vcmoe_path_observed(influence, vapply(grids, function(j) observed$estimate[[as.character(j)]][[key]], numeric(1L)))
    })
    zero_draws <- do.call(cbind, lapply(keys, function(key) maxima[[paste(c(grids, key), collapse = ",")]]))
    attr(out[[set]], "tests") <- .vcmoe_path_test_frame(keys, vapply(stat, `[[`, numeric(1L), "zero"), zero_draws,
                                                          vapply(stat, `[[`, numeric(1L), "constant"),
                                                          as.matrix(constant_draws[[set]]))
  }
  out
}

.vcmoe_cluster_method <- function(control, is_cluster) {
  method <- toupper(control$covariance_adjustment %||% if (is_cluster) "CR2" else "HC0")
  if (identical(method, "NONE")) method <- "HC0"
  allowed <- if (is_cluster) c("CR0", "CR1", "CR2", "CR3") else "HC0"
  path <- .vcmoe_dev_simultaneous_method(control$simultaneous_method %||%
    if (is_cluster) "wild_cluster_path" else "multiplier_coefficient_path")
  paths <- c("multiplier_coefficient_path", if (is_cluster) "wild_cluster_path")
  if (!method %in% allowed || !path %in% paths) {
    stop(if (is_cluster) {
      "`cluster` requires CR0, CR1, CR2 or CR3 with `multiplier_coefficient_path` or `wild_cluster_path`; legacy critical values are not cluster calibrated."
    } else "Cell-level path inference is HC0-only with `multiplier_coefficient_path`.", call. = FALSE)
  }
  list(method = method, path = path)
}

.vcmoe_cluster_path_sets <- function(sets, n_grid) {
  sets <- sets %||% list()
  if (!is.list(sets) || (length(sets) && (is.null(names(sets)) || any(!nzchar(names(sets))) ||
                                            anyDuplicated(names(sets)) || "full" %in% names(sets)))) {
    stop("`path_sets` must be a uniquely named list of grid IDs (not named \"full\").", call. = FALSE)
  }
  for (grids in sets) {
    if (!is.numeric(grids) || !length(grids) || anyNA(grids) || any(grids != as.integer(grids)) ||
        any(grids < 1L | grids > n_grid) || anyDuplicated(grids)) {
      stop("Each path set must contain distinct grid IDs.", call. = FALSE)
    }
  }
  c(list(full = seq_len(n_grid)), lapply(sets, function(x) sort(as.integer(x))))
}

.vcmoe_cluster_intervals <- function(fit, cluster = NULL, level = .95, strict = TRUE,
                                     control = list(), moments = NULL) {
  .vcmoe_cluster_require_fit(fit)
  group <- .vcmoe_cluster_index(fit, cluster)
  is_cluster <- !is.null(cluster)
  choice <- .vcmoe_cluster_method(control, is_cluster)
  method <- choice$method
  wild <- identical(choice$path, "wild_cluster_path")
  # CR0 with Gaussian multipliers keeps the original bread and score algebra.
  adjust <- method != "HC0" && (method != "CR0" || wild)
  path_sets <- .vcmoe_cluster_path_sets(control$path_sets, length(fit$u_grid))
  ctl <- .vcmoe_dev_control(utils::modifyList(control[setdiff(names(control), "path_sets")], list(
    covariance_adjustment = "HC0", simultaneous_method = choice$path)))
  metadata_control <- ctl
  metadata_control$covariance_adjustment <- method
  .vcmoe_check_multiplier_B(ctl$multiplier_B)
  interval_rows <- diagnostic_rows <- influence <- support_rows <- adjustments <- active_sets <- list()
  z <- stats::qnorm(1 - (1 - level) / 2)
  members <- split(seq_along(group), group)
  for (j in seq_along(fit$u_grid)) {
    spec <- .vcmoe_dev_parameter_spec(fit, j)
    par <- .vcmoe_dev_parameter_vector(fit, j)
    design <- .vcmoe_cluster_design(spec$weights, group)
    few <- is_cluster && (design$active_clusters < 20 ||
                            !is.finite(design$effective_clusters) || design$effective_clusters < 10)
    pre_reason <- .vcmoe_dev_grid_block_reason(fit, j, ctl)
    reason <- if (strict) pre_reason else ""
    if (design$active_clusters < 2L) reason <- paste(c(reason, "fewer_than_two_active_clusters"), collapse = ";")
    reason <- sub("^;", "", reason)
    sh <- list()
    meat <- covariance <- adjustment <- NULL
    component_support <- NULL
    if (!nzchar(reason)) {
      sh <- if (!is.null(moments)) moments[[j]] else if (adjust) {
        .vcmoe_dev_score_hessian(fit, grid_id = j, control = ctl, group = group)
      } else .vcmoe_dev_score_hessian(fit, grid_id = j, control = ctl)
      if (is.null(sh)) stop("Missing cached local moments.", call. = FALSE)
      reason <- if (identical(sh$status, "ok")) "" else sh$block_reason
      if (!nzchar(reason) && adjust && is.null(sh$cluster_hessian)) {
        stop("Small-cluster adjustment needs moments with cluster Hessian blocks.", call. = FALSE)
      }
    }
    # Keep CR0 diagnostic meat honest even when the bread is blocked.
    if (!is.null(sh$score)) {
      summed <- if (is_cluster) rowsum(sh$score, group, reorder = FALSE) else sh$score
      meat <- crossprod(summed)
    }
    if (!nzchar(reason) && adjust) {
      active <- which(rowsum(spec$weights, group, reorder = FALSE)[, 1L] > 0)
      # Working model: independent observations within a cluster.
      omega <- lapply(active, function(g) crossprod(sh$score[members[[g]], , drop = FALSE]))
      adjustment <- tryCatch(.vcmoe_cluster_adjust(summed[active, , drop = FALSE],
        sh$cluster_hessian[active, , , drop = FALSE], omega, sh$bread, method), error = function(e) e)
      if (inherits(adjustment, "error")) {
        adjustment <- NULL
        reason <- "small_cluster_adjustment_failed"
      }
    }
    rows <- .vcmoe_dev_blocked_interval_rows(par, spec, level, reason, sh,
                                            metadata_control)
    rows$pointwise_status <- "blocked"
    rows$pointwise_block_reason <- reason
    rows$simultaneous_status <- "blocked"
    rows$simultaneous_block_reason <- "incomplete_grid_path"
    rows$path_complete <- FALSE
    rows$fit_warning <- if (nzchar(pre_reason)) pre_reason else NA_character_
    rows$influence_effective_clusters <- rows$max_cluster_influence_share <- NA_real_
    rows$df <- rows$pointwise_critical <- NA_real_
    if (!nzchar(reason)) {
      if (is.null(adjustment)) {
        inv <- solve(sh$bread)
        contributions <- summed %*% t(inv)
        df <- rep(Inf, length(par))
      } else {
        contributions <- matrix(0, nrow(summed), length(par), dimnames = list(NULL, names(par)))
        contributions[active, ] <- adjustment$contributions
        df <- adjustment$df
        adjustments[[j]] <- adjustment
        active_sets[[j]] <- active
      }
      covariance <- crossprod(contributions)
      variance <- diag(covariance)
      critical <- rep(z, length(df))
      finite_df <- is.finite(df) & df > 0
      critical[finite_df] <- stats::qt(1 - (1 - level) / 2, df[finite_df])
      critical[is.na(df) | df <= 0] <- NA_real_
      valid <- is.finite(variance) & variance > 0 & is.finite(critical)
      se <- sqrt(pmax(variance, 0))
      rows$se[valid] <- se[valid]
      rows$df[valid] <- df[valid]
      rows$pointwise_critical[valid] <- critical[valid]
      rows$pointwise_lower[valid] <- as.numeric(par[valid]) - critical[valid] * se[valid]
      rows$pointwise_upper[valid] <- as.numeric(par[valid]) + critical[valid] * se[valid]
      rows$pointwise_status[valid] <- "ok"
      rows$pointwise_block_reason[valid] <- NA_character_
      rows$pointwise_block_reason[!valid] <- "invalid_scalar_variance"
      rows$status <- rows$pointwise_status
      rows$block_reason <- rows$pointwise_block_reason
      rows$covariance_factor <- switch(method, CR1 = design$active_clusters / (design$active_clusters - 1),
                                       CR2 = NA_real_, CR3 = NA_real_, 1)
      rows$local_effective_n <- .vcmoe_dev_effective_n(spec$weights)
      squares <- contributions^2
      rows$influence_effective_clusters[valid] <- variance[valid]^2 / colSums(squares^2)[valid]
      rows$max_cluster_influence_share[valid] <- apply(squares[, valid, drop = FALSE], 2L, max) / variance[valid]
      colnames(contributions) <- rows$parameter
      influence[[as.character(j)]] <- contributions
      params <- .vcmoe_dev_unpack_parameter(par, spec)
      log_terms <- .vcmoe_dev_family_loglik_matrix(params, spec)
      posterior <- exp(log_terms - .row_log_sum_exp(log_terms))
      mass <- rowsum(posterior * spec$weights, group, reorder = FALSE)
      component_support <- data.frame(grid_id = j, u = spec$u0,
        component = paste0("component", seq_len(fit$k)),
        effective_clusters = apply(mass, 2L, .vcmoe_dev_effective_n),
        max_cluster_weight_share = apply(mass, 2L, max) / colSums(mass))
    }
    for (name in names(design)) rows[[name]] <- design[[name]]
    rows$few_clusters <- few
    rows$cluster_warning <- if (few) "few_clusters" else NA_character_
    rows$meat_rank <- .vcmoe_cluster_rank(meat)
    rows$covariance_rank <- .vcmoe_cluster_rank(covariance)
    rows$inference_scope <- if (is_cluster) "independent_cluster_sensitivity" else "independent_cell_sensitivity"
    interval_rows[[j]] <- rows
    diagnostic_rows[[j]] <- data.frame(grid_id = j, u = spec$u0,
      status = if (all(rows$pointwise_status == "ok")) "ok" else "blocked",
      block_reason = paste(unique(stats::na.omit(rows$pointwise_block_reason)), collapse = ";"),
      fit_warning = if (nzchar(pre_reason)) pre_reason else NA_character_,
      parameter_count = length(par), hessian_condition = sh$condition %||% NA_real_,
      min_eigenvalue = sh$min_eigenvalue %||% NA_real_,
      covariance_adjustment = method, simultaneous_method = choice$path,
      meat_rank = .vcmoe_cluster_rank(meat), covariance_rank = .vcmoe_cluster_rank(covariance),
      n_clusters = design$n_clusters, active_clusters = design$active_clusters,
      effective_clusters = design$effective_clusters,
      max_cluster_weight_share = design$max_cluster_weight_share, few_clusters = few,
      min_influence_effective_clusters = if (all(is.na(rows$influence_effective_clusters))) NA_real_ else min(rows$influence_effective_clusters, na.rm = TRUE),
      max_cluster_influence_share = if (all(is.na(rows$max_cluster_influence_share))) NA_real_ else max(rows$max_cluster_influence_share, na.rm = TRUE),
      min_df = if (all(is.na(rows$df))) NA_real_ else min(rows$df, na.rm = TRUE),
      max_cluster_leverage = adjustment$max_cluster_leverage %||% NA_real_,
      cr2_min_eigenvalue = adjustment$cr2_min_eigenvalue %||% NA_real_,
      cr2_dropped_directions = adjustment$cr2_dropped_directions %||% NA_integer_)
    if (!is.null(component_support)) support_rows[[length(support_rows) + 1L]] <- component_support
  }
  intervals <- do.call(rbind, interval_rows)
  keys <- unique(intervals$parameter)
  ok_grids <- lapply(keys, function(key) {
    rows <- intervals[intervals$parameter == key, ]
    if (nrow(rows) != length(fit$u_grid) || anyDuplicated(rows$grid_id)) return(integer())
    rows$grid_id[rows$pointwise_status == "ok"]
  })
  path_keys <- lapply(path_sets, function(grids) {
    keys[vapply(ok_grids, function(ok) all(grids %in% ok), logical(1L))]
  })
  names(path_keys) <- vapply(path_sets, paste, character(1L), collapse = ",")
  testing <- isTRUE(control$path_tests)
  estimates <- if (testing) lapply(split(intervals, intervals$grid_id), function(x) {
    stats::setNames(x$estimate, x$parameter)
  })
  critical <- if (wild) {
    .vcmoe_wild_path_critical(adjustments, active_sets, path_keys, path_sets, level,
                              ctl$multiplier_B, ctl$seed, max(group),
                              observed = if (testing) list(influence = influence, estimate = estimates))
  } else lapply(path_sets, function(grids) {
    selected <- path_keys[[paste(grids, collapse = ",")]]
    if (!length(selected)) return(numeric())
    combined <- do.call(cbind, lapply(grids, function(j) influence[[as.character(j)]][, selected, drop = FALSE]))
    .vcmoe_path_multiplier(combined, keys = rep(selected, times = length(grids)), level,
                           ctl$multiplier_B, ctl$seed,
                           estimate = if (testing) unlist(lapply(grids, function(j) estimates[[as.character(j)]][selected])))
  })
  tests <- if (testing) {
    labels <- unique(intervals[c("parameter", "coefficient_set", "component", "term", "block")])
    rows <- lapply(names(path_sets), function(set) {
      x <- attr(critical[[set]], "tests")
      if (is.null(x) || !nrow(x)) return(NULL)
      x <- merge(labels, x, by = "parameter", sort = FALSE)
      data.frame(x, path_set = set, n_grids = length(path_sets[[set]]), covariance_adjustment = method,
                 simultaneous_method = choice$path, multiplier_B = ctl$multiplier_B, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }
  for (set in names(path_sets)) {
    suffix <- if (identical(set, "full")) "" else paste0("_", set)
    columns <- paste0("simultaneous_", c("critical", "lower", "upper", "status"), suffix)
    if (nzchar(suffix)) {
      intervals[columns[1:3]] <- NA_real_
      intervals[[columns[4L]]] <- "blocked"
    }
    keep <- intervals$parameter %in% names(critical[[set]]) & intervals$grid_id %in% path_sets[[set]]
    intervals[[columns[1L]]][keep] <- critical[[set]][intervals$parameter[keep]]
    intervals[[columns[2L]]][keep] <- intervals$estimate[keep] - intervals[[columns[1L]]][keep] * intervals$se[keep]
    intervals[[columns[3L]]][keep] <- intervals$estimate[keep] + intervals[[columns[1L]]][keep] * intervals$se[keep]
    intervals[[columns[4L]]][keep] <- "ok"
    if (!nzchar(suffix)) {
      intervals$simultaneous_block_reason[keep] <- NA_character_
      intervals$path_complete[keep] <- TRUE
    }
  }
  diagnostics <- do.call(rbind, diagnostic_rows)
  diagnostics$complete_paths <- length(path_keys[[1L]])
  diagnostics$total_paths <- length(keys)
  if (is_cluster && any(diagnostics$few_clusters)) warning(
    "Few locally supported independent clusters: intervals are sensitivity analyses, not small-cluster calibrated inference.", call. = FALSE)
  out <- list(intervals = intervals, diagnostics = diagnostics,
              cluster_support = if (length(support_rows)) do.call(rbind, support_rows) else data.frame(),
              tests = tests)
  out
}

# Paired comparison of covariance choices from one numerical score/Hessian
# pass (HC0 analytic, HC0 multiplier, CR0 multiplier). `small_sample` adds
# CR0-CR3 intervals with studentized wild cluster paths; their moments carry
# cluster Hessian blocks, from the same function values.
.vcmoe_confband_comparison <- function(fit, cluster, level = .95,
                                        coefficient_set = c("expert", "gating", "sigma", "theta"),
                                        strict = TRUE, control = list(), small_sample = character()) {
  .vcmoe_cluster_require_fit(fit)
  group <- .vcmoe_cluster_index(fit, cluster)
  if (is.null(cluster)) stop("Comparison requires cluster IDs.", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) stop("Invalid confidence level.", call. = FALSE)
  small_sample <- toupper(small_sample)
  if (length(small_sample) && (!all(small_sample %in% c("CR0", "CR1", "CR2", "CR3")) || anyDuplicated(small_sample))) {
    stop("`small_sample` must name distinct methods among CR0, CR1, CR2 and CR3.", call. = FALSE)
  }
  ctl <- .vcmoe_dev_control(utils::modifyList(control[setdiff(names(control), "path_sets")],
                                              list(covariance_adjustment = "HC0")))
  moments <- .vcmoe_inference_moments(fit, strict, ctl, group = if (length(small_sample)) group)
  dl <- fit$u_scaling$domain_length %||% vcmoe_parameterization(fit)$u_scaling$domain_length %||% 1
  analytic_control <- utils::modifyList(ctl, list(simultaneous_method = "analytic_epanechnikov_path"))
  analytic <- .vcmoe_dev_intervals(fit, level = level, strict = strict,
                                  control = analytic_control, moments = moments)
  analytic$intervals$pointwise_status <- analytic$intervals$status
  analytic$intervals$pointwise_block_reason <- analytic$intervals$block_reason
  analytic$intervals$simultaneous_status <- analytic$intervals$status
  analytic$intervals$simultaneous_block_reason <- analytic$intervals$block_reason
  for (key in unique(analytic$intervals$parameter)) {
    ids <- which(analytic$intervals$parameter == key)
    if (any(!analytic$intervals$status[ids] %in% c("ok", "warning"))) {
      analytic$intervals$simultaneous_status[ids] <- "blocked"
      analytic$intervals$simultaneous_block_reason[ids] <- "incomplete_grid_path"
      analytic$intervals$simultaneous_lower[ids] <- analytic$intervals$simultaneous_upper[ids] <- NA_real_
    }
  }
  path_control <- utils::modifyList(ctl, list(simultaneous_method = "multiplier_coefficient_path",
                                              path_sets = control$path_sets))
  cell <- .vcmoe_cluster_intervals(fit, level = level, strict = strict,
                                  control = path_control, moments = moments)
  clustered <- .vcmoe_cluster_intervals(fit, cluster, level, strict,
                                       utils::modifyList(path_control, list(covariance_adjustment = "CR0")), moments)
  out <- list(hc0_analytic = .vcmoe_confband_pack(fit, analytic, level, "simultaneous", coefficient_set, analytic_control, dl),
       hc0_multiplier = .vcmoe_confband_pack(fit, cell, level, "simultaneous", coefficient_set, path_control, dl),
       cluster_multiplier = .vcmoe_confband_pack(fit, clustered, level, "simultaneous", coefficient_set, path_control, dl))
  for (method in small_sample) {
    wild_control <- utils::modifyList(path_control, list(covariance_adjustment = method,
                                                         simultaneous_method = "wild_cluster_path"))
    result <- .vcmoe_cluster_intervals(fit, cluster, level, strict, wild_control, moments)
    out[[paste0(tolower(method), "_wild")]] <- .vcmoe_confband_pack(fit, result, level, "simultaneous",
                                                                    coefficient_set, wild_control, dl)
  }
  out
}
