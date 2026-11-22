# Null-anchored, cluster-robust gene score screen.
#
# The null VCMoE is fitted once. For each candidate covariate x_j, every grid
# point gets the local model augmented by x_j b, b = (1, du) the local-linear
# basis, with new coefficients at zero: in each expert's linear predictor and,
# for the expert-side test, also in each gating contrast. The tested
# coefficients are the expert ones; the gene's gating coefficients are
# nuisance, so a gene that only marks the latent components (and therefore
# predicts the response through composition alone) is absorbed by its gating
# direction instead of being flagged as a within-component association.
#
# Blocks are exact derivatives of the local likelihood the CR2 bands use:
#   d(w_i l_i)/d psi_{a,b}                 =  w_i b_i x_ij G_ia,
#   -d2(w_i l_i)/d psi_{a,b} d par         = -w_i b_i x_ij dG_ia/d par,
#   -d2(w_i l_i)/d psi_{a,b} d psi_{a',b'} = -w_i b_i b'_i x_ij^2 Q_iaa'.
# Expert family (component k): G = tau_k r_k, with r and c the NB eta score and
# curvature. Gating family (contrast column a): G = tau_a - pi_a. Second
# derivatives: expert-expert 1(k = k') tau_k (c_k + r_k^2) - tau_k tau_k' r_k r_k';
# gating-gating tau_a (1(a = a') - tau_a') - pi_a (1(a = a') - pi_a');
# expert-gating r_k tau_k (1(k = a) - tau_a). G and Q depend on the null fit
# only, so all genes' blocks come from one matrix product per cluster.
#
# Under an alternative the observed information of a mixture in the tested
# directions can be indefinite at the null fit (a saddle), and a Newton step
# from there is meaningless. information = "complete" (default) sets the
# efficient information of the tested coefficients given everything else
# (null parameters and the gene's nuisance directions) to their complete-data
# information J^c + ridge (expert: 1(k = k') tau_k c_k; gating:
# -pi_a (1(a = a') - pi_a')), by taking the tested block as
# J^c + A_tn A_nn^{-1} A_nt, split symmetrically over clusters. Scores and
# every other block stay exact, so the nuisance projection of the efficient
# score is unchanged; the step is a Fisher-scoring/EM step in the tested
# coefficients. information = "observed" keeps the exact block.
#
# The per-gene statistic is the one-step estimate Delta = A^{-1} (sum_g U_g +
# penalty gradient) of the tested level coefficients, studentized per grid
# under the CR2 small-cluster adjustment, with the studentized wild cluster
# bootstrap (Webb weights shared by every gene, grid and component) as in the
# CR2 path tests. studentize = "onestep" (default) feeds CR2 the cluster
# scores linearized to the one-step point, U_g - J_g Delta, which sum to the
# penalty gradient there as fitted scores do (the Wald form of the path tests).
# studentize = "restricted" keeps the null-fit scores (score form); those are
# uncentred under an alternative, so by Cauchy-Schwarz |t| <= sqrt(active
# clusters) at every grid point, which with few clusters caps the statistic
# below any useful critical value. The working covariance Omega_g stays at the
# null fit in both forms.

.vcmoe_screen_side <- function(side) match.arg(side, c("expert", "gating"))

.vcmoe_screen_contrast_columns <- function(k) if (k == 2L) 1L else seq.int(2L, k)

.vcmoe_screen_contrast_labels <- function(k) {
  if (k == 2L) "component1_vs_component2" else paste0("component", seq.int(2L, k), "_vs_component1")
}

# Per-cell G (n x m), observed Q and complete-data Qc (n x m x m) for the
# requested families, expert columns first.
.vcmoe_screen_cell_terms <- function(par, spec, families) {
  params <- .vcmoe_dev_unpack_parameter(par, spec)
  log_terms <- .vcmoe_dev_family_loglik_matrix(params, spec)
  tau <- exp(log_terms - .row_log_sum_exp(log_terms))
  n <- nrow(tau)
  k <- spec$k
  expert <- "expert" %in% families
  gating <- "gating" %in% families
  cols <- .vcmoe_screen_contrast_columns(k)
  G <- NULL
  if (expert) {
    design <- .make_expert_design(spec$z_design, spec$du)
    ctl <- spec$control
    theta <- pmin(pmax(params$theta, ctl$negbin_theta_min), ctl$negbin_theta_max)
    r <- curvature <- matrix(0, n, k)
    for (comp in seq_len(k)) {
      eta <- spec$offset + as.numeric(design %*% params$expert_coef[comp, ])
      inside <- as.numeric(eta > ctl$negbin_eta_min & eta < ctl$negbin_eta_max)
      mu <- .bounded_exp(eta, ctl)
      r[, comp] <- inside * (spec$y - mu) * theta[[comp]] / (theta[[comp]] + mu)
      curvature[, comp] <- -inside * mu * theta[[comp]] * (theta[[comp]] + spec$y) / (theta[[comp]] + mu)^2
    }
    G <- tau * r
  }
  if (gating) {
    prior <- .gating_prob(.make_expert_design(spec$x_design, spec$du), params$gating_coef)
    G <- cbind(G, tau[, cols, drop = FALSE] - prior[, cols, drop = FALSE])
  }
  family <- c(if (expert) rep("expert", k), if (gating) rep("gating", length(cols)))
  index <- c(if (expert) seq_len(k), if (gating) cols)
  m <- length(family)
  Q <- Qc <- array(0, c(n, m, m))
  for (a in seq_len(m)) for (b in seq_len(m)) {
    i <- index[[a]]
    j <- index[[b]]
    same <- as.numeric(i == j)
    if (family[[a]] == "expert" && family[[b]] == "expert") {
      Q[, a, b] <- same * tau[, i] * (curvature[, i] + r[, i]^2) - tau[, i] * tau[, j] * r[, i] * r[, j]
      Qc[, a, b] <- same * tau[, i] * curvature[, i]
    } else if (family[[a]] == "gating" && family[[b]] == "gating") {
      Q[, a, b] <- tau[, i] * (same - tau[, j]) - prior[, i] * (same - prior[, j])
      Qc[, a, b] <- -prior[, i] * (same - prior[, j])
    } else if (family[[a]] == "expert") {
      Q[, a, b] <- r[, i] * tau[, i] * (same - tau[, j])
    } else {
      Q[, a, b] <- r[, j] * tau[, j] * (same - tau[, i])
    }
  }
  labels <- c(if (expert) spec$component_labels, if (gating) .vcmoe_screen_contrast_labels(k))
  list(G = G, Q = Q, Qc = Qc, family = family, labels = labels)
}

# Central differences of G in every null parameter: n x m x d.
.vcmoe_screen_cell_derivative <- function(par, spec, families, eps) {
  steps <- eps * pmax(1, abs(par))
  base <- .vcmoe_screen_cell_terms(par, spec, families)$G
  out <- array(NA_real_, c(dim(base), length(par)))
  for (a in seq_along(par)) {
    e <- replace(numeric(length(par)), a, steps[[a]])
    out[, , a] <- (.vcmoe_screen_cell_terms(par + e, spec, families)$G -
                     .vcmoe_screen_cell_terms(par - e, spec, families)$G) / (2 * steps[[a]])
  }
  out
}

# All-gene moments at one grid point. New directions are ordered level
# (families 1..m) then slope (families 1..m); block arrays are indexed
# [active cluster, gene, ...] with the new direction outermost. `tested` is
# the family whose coefficients are tested; the others are nuisance.
.vcmoe_screen_grid_moments <- function(fit, grid_id, X, group, families, tested, ctl, information = "complete") {
  reason <- if (isTRUE(ctl$strict)) .vcmoe_dev_grid_block_reason(fit, grid_id, ctl) else ""
  if (nzchar(reason)) return(list(status = "blocked", reason = reason))
  sh <- .vcmoe_dev_score_hessian(fit, grid_id = grid_id, control = ctl, group = group)
  if (!identical(sh$status, "ok")) return(list(status = "blocked", reason = sh$block_reason))
  spec <- sh$spec
  par <- sh$parameter
  d <- length(par)
  totals <- rowsum(spec$weights, group, reorder = FALSE)[, 1L]
  clusters <- which(totals > 0)
  if (length(clusters) < 2L) return(list(status = "blocked", reason = "fewer_than_two_active_clusters"))
  active <- which(spec$weights > 0)
  # Per-cell terms are needed on the kernel support only.
  spec_active <- .vcmoe_dev_spec_rows(spec, active)
  terms <- .vcmoe_screen_cell_terms(par, spec_active, families)
  D <- .vcmoe_screen_cell_derivative(par, spec_active, families, ctl$finite_diff_eps)
  m <- ncol(terms$G)
  e <- 2L * m
  na <- length(active)
  w <- spec$weights[active]
  basis <- cbind(1, spec$du[active])
  dir_b <- rep(1:2, each = m)
  dir_m <- rep(seq_len(m), 2L)
  V <- matrix(vapply(seq_len(e), function(a) w * basis[, dir_b[a]] * terms$G[, dir_m[a]], numeric(na)), nrow = na)
  S <- sh$score[active, , drop = FALSE]
  Wd <- do.call(cbind, lapply(seq_len(e), function(a) {
    -(w * basis[, dir_b[a]]) * matrix(D[, dir_m[a], ], nrow = na)
  }))
  VS <- do.call(cbind, lapply(seq_len(e), function(a) V[, a] * S))
  pairs <- expand.grid(second = seq_len(e), first = seq_len(e))
  second_order <- function(Q) {
    matrix(vapply(seq_len(nrow(pairs)), function(r) {
      a <- pairs$first[r]
      b <- pairs$second[r]
      -w * basis[, dir_b[a]] * basis[, dir_b[b]] * Q[, dir_m[a], dir_m[b]]
    }, numeric(na)), nrow = na)
  }
  Wq <- second_order(terms$Q)
  Wc <- if (identical(information, "complete")) second_order(terms$Qc) else NULL
  VV <- matrix(vapply(seq_len(nrow(pairs)), function(r) V[, pairs$first[r]] * V[, pairs$second[r]], numeric(na)), nrow = na)
  Xa <- X[active, , drop = FALSE]
  Xa2 <- Xa * Xa
  ga <- group[active]
  p <- ncol(X)
  Ga <- length(clusters)
  U <- array(0, c(Ga, p, e))
  Hx <- Ox <- array(0, c(Ga, p, e * d))
  Hn <- On <- array(0, c(Ga, p, e * e))
  Hc <- if (is.null(Wc)) NULL else array(0, c(Ga, p, e * e))
  omega_null <- vector("list", Ga)
  for (i in seq_len(Ga)) {
    rows <- which(ga == clusters[i])
    x <- Xa[rows, , drop = FALSE]
    x2 <- Xa2[rows, , drop = FALSE]
    U[i, , ] <- crossprod(x, V[rows, , drop = FALSE])
    Hx[i, , ] <- crossprod(x, Wd[rows, , drop = FALSE])
    Ox[i, , ] <- crossprod(x, VS[rows, , drop = FALSE])
    Hn[i, , ] <- crossprod(x2, Wq[rows, , drop = FALSE])
    if (!is.null(Hc)) Hc[i, , ] <- crossprod(x2, Wc[rows, , drop = FALSE])
    On[i, , ] <- crossprod(x2, VV[rows, , drop = FALSE])
    omega_null[[i]] <- crossprod(S[rows, , drop = FALSE])
  }
  ridge_family <- c(expert = if (identical(spec$family, "negative-binomial")) spec$control$negbin_ridge else spec$control$ridge,
                    gating = spec$control$ridge)
  dir_family <- terms$family[dir_m]
  set <- c(expert = "expert", gating = "gating_contrast")
  list(status = "ok", reason = "", grid_id = grid_id, u0 = spec$u0, par = par, d = d, m = m, e = e,
       information = information, tested = dir_family == tested, dir_b = dir_b, dir_family = dir_family,
       dir_names = paste(set[dir_family], terms$labels[dir_m], "%s", c("intercept", "slope")[dir_b], sep = ":"),
       labels = terms$labels, family = terms$family,
       bread_null = sh$bread, clusters = clusters,
       score_null = rowsum(sh$score, group, reorder = FALSE)[clusters, , drop = FALSE],
       hessian_null = sh$cluster_hessian[clusters, , , drop = FALSE], omega_null = omega_null,
       max_condition = ctl$max_condition, min_eigenvalue = ctl$min_eigenvalue,
       penalty = -spec$ridge * (par - spec$ridge_target), ridge_new = unname(ridge_family[dir_family]),
       U = U, Hx = Hx, Ox = Ox, Hn = Hn, Hc = Hc, On = On,
       newton_step = max(abs(solve(sh$bread, colSums(sh$score) - spec$ridge * (par - spec$ridge_target)))),
       effective_clusters = .vcmoe_dev_effective_n(totals[clusters]))
}

# Augmented (null + every new direction of gene j) score, cluster Hessian,
# working covariance and bread at one grid point, in the layout
# .vcmoe_cluster_adjust expects. Null parameters come first.
.vcmoe_screen_gene_blocks <- function(mom, j) {
  d <- mom$d
  e <- mom$e
  null <- seq_len(d)
  new <- d + seq_len(e)
  Ga <- length(mom$clusters)
  score <- cbind(mom$score_null, matrix(mom$U[, j, ], nrow = Ga))
  hx <- lapply(seq_len(Ga), function(i) matrix(mom$Hx[i, j, ], e, d, byrow = TRUE))
  hn <- lapply(seq_len(Ga), function(i) {
    x <- matrix(mom$Hn[i, j, ], e, e, byrow = TRUE)
    (x + t(x)) / 2
  })
  if (identical(mom$information, "complete")) {
    t_dir <- which(mom$tested)
    n_dir <- which(!mom$tested)
    # A_nn over the null parameters and the gene's nuisance directions.
    cross_total <- Reduce(`+`, hx)
    hn_total <- Reduce(`+`, hn)
    nuisance <- rbind(cbind(mom$bread_null, t(cross_total[n_dir, , drop = FALSE])),
                      cbind(cross_total[n_dir, , drop = FALSE],
                            hn_total[n_dir, n_dir, drop = FALSE] + diag(mom$ridge_new[n_dir], length(n_dir))))
    tested_cross <- cbind(cross_total[t_dir, , drop = FALSE], hn_total[t_dir, n_dir, drop = FALSE])
    projection <- solve(nuisance, t(tested_cross))
    for (i in seq_len(Ga)) {
      jc <- matrix(mom$Hc[i, j, ], e, e, byrow = TRUE)[t_dir, t_dir, drop = FALSE]
      share <- cbind(hx[[i]][t_dir, , drop = FALSE], hn[[i]][t_dir, n_dir, drop = FALSE]) %*% projection
      hn[[i]][t_dir, t_dir] <- (jc + t(jc)) / 2 + (share + t(share)) / 2
    }
  }
  hessian <- array(0, c(Ga, d + e, d + e))
  omega <- vector("list", Ga)
  for (i in seq_len(Ga)) {
    h <- matrix(0, d + e, d + e)
    h[null, null] <- mom$hessian_null[i, , ]
    h[new, null] <- hx[[i]]
    h[null, new] <- t(hx[[i]])
    h[new, new] <- hn[[i]]
    hessian[i, , ] <- h
    ox <- matrix(mom$Ox[i, j, ], e, d, byrow = TRUE)
    on <- matrix(mom$On[i, j, ], e, e, byrow = TRUE)
    o <- matrix(0, d + e, d + e)
    o[null, null] <- mom$omega_null[[i]]
    o[new, null] <- ox
    o[null, new] <- t(ox)
    o[new, new] <- (on + t(on)) / 2
    omega[[i]] <- o
  }
  bread <- matrix(0, d + e, d + e)
  bread[null, null] <- mom$bread_null
  bread[new, null] <- Reduce(`+`, hx)
  bread[null, new] <- t(bread[new, null])
  bread[new, new] <- Reduce(`+`, hn) + diag(mom$ridge_new, e)
  list(score = score, hessian = hessian, omega = omega, bread = bread, new = new)
}

.vcmoe_screen_new_names <- function(mom, gene) sprintf(mom$dir_names, gene)

# One gene at one grid: one-step estimates, CR2 influence and BM df of the
# tested coefficients (level then slope); with `weights`, the studentized wild
# numerators and per-cluster bootstrap residual rows of the tested levels.
.vcmoe_screen_gene_grid <- function(mom, j, weights = NULL, studentize = "onestep") {
  tested <- which(mom$tested)
  level <- which(mom$tested & mom$dir_b == 1L)
  blocked <- function() {
    list(blocked = TRUE, estimate = rep(NA_real_, length(tested)), se = rep(NA_real_, length(tested)),
         df = rep(NA_real_, length(tested)))
  }
  blocks <- .vcmoe_screen_gene_blocks(mom, j)
  condition <- .vcmoe_dev_matrix_condition(blocks$bread)
  if (!is.finite(condition$condition) || condition$condition > mom$max_condition ||
      !is.finite(condition$min_eigenvalue) || condition$min_eigenvalue <= mom$min_eigenvalue) return(blocked())
  total <- colSums(blocks$score) + c(mom$penalty, numeric(mom$e))
  score <- blocks$score
  if (identical(studentize, "onestep")) {
    delta <- solve(blocks$bread, total)
    score <- score - t(vapply(seq_len(nrow(score)), function(i) as.numeric(blocks$hessian[i, , ] %*% delta),
                              numeric(ncol(score))))
  }
  adj <- tryCatch(.vcmoe_cluster_adjust(score, blocks$hessian, blocks$omega, blocks$bread, "CR2"),
                  error = function(e) NULL)
  if (is.null(adj)) return(blocked())
  cols <- mom$d + tested
  influence <- adj$contributions[, cols, drop = FALSE]
  out <- list(blocked = FALSE, estimate = as.numeric(adj$inv[cols, , drop = FALSE] %*% total),
              influence = influence, se = sqrt(colSums(influence^2)), df = adj$df[cols])
  if (is.null(weights)) return(out)
  lc <- mom$d + level
  Wg <- weights[, mom$clusters, drop = FALSE]
  Sw <- Wg %*% adj$score
  out$numerator <- Sw %*% t(adj$inv[lc, , drop = FALSE])
  out$rows <- lapply(seq_along(mom$clusters), function(i) {
    outer(Wg[, i], as.numeric(adj$R[[i]] %*% adj$score[i, ])[lc]) - Sw %*% t(adj$K[[i]][lc, , drop = FALSE])
  })
  out
}

#' Null-anchored gene score screen (internal)
#'
#' @param fit A single-response NB local-grid VCMoE fit with `keep_data = TRUE`
#'   (the null model).
#' @param X Numeric matrix of candidate covariates, one row per fitted row and
#'   one named column per gene. Columns are centred internally.
#' @param cluster One cluster id per fitted row.
#' @param grids Grid ids that enter the path statistics.
#' @param side "expert" (tested) or "gating" (diagnostic).
#' @param B Webb wild draws; 0 skips the bootstrap.
#' @param seed Seed of the shared Webb weights (B x clusters).
#' @param studentize "onestep" (Wald form, default) or "restricted" (score form).
#' @param information "complete" (default) or "observed" tested block.
#' @param adjust_gating Expert side: also augment the gating by the gene and
#'   treat those coefficients as nuisance (default).
#' @noRd
.vcmoe_score_screen <- function(fit, X, cluster, grids = seq_along(fit$u_grid), side = "expert",
                                B = 9999L, seed = NULL, control = list(), center = TRUE,
                                studentize = c("onestep", "restricted"), information = c("complete", "observed"),
                                adjust_gating = TRUE, weights = NULL) {
  .vcmoe_cluster_require_fit(fit)
  side <- .vcmoe_screen_side(side)
  studentize <- match.arg(studentize)
  information <- match.arg(information)
  if (identical(side, "expert") && !identical(fit$family, "negative-binomial")) {
    stop("The expert-side screen is implemented for negative-binomial experts.", call. = FALSE)
  }
  X <- as.matrix(X)
  if (!is.numeric(X) || nrow(X) != length(fit$fitted$u) || any(!is.finite(X)) ||
      is.null(colnames(X)) || anyDuplicated(colnames(X))) {
    stop("`X` must be a finite numeric matrix with one row per fitted row and unique column names.", call. = FALSE)
  }
  if (identical(side, "expert") && isTRUE(adjust_gating)) {
    # A gene already in the null gating keeps its fitted gating effect there;
    # augmenting the gating again would duplicate that column.
    in_gating <- .vcmoe_screen_in_gating(fit, X)
    if (all(in_gating)) {
      adjust_gating <- FALSE
    } else if (any(in_gating)) {
      args <- list(fit = fit, cluster = cluster, grids = grids, side = side, B = B, seed = seed, control = control,
                   center = center, studentize = studentize, information = information,
                   weights = if (as.integer(B) > 0L) weights %||% .vcmoe_seeded(seed, .vcmoe_webb_weights(B, max(.vcmoe_cluster_index(fit, cluster)))))
      parts <- list(do.call(.vcmoe_score_screen, c(args, list(X = X[, !in_gating, drop = FALSE], adjust_gating = TRUE))),
                    do.call(.vcmoe_score_screen, c(args, list(X = X[, in_gating, drop = FALSE], adjust_gating = FALSE))))
      genes <- rbind(parts[[1L]]$genes, parts[[2L]]$genes)
      out <- parts[[1L]]
      out$genes <- genes[match(colnames(X), genes$gene), , drop = FALSE]
      rownames(out$genes) <- NULL
      out$pointwise <- rbind(parts[[1L]]$pointwise, parts[[2L]]$pointwise)
      if (!is.null(out$draws)) out$draws <- cbind(parts[[1L]]$draws, parts[[2L]]$draws)[, colnames(X), drop = FALSE]
      out$adjust_gating <- stats::setNames(!in_gating, colnames(X))
      return(out)
    }
  }
  families <- if (identical(side, "gating")) "gating" else if (isTRUE(adjust_gating)) c("expert", "gating") else "expert"
  if (center) X <- sweep(X, 2L, colMeans(X), "-")
  group <- .vcmoe_cluster_index(fit, cluster)
  grids <- sort(unique(as.integer(grids)))
  if (!length(grids) || any(grids < 1L | grids > length(fit$u_grid))) stop("Invalid `grids`.", call. = FALSE)
  B <- as.integer(B)
  if (B > 0L) .vcmoe_check_multiplier_B(B)
  ctl <- .vcmoe_dev_control(utils::modifyList(control, list(covariance_adjustment = "HC0")))
  if (B > 0L && is.null(weights)) weights <- .vcmoe_seeded(seed, .vcmoe_webb_weights(B, max(group)))
  if (B > 0L) {
    if (!is.matrix(weights) || ncol(weights) != max(group)) stop("`weights` must have one column per cluster.", call. = FALSE)
    B <- nrow(weights)
  } else {
    weights <- NULL
  }
  moments <- lapply(grids, function(j) .vcmoe_screen_grid_moments(fit, j, X, group, families, side, ctl, information))
  ok <- vapply(moments, function(x) identical(x$status, "ok"), logical(1L))
  grid_rows <- data.frame(grid_id = grids, u0 = fit$u_grid[grids],
                          status = vapply(moments, `[[`, character(1L), "status"),
                          block_reason = vapply(moments, function(x) x$reason %||% "", character(1L)),
                          active_clusters = vapply(moments, function(x) length(x$clusters %||% integer()), integer(1L)),
                          effective_clusters = vapply(moments, function(x) x$effective_clusters %||% NA_real_, numeric(1L)),
                          null_newton_step = vapply(moments, function(x) x$newton_step %||% NA_real_, numeric(1L)),
                          stringsAsFactors = FALSE)
  if (!any(ok)) stop("No usable grid point for the screen.", call. = FALSE)
  moments <- moments[ok]
  first <- moments[[1L]]
  components <- first$labels[first$family == side]
  m <- length(components)
  n_grid <- length(moments)
  G <- max(group)
  genes <- colnames(X)
  gene_rows <- pointwise <- vector("list", length(genes))
  draws <- if (B > 0L) matrix(NA_real_, B, length(genes), dimnames = list(NULL, genes)) else NULL
  for (j in seq_along(genes)) {
    per_grid <- lapply(moments, function(mom) .vcmoe_screen_gene_grid(mom, j, weights, studentize))
    blocked <- vapply(per_grid, `[[`, logical(1L), "blocked")
    est <- matrix(t(vapply(per_grid, function(x) x$estimate, numeric(2L * m))), nrow = n_grid)
    se <- matrix(t(vapply(per_grid, function(x) x$se, numeric(2L * m))), nrow = n_grid)
    df <- matrix(t(vapply(per_grid, function(x) x$df, numeric(2L * m))), nrow = n_grid)
    names_tested <- .vcmoe_screen_new_names(first, genes[[j]])[first$tested]
    pointwise[[j]] <- data.frame(gene = genes[[j]], grid_id = rep(grid_rows$grid_id[ok], 2L * m),
      u0 = rep(grid_rows$u0[ok], 2L * m), parameter = rep(names_tested, each = n_grid),
      component = rep(rep(components, 2L), each = n_grid),
      block = rep(rep(c("intercept", "slope"), each = m), each = n_grid),
      estimate = as.numeric(est), se = as.numeric(se), t = as.numeric(est / se), df = as.numeric(df),
      blocked = rep(blocked, 2L * m), stringsAsFactors = FALSE)
    use <- which(!blocked)
    row <- data.frame(gene = genes[[j]], side = side, adjust_gating = side == "expert" && isTRUE(adjust_gating),
                      n_grids = length(use), n_blocked = sum(blocked), statistic = NA_real_,
                      argmax_component = NA_character_, argmax_grid = NA_integer_, min_df = NA_real_,
                      stringsAsFactors = FALSE)
    for (comp in components) row[[paste0("statistic_", comp)]] <- row[[paste0("sign_", comp)]] <- NA_real_
    if (B > 0L) {
      row$p_value <- NA_real_
      for (comp in components) {
        row[[paste0("p_value_", comp)]] <- NA_real_
        row[[paste0("constancy_statistic_", comp)]] <- row[[paste0("constancy_p_value_", comp)]] <- NA_real_
      }
    }
    if (length(use)) {
      grid_use <- grid_rows$grid_id[ok][use]
      level_t <- abs(est[use, seq_len(m), drop = FALSE] / se[use, seq_len(m), drop = FALSE])
      row$statistic <- max(level_t)
      row$argmax_component <- components[[which.max(apply(level_t, 2L, max))]]
      row$argmax_grid <- grid_use[[which.max(apply(level_t, 1L, max))]]
      row$min_df <- min(df[use, seq_len(m)])
      for (a in seq_len(m)) {
        row[[paste0("statistic_", components[[a]])]] <- max(level_t[, a])
        row[[paste0("sign_", components[[a]])]] <- sign(est[use[which.max(level_t[, a])], a])
      }
    }
    if (B > 0L && length(use)) {
      pg <- per_grid[use]
      mu <- moments[use]
      t_star <- lapply(seq_len(m), function(a) {
        do.call(pmax, lapply(pg, function(x) {
          abs(x$numerator[, a]) / sqrt(Reduce(`+`, lapply(x$rows, function(r) r[, a]^2)))
        }))
      })
      gene_draws <- do.call(pmax, t_star)
      draws[, j] <- gene_draws
      row$p_value <- .vcmoe_path_p_value(gene_draws, row$statistic)
      for (a in seq_len(m)) {
        comp <- components[[a]]
        row[[paste0("p_value_", comp)]] <- .vcmoe_path_p_value(t_star[[a]], row[[paste0("statistic_", comp)]])
        # Constancy of the one-step level curve (descriptive): per-cluster
        # influence and bootstrap rows, zero where a cluster is inactive.
        if (length(use) < 2L) next
        influence <- matrix(0, G, length(use))
        rows_full <- vector("list", length(use))
        for (p in seq_along(use)) {
          influence[mu[[p]]$clusters, p] <- pg[[p]]$influence[, a]
          rows_full[[p]] <- matrix(0, B, G)
          for (i in seq_along(mu[[p]]$clusters)) rows_full[[p]][, mu[[p]]$clusters[i]] <- pg[[p]]$rows[[i]][, a]
        }
        observed <- .vcmoe_path_observed(influence, est[use, a])
        if (!is.finite(observed$constant)) next
        num <- vapply(pg, function(x) x$numerator[, a], numeric(B))
        center_num <- rowMeans(num)
        center_rows <- Reduce(`+`, rows_full) / length(rows_full)
        stat <- vapply(seq_along(use), function(p) {
          abs(num[, p] - center_num) / sqrt(rowSums((rows_full[[p]] - center_rows)^2))
        }, numeric(B))
        row[[paste0("constancy_statistic_", comp)]] <- observed$constant
        row[[paste0("constancy_p_value_", comp)]] <- .vcmoe_path_p_value(apply(stat, 1L, max), observed$constant)
      }
    }
    gene_rows[[j]] <- row
  }
  list(genes = do.call(rbind, gene_rows), pointwise = do.call(rbind, pointwise), grids = grid_rows,
       draws = draws, side = side, B = B, seed = seed, components = components, studentize = studentize,
       information = information, adjust_gating = side == "expert" && isTRUE(adjust_gating))
}

# Columns of X that reproduce a non-intercept column of the null gating design.
.vcmoe_screen_in_gating <- function(fit, X) {
  Z <- fit$fitted$x_design
  Z <- Z[, colnames(Z) != "(Intercept)" & apply(Z, 2L, stats::sd) > 0, drop = FALSE]
  if (!ncol(Z)) return(rep(FALSE, ncol(X)))
  vapply(seq_len(ncol(X)), function(j) {
    if (stats::sd(X[, j]) == 0) return(FALSE)
    any(abs(suppressWarnings(stats::cor(X[, j], Z))) > 1 - 1e-8)
  }, logical(1L))
}

# Single-step Westfall-Young minP over genes from the shared-weight sup-t
# draws (one column per gene). Each draw is turned into its within-gene
# bootstrap p-value first: with few clusters the studentized draws are
# heavy-tailed to a gene-specific degree, so a raw maxT would be governed by
# the heaviest-tailed genes rather than by the evidence of each gene.
.vcmoe_screen_minp <- function(statistic, draws) {
  if (length(statistic) != ncol(draws) || anyNA(draws)) stop("One complete draw column per statistic is required.", call. = FALSE)
  B <- nrow(draws)
  raw <- (1 + colSums(draws >= matrix(statistic, B, length(statistic), byrow = TRUE))) / (B + 1)
  # p*_bj = #{b': T*_b'j >= T*_bj} / B.
  draw_p <- apply(draws, 2L, function(x) (B - rank(x, ties.method = "min") + 1) / B)
  minimum <- apply(matrix(draw_p, nrow = B), 1L, min)
  pmax(raw, vapply(raw, function(p) (1 + sum(minimum <= p)) / (B + 1), numeric(1L)))
}

# Posterior and local log-likelihood caches of a local-grid fit after its
# coefficients changed.
.vcmoe_screen_refresh_posterior <- function(fit) {
  for (grid_id in seq_along(fit$u_grid)) {
    spec <- .vcmoe_dev_parameter_spec(fit, grid_id)
    par <- .vcmoe_dev_parameter_vector(fit, grid_id)
    log_terms <- .vcmoe_dev_family_loglik_matrix(.vcmoe_dev_unpack_parameter(par, spec), spec)
    log_norm <- .row_log_sum_exp(log_terms)
    posterior <- exp(log_terms - log_norm)
    fit$posterior[, , grid_id] <- posterior / rowSums(posterior)
    fit$diagnostics$loglik[[grid_id]] <- sum(spec$weights * log_norm)
    fit$diagnostics$posterior_entropy[[grid_id]] <- .weighted_entropy(fit$posterior[, , grid_id], spec$weights)
  }
  fit
}

# The null fit with x appended as a zero-coefficient covariate on the given
# sides; the reference for checking the closed-form screen moments.
.vcmoe_screen_augment_fit <- function(fit, x, name, sides = "expert") {
  append_slice <- function(a) {
    dn <- dimnames(a)
    out <- array(c(a, numeric(prod(dim(a)[1:2]))), c(dim(a)[1:2], dim(a)[3L] + 1L))
    dn[[3L]] <- c(dn[[3L]], name)
    dimnames(out) <- dn
    out
  }
  if ("expert" %in% sides) {
    z <- fit$fitted$z_design
    fit$fitted$z_design <- cbind(z, x)
    colnames(fit$fitted$z_design) <- c(colnames(z), name)
    fit$coefficients$expert <- append_slice(fit$coefficients$expert)
    fit$coefficients$expert_slope <- append_slice(fit$coefficients$expert_slope)
  }
  if ("gating" %in% sides) {
    z <- fit$fitted$x_design
    fit$fitted$x_design <- cbind(z, x)
    colnames(fit$fitted$x_design) <- c(colnames(z), name)
    fit$coefficients$gating <- append_slice(fit$coefficients$gating)
    fit$coefficients$gating_slope <- append_slice(fit$coefficients$gating_slope)
  }
  fit
}

# Gating-profiled screen. The one-step screen is a Neyman C(alpha)
# step, valid when every nuisance estimate is root-n consistent under H0. A
# gene's own gating coefficient is nuisance, but the null fit holds it at
# zero: for a gene that marks the latent subtypes that is far from the truth,
# and the first-order projection (adjust_gating) cannot absorb an O(1) gating
# effect. Here the gene first enters the gating of the null and every screened
# grid point is refitted by local EM from the null fit: the restricted local
# MLE under H0 (no expert effect of the gene, gating effect free). The
# expert-side screen then runs on that restricted fit, where the gene is
# already in the gating, so only its expert directions are added.
.vcmoe_screen_profile_fit <- function(fit, x, name, grids, maxit = 300L, tol = 1e-9) {
  aug <- .vcmoe_screen_augment_fit(fit, x, name, sides = "gating")
  grid_parameters <- lapply(seq_along(aug$u_grid), function(g) .vcmoe_dev_parameter_vector(aug, g))
  rows <- vector("list", length(grids))
  for (i in seq_along(grids)) {
    g <- grids[[i]]
    spec <- .vcmoe_dev_parameter_spec(aug, g)
    start <- grid_parameters[[g]]
    active <- spec$weights > 0
    # Density-normalised kernel weights do not depend on the row set, so the
    # local EM on the kernel support is the full local EM.
    em_rows <- if (identical(spec$estimation_spec$weight_normalization, "density")) which(active) else seq_along(active)
    em <- .vcmoe_dev_spec_rows(spec, em_rows)
    em_control <- utils::modifyList(em$control, list(maxit = as.integer(maxit), tol = tol))
    local <- .local_em(em$y, em$trials, aug$family, NULL, em$z_design, em$x_design, em$offset,
                       em$u, em$u0, aug$bandwidth, aug$k, em_control, em$estimation_spec,
                       init_params = .vcmoe_dev_unpack_parameter(start, spec))
    par <- .vcmoe_glrt_params_to_parameter_vector(.local_to_params(local), spec)
    names(par) <- names(start)
    grid_parameters[[g]] <- par
    tau0 <- fit$posterior[active, , g, drop = FALSE][, , 1L]
    support <- .vcmoe_dev_spec_rows(spec, which(active))
    terms <- .vcmoe_dev_family_loglik_matrix(.vcmoe_dev_unpack_parameter(par, support), support)
    rows[[i]] <- data.frame(
      # A zero-length EM loop (maxit = 0) leaves no iteration count.
      grid_id = g, iterations = if (length(local$iterations)) as.integer(local$iterations) else 0L,
      converged = isTRUE(local$converged),
      objective_gain = .vcmoe_dev_local_loglik(aug, grid_id = g, par = par, penalized = TRUE) -
        .vcmoe_dev_local_loglik(aug, grid_id = g, par = start, penalized = TRUE),
      map_agreement = stats::weighted.mean(max.col(tau0, ties.method = "first") ==
                                             max.col(terms, ties.method = "first"), spec$weights[active]))
  }
  profiled <- .vcmoe_screen_refresh_posterior(.vcmoe_glrt_apply_grid_parameters(aug, grid_parameters))
  k <- profiled$k
  contrast <- profiled$coefficients$gating[grids, 1L, name] - profiled$coefficients$gating[grids, k, name]
  info <- do.call(rbind, rows)
  info$gating_contrast <- as.numeric(contrast)
  list(fit = profiled, grids = info)
}

#' Gating-profiled gene screen (internal)
#'
#' Per gene: the restricted local MLE with the gene in the null gating
#' (`.vcmoe_screen_profile_fit`), then the expert-side screen on it with the
#' shared Webb weights of `seed`, so draws stay comparable across genes and
#' with `.vcmoe_score_screen`. Genes already in the null gating are screened
#' on the null fit unchanged.
#' @noRd
.vcmoe_score_screen_profiled <- function(fit, X, cluster, grids = seq_along(fit$u_grid), B = 9999L, seed = NULL,
                                         control = list(), center = TRUE, studentize = "onestep",
                                         information = "complete", maxit = 300L, tol = 1e-9, weights = NULL) {
  .vcmoe_cluster_require_fit(fit)
  if (!identical(fit$family, "negative-binomial")) {
    stop("The expert-side screen is implemented for negative-binomial experts.", call. = FALSE)
  }
  X <- as.matrix(X)
  if (!is.numeric(X) || nrow(X) != length(fit$fitted$u) || any(!is.finite(X)) ||
      is.null(colnames(X)) || anyDuplicated(colnames(X))) {
    stop("`X` must be a finite numeric matrix with one row per fitted row and unique column names.", call. = FALSE)
  }
  grids <- sort(unique(as.integer(grids)))
  if (as.integer(B) > 0L && is.null(weights)) {
    weights <- .vcmoe_seeded(seed, .vcmoe_webb_weights(B, max(.vcmoe_cluster_index(fit, cluster))))
  }
  in_gating <- .vcmoe_screen_in_gating(fit, X)
  parts <- profiles <- vector("list", ncol(X))
  for (j in seq_len(ncol(X))) {
    name <- colnames(X)[[j]]
    target <- fit
    if (!in_gating[[j]]) {
      prof <- .vcmoe_screen_profile_fit(fit, X[, j], name, grids, maxit = maxit, tol = tol)
      target <- prof$fit
      profiles[[j]] <- cbind(gene = name, prof$grids, stringsAsFactors = FALSE)
    }
    res <- .vcmoe_score_screen(target, X[, j, drop = FALSE], cluster, grids = grids, side = "expert", B = B,
                               seed = seed, control = control, center = center, studentize = studentize,
                               information = information, weights = weights)
    p <- profiles[[j]]
    res$genes$profiled <- !in_gating[[j]]
    res$genes$profile_max_iterations <- if (is.null(p)) NA_integer_ else max(p$iterations)
    res$genes$profile_all_converged <- if (is.null(p)) NA else all(p$converged)
    res$genes$profile_min_map_agreement <- if (is.null(p)) NA_real_ else min(p$map_agreement)
    res$genes$profile_max_abs_gating <- if (is.null(p)) NA_real_ else max(abs(p$gating_contrast))
    res$genes$profile_objective_gain <- if (is.null(p)) NA_real_ else sum(p$objective_gain)
    parts[[j]] <- res
  }
  first <- parts[[1L]]
  draws <- if (B > 0L) do.call(cbind, lapply(parts, `[[`, "draws")) else NULL
  list(genes = do.call(rbind, lapply(parts, `[[`, "genes")),
       pointwise = do.call(rbind, lapply(parts, `[[`, "pointwise")),
       grids = first$grids, profile = do.call(rbind, profiles), draws = draws, side = "expert", B = B,
       seed = seed, components = first$components, studentize = studentize, information = information,
       adjust_gating = FALSE, profiled = TRUE)
}

#' Cluster-robust gene screen for Negative-Binomial expert effects
#'
#' Screens candidate covariates (typically genes) for an association with the
#' response inside the latent components, while the fitted VCMoE (the null
#' model) is held fixed. For each candidate, every grid point of the null fit
#' is augmented by the covariate, entering each expert's local linear predictor
#' (the tested coefficients) and each gating contrast (nuisance), with
#' coefficients starting at zero. A one-step estimate of the expert
#' coefficients is studentized with the CR2 cluster-robust covariance and
#' calibrated over the grid by a studentized wild cluster bootstrap, as in the
#' path tests of [vcmoe_confband()]. The null model is fitted once, so all
#' candidates share its moments and the screen scales to thousands of genes.
#'
#' @param fit The null model: a single-response Negative-Binomial local-grid
#'   `vcmoe` fit that retains its data.
#' @param X Numeric matrix of candidate covariates, one row per retained fitted
#'   row (in fitted-row order) and one uniquely named column per candidate.
#'   Columns are centred internally.
#' @param cluster Independent-cluster identifiers (for example donors), one per
#'   retained fitted row. Required.
#' @param profile_genes Optional column names of `X` to screen after profiling
#'   their gating effect: the candidate first enters the gating of the null
#'   model and every screened grid point is refitted by local EM, giving the
#'   restricted local fit under the null (no expert effect, gating effect
#'   free). Use this for strong markers of the latent components, whose gating
#'   effect is far from zero; the one-step screen holds it at zero.
#' @param grids Grid indices that form the tested path. Defaults to every grid
#'   point. A candidate is reported with an `NA` p-value when any of these grid
#'   points is blocked for it.
#' @param B Number of wild bootstrap draws (Webb weights). One set of weights
#'   is generated per call and shared by every candidate, chunk and profiled
#'   refit, including when `seed = NULL`. When `6^G <= B` for `G` clusters,
#'   all `6^G` weight vectors are used.
#' @param seed Optional seed for the bootstrap weights.
#' @param adjust Multiplicity adjustment over the full candidate set:
#'   `"BH"` (Benjamini-Hochberg, default), `"minP"` (single-step
#'   Westfall-Young min-p from the shared draws) or `"both"`. Candidates with
#'   an `NA` p-value stay in the family size of the BH adjustment.
#' @param chunk_size Number of candidates processed together; smaller values
#'   use less memory and do not change the results.
#' @param control Optional inference controls passed to the local moment
#'   computations (for example `strict`, `max_hessian_dim`).
#' @param profile_control List with `maxit` (default 300) and `tol` (default
#'   `1e-9`) for the profiled refits.
#' @details The screen tests expert effects, so the null model should explain
#' component membership well: include covariates that distinguish the latent
#' components (for example known subtype markers) in its gating formula.
#' Otherwise strong markers can leak into the expert test. Results describe
#' associations on the log-mean scale of each expert; `effect_per_sd` gives
#' the one-step level coefficient at the maximising grid point per standard
#' deviation of the candidate.
#' @return A list of class `vcmoe_gene_screen` with `genes` (one row per
#'   candidate: sup-t statistic, raw and adjusted p-values, maximising
#'   component and grid, per-component statistics and directions, per-SD
#'   effects, grid usage and status), `pointwise` (per grid and coefficient
#'   one-step estimates, standard errors and Bell-McCaffrey degrees of
#'   freedom), `grids` (grid-level status), `profile` (profiled refits),
#'   `draws` (shared bootstrap maxima) and `settings`.
#' @export
vcmoe_gene_screen <- function(fit, X, cluster, profile_genes = NULL, grids = NULL, B = 9999L, seed = NULL,
                              adjust = c("BH", "minP", "both"), chunk_size = 200L, control = list(),
                              profile_control = list()) {
  if (!inherits(fit, "vcmoe")) stop("`fit` must be a VCMoE fit.", call. = FALSE)
  .vcmoe_cluster_require_fit(fit)
  if (!identical(fit$family, "negative-binomial")) {
    stop("`vcmoe_gene_screen()` supports Negative-Binomial fits.", call. = FALSE)
  }
  if (missing(cluster) || is.null(cluster)) {
    stop("`cluster` is required: give one independent-cluster ID per retained fitted row.", call. = FALSE)
  }
  adjust <- match.arg(adjust)
  X <- as.matrix(X)
  n <- length(fit$fitted$u)
  if (!is.numeric(X) || nrow(X) != n || any(!is.finite(X)) || is.null(colnames(X)) || anyDuplicated(colnames(X))) {
    stop("`X` must be a finite numeric matrix with one row per retained fitted row and unique column names.",
         call. = FALSE)
  }
  if (!is.null(rownames(X))) {
    expected <- rownames(fit$fitted$z_design) %||% as.character(fit$rows_used)
    if (!identical(rownames(X), expected)) {
      stop("Row names of `X` must match the retained fitted rows in order.", call. = FALSE)
    }
  }
  group <- .vcmoe_cluster_index(fit, cluster)
  G <- max(group)
  if (G > n / 2) {
    warning("More than half of the rows are their own cluster; `cluster` should identify independent units such as donors.",
            call. = FALSE)
  }
  B <- as.integer(B)
  .vcmoe_check_multiplier_B(B)
  if (min(B, 6^G) * G > 5e7) {
    stop("The shared bootstrap weight matrix (B x clusters) would exceed 5e7 entries; reduce `B`.", call. = FALSE)
  }
  grids <- sort(unique(as.integer(grids %||% seq_along(fit$u_grid))))
  if (!length(grids) || anyNA(grids) || any(grids < 1L | grids > length(fit$u_grid))) {
    stop(sprintf("`grids` must be grid indices in 1..%d.", length(fit$u_grid)), call. = FALSE)
  }
  profile_genes <- as.character(profile_genes %||% character())
  unknown <- setdiff(profile_genes, colnames(X))
  if (length(unknown)) stop("`profile_genes` not in `X`: ", paste(unknown, collapse = ", "), call. = FALSE)
  chunk_size <- as.integer(chunk_size)
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("`chunk_size` must be a positive integer.", call. = FALSE)
  }
  profile_control <- utils::modifyList(list(maxit = 300L, tol = 1e-9), profile_control %||% list())
  weights <- .vcmoe_seeded(seed, .vcmoe_webb_weights(B, G))
  B_used <- nrow(weights)
  genes <- colnames(X)
  plain <- setdiff(genes, profile_genes)
  parts <- list()
  for (ids in split(plain, ceiling(seq_along(plain) / chunk_size))) {
    parts[[length(parts) + 1L]] <- .vcmoe_score_screen(fit, X[, ids, drop = FALSE], cluster, grids = grids,
                                                       side = "expert", B = B_used, control = control,
                                                       weights = weights)
  }
  profile_rows <- NULL
  if (length(profile_genes)) {
    prof <- .vcmoe_score_screen_profiled(fit, X[, profile_genes, drop = FALSE], cluster, grids = grids,
                                         B = B_used, control = control, maxit = profile_control$maxit,
                                         tol = profile_control$tol, weights = weights)
    parts[[length(parts) + 1L]] <- prof
    profile_rows <- prof$profile
  }
  bind <- function(rows) {
    rows <- Filter(function(x) !is.null(x) && NROW(x), rows)
    if (!length(rows)) return(data.frame())
    cols <- unique(unlist(lapply(rows, names)))
    do.call(rbind, lapply(rows, function(r) {
      r[setdiff(cols, names(r))] <- NA
      r[cols]
    }))
  }
  gene_table <- bind(lapply(parts, `[[`, "genes"))
  gene_table <- gene_table[match(genes, gene_table$gene), , drop = FALSE]
  pointwise <- bind(lapply(parts, `[[`, "pointwise"))
  draws <- do.call(cbind, lapply(parts, `[[`, "draws"))[, genes, drop = FALSE]
  gene_table$profiled <- gene_table$gene %in% profile_genes
  # Status: a candidate whose declared path is incomplete, or whose profiled
  # refit did not converge, keeps its row but gets an NA p-value.
  status <- rep("ok", nrow(gene_table))
  status[gene_table$n_blocked > 0L] <- "incomplete_path"
  if (!is.null(gene_table$profile_all_converged)) {
    status[gene_table$profiled & !gene_table$profile_all_converged %in% TRUE] <- "profile_not_converged"
  }
  gene_table$status <- status
  p_cols <- grep("^(p_value|constancy_p_value)", names(gene_table), value = TRUE)
  gene_table$p_value_raw_path <- gene_table$p_value
  for (col in p_cols) gene_table[[col]][status != "ok"] <- NA_real_
  # Per-SD effects at each component's maximising grid point.
  sds <- apply(X, 2L, stats::sd)
  components <- parts[[1L]]$components
  gene_table$effect_per_sd <- NA_real_
  for (comp in components) {
    col <- paste0("effect_per_sd_", comp)
    gene_table[[col]] <- NA_real_
    for (j in seq_len(nrow(gene_table))) {
      rows <- pointwise[pointwise$gene == gene_table$gene[j] & pointwise$component == comp &
                          pointwise$block == "intercept" & !pointwise$blocked, , drop = FALSE]
      if (!nrow(rows)) next
      best <- which.max(abs(rows$t))
      gene_table[[col]][j] <- rows$estimate[best] * sds[[gene_table$gene[j]]]
      if (identical(comp, gene_table$argmax_component[j])) gene_table$effect_per_sd[j] <- gene_table[[col]][j]
    }
  }
  if (adjust %in% c("BH", "both")) {
    gene_table$q_value <- stats::p.adjust(gene_table$p_value, method = "BH", n = nrow(gene_table))
  }
  if (adjust %in% c("minP", "both")) {
    gene_table$p_minp <- NA_real_
    usable <- is.finite(gene_table$statistic) & colSums(is.na(draws)) == 0
    if (any(usable)) {
      adjusted <- .vcmoe_screen_minp(gene_table$statistic[usable], draws[, usable, drop = FALSE])
      gene_table$p_minp[usable] <- adjusted
      gene_table$p_minp[status != "ok"] <- NA_real_
    }
  }
  rownames(gene_table) <- NULL
  structure(list(
    genes = gene_table, pointwise = pointwise, grids = parts[[1L]]$grids, profile = profile_rows,
    draws = draws,
    settings = list(grids = grids, B = B_used, B_requested = B, enumerated = isTRUE(attr(weights, "enumerated")),
                    seed = seed, n_clusters = G, adjust = adjust, n_candidates = length(genes),
                    profile_genes = profile_genes, components = components, chunk_size = chunk_size)
  ), class = "vcmoe_gene_screen")
}

#' @export
print.vcmoe_gene_screen <- function(x, ...) {
  g <- x$genes
  cat("VCMoE gene screen\n")
  cat("  candidates: ", nrow(g), " (profiled: ", sum(g$profiled), ")\n", sep = "")
  cat("  grid points: ", length(x$settings$grids), "; clusters: ", x$settings$n_clusters,
      "; bootstrap draws: ", x$settings$B, "\n", sep = "")
  cat("  status ok: ", sum(g$status == "ok"), "/", nrow(g), "\n", sep = "")
  if (!is.null(g$q_value)) cat("  BH q <= 0.05: ", sum(g$q_value <= 0.05, na.rm = TRUE), "\n", sep = "")
  if (!is.null(g$p_minp)) cat("  minP <= 0.05: ", sum(g$p_minp <= 0.05, na.rm = TRUE), "\n", sep = "")
  invisible(x)
}
