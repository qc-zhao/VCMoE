.vcmoe_confband_filter <- function(intervals, coefficient_set) {
  keep <- rep(FALSE, nrow(intervals))
  if ("expert" %in% coefficient_set) {
    keep <- keep | intervals$coefficient_set == "expert"
  }
  if ("gating" %in% coefficient_set) {
    keep <- keep | intervals$coefficient_set == "gating_contrast"
  }
  if ("sigma" %in% coefficient_set) {
    keep <- keep | intervals$coefficient_set == "nuisance" & intervals$term == "log_sigma"
  }
  if ("theta" %in% coefficient_set) {
    keep <- keep | intervals$coefficient_set == "nuisance" & intervals$term == "log_theta"
  }
  intervals[keep, , drop = FALSE]
}

#' Confidence bands for a VCMoE fit
#'
#' Pointwise intervals and finite-grid simultaneous bands for the coefficient
#' functions of a fitted VCMoE model, from the local sandwich covariance.
#'
#' @param fit A `vcmoe` fit with `k = 2:10`.
#' @param data Optional original data frame. The current implementation uses
#'   the data stored in `fit$fitted`; refit with `keep_data = TRUE` if needed.
#' @param level Confidence level.
#' @param type Interval columns to expose as `lower` and `upper`.
#' @param coefficient_set Coefficient blocks to return.
#' @param strict Whether weak local fits should return blocked intervals.
#' @param control Optional inference controls. Without `cluster`, the
#'   covariance is HC0 and the simultaneous band uses the analytic
#'   Epanechnikov path critical value; `simultaneous_method =
#'   "multiplier_coefficient_path"` calibrates each coefficient over the
#'   complete fitted grid by a Gaussian multiplier bootstrap instead. With
#'   `cluster`, `covariance_adjustment` may be `"CR0"`, `"CR1"`, `"CR2"`
#'   (default) or `"CR3"`, and `simultaneous_method` `"wild_cluster_path"`
#'   (default) or `"multiplier_coefficient_path"`. `multiplier_B` (default
#'   199) sets the number of bootstrap draws and `seed` their random seed.
#'   `path_tests = TRUE` adds the zero and constancy path tests described
#'   below; `path_sets`, a named list of grid-index vectors, adds bands and
#'   tests over sub-paths (for example the interior grid points).
#' @param cluster Optional atomic vector of independent-cluster identifiers
#'   (for example donors), one per retained fitted row, in fitted-row order
#'   (not the original unfiltered order). If named, names must match the
#'   retained design-matrix row names exactly. Missing IDs are rejected.
#' @details
#' **Independent observations.** HC0 sandwich intervals with normal critical
#' values. For `engine = "joint_path_em"`, the covariance is the observed
#' local-likelihood sandwich evaluated at the joint-path estimate; it does not
#' include shared-path, label-selection, or finite-grid cross-grid
#' responsibility uncertainty, so joint-path convergence and the returned
#' score-imbalance diagnostics should be inspected.
#'
#' **Clustered observations.** With `cluster`, kernel-weighted scores are
#' summed within clusters. Cluster inference supports unrestricted local-grid
#' fits that retain their data. Coefficients and the penalized bread are
#' unchanged. CR0 uses normal critical values. CR1 scales by `G / (G - 1)` with
#' `t(G - 1)`, for `G` clusters with positive kernel weight at the grid point.
#' CR2 (bias-reduced linearization in score space) and CR3 (one-step cluster
#' jackknife) use Bell-McCaffrey degrees of freedom under a working model of
#' independent observations within clusters; both use cluster Hessian
#' blocks. The default simultaneous band studentizes Webb-weighted wild
#' cluster score draws with the same adjustment at every grid point; one
#' weight per cluster is shared by all grid points and coefficients. When the
#' number of distinct Webb weight vectors, `6^G`, is at most `multiplier_B`,
#' all of them are used. Clusters must be independent; smoothing bias and
#' bandwidth selection are not corrected. With fewer than 20 clusters, or
#' fewer than 10 effective clusters, at a grid point a warning flags the
#' intervals as sensitivity analyses.
#'
#' **Paths and tests.** A band covers one coefficient's complete fitted grid,
#' not a continuous domain or all coefficients jointly. A blocked grid point
#' blocks that coefficient's complete-path band, not its other pointwise
#' intervals; a band over a sub-path in `path_sets` remains available. With
#' `path_tests = TRUE` each complete path also gets two sup-t tests from the
#' same bootstrap draws as its band: `beta(u) = 0` (the band-inversion test)
#' and `beta(u)` constant, which contrasts each grid estimate with the path
#' mean and studentizes the contrast with the same small-cluster adjustment.
#' These tests are the recommended constancy tests for clustered data; the
#' likelihood-ratio test of [vcmoe_glrt()] assumes independent observations.
#' @return A `vcmoe_confband` object with interval and diagnostic data frames;
#'   with `cluster`, `cluster_support` reports per-component cluster support,
#'   and with `path_tests = TRUE`, `tests` holds the path-test statistics and
#'   p-values.
#' @export
vcmoe_confband <- function(fit, data = NULL, level = 0.95,
                           type = c("pointwise", "simultaneous"),
                           coefficient_set = c("expert", "gating", "sigma", "theta"),
                           strict = TRUE,
                           control = list(), cluster = NULL) {
  if (!inherits(fit, "vcmoe")) {
    stop("`fit` must be a VCMoE fit.", call. = FALSE)
  }
  type <- match.arg(type)
  coefficient_set <- match.arg(coefficient_set, c("expert", "gating", "sigma", "theta"), several.ok = TRUE)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be a number between 0 and 1.", call. = FALSE)
  }
  if (!is.null(cluster)) .vcmoe_cluster_require_fit(fit)
  joint_path <- identical(fit$engine_id %||% "local_grid_em", "joint_path_em")
  if (joint_path) {
    warning(
      paste0(
        "Joint-path analytic-style bands use the observed local-likelihood ",
        "sandwich plug-in and exclude shared-path, label-selection, and finite-grid ",
        "cross-grid responsibility uncertainty; inspect convergence and ",
        "`score_imbalance_max`."
      ),
      call. = FALSE
    )
  }
  domain_length <- fit$u_scaling$domain_length %||%
    vcmoe_parameterization(fit)$u_scaling$domain_length %||% 1
  inference_control <- utils::modifyList(
    list(
      simultaneous_method = if (is.null(cluster)) "analytic_epanechnikov_path" else "wild_cluster_path",
      covariance_adjustment = if (is.null(cluster)) "HC0" else "CR2",
      scb_domain_length = domain_length
    ),
    control %||% list()
  )
  inference_control$simultaneous_method <- .vcmoe_dev_simultaneous_method(inference_control$simultaneous_method)
  if (is.null(cluster) && toupper(inference_control$covariance_adjustment) %in% c("CR0", "CR1", "CR2", "CR3")) {
    stop(toupper(inference_control$covariance_adjustment),
         " requires explicit `cluster` IDs in retained fitted-row order.", call. = FALSE)
  }
  if (is.null(cluster) && identical(inference_control$simultaneous_method, "wild_cluster_path")) {
    stop("`wild_cluster_path` requires `cluster`.", call. = FALSE)
  }
  if (!is.null(cluster)) {
    .vcmoe_cluster_method(inference_control, is_cluster = TRUE)
    result <- .vcmoe_cluster_intervals(fit, cluster = cluster, level = level,
                                       strict = strict, control = inference_control)
  } else {
    result <- .vcmoe_dev_intervals(
      fit = fit,
      data = data,
      level = level,
      strict = strict,
      control = inference_control
    )
  }
  .vcmoe_confband_pack(fit, result, level, type, coefficient_set, inference_control, domain_length)
}

.vcmoe_confband_pack <- function(fit, result, level, type, coefficient_set, control, domain_length) {
  joint_path <- identical(fit$engine_id %||% "local_grid_em", "joint_path_em")
  intervals <- .vcmoe_confband_filter(result$intervals, coefficient_set)
  if (identical(type, "pointwise")) {
    intervals$lower <- intervals$pointwise_lower
    intervals$upper <- intervals$pointwise_upper
  } else {
    intervals$lower <- intervals$simultaneous_lower
    intervals$upper <- intervals$simultaneous_upper
  }
  if ("pointwise_status" %in% names(intervals)) {
    intervals$status <- intervals[[paste0(type, "_status")]]
    intervals$block_reason <- intervals[[paste0(type, "_block_reason")]]
  }
  intervals$type <- type
  rownames(intervals) <- NULL
  out <- list(
    fit = fit,
    intervals = intervals,
    diagnostics = result$diagnostics,
    cluster_support = result$cluster_support,
    settings = list(
      family = fit$family,
      k = fit$k,
      engine_id = fit$engine_id %||% "local_grid_em",
      estimating_equation = if (joint_path) {
        "joint_path_observed_local_likelihood_plugin"
      } else {
        "observed_local_likelihood"
      },
      covariance_scope = if (joint_path) {
        "local_asymptotic_plugin_excludes_finite_grid_cross_grid_coupling"
      } else {
        "local_sandwich"
      },
      covariance_target = "observed_local_marginal_likelihood",
      estimator_covariance_match = if (joint_path) {
        "asymptotic_plugin_not_exact_finite_grid_match"
      } else {
        "matched_when_local_fit_converged_and_penalty_negligible"
      },
      shared_path_uncertainty_accounted = FALSE,
      label_uncertainty_accounted = FALSE,
      coverage_theory = if (joint_path) {
        "local_likelihood_asymptotic_plugin_not_finite_grid_joint_path_theory"
      } else {
        "local_likelihood_asymptotic_with_bias_and_boundary_caveats"
      },
      level = level,
      type = type,
      coefficient_set = coefficient_set,
      covariance_adjustment = unique(result$intervals$covariance_adjustment),
      simultaneous_method = unique(result$intervals$simultaneous_method),
      multiplier_B = if (control$simultaneous_method %in% c("multiplier_coefficient_path", "wild_cluster_path")) {
        control$multiplier_B %||% 199L
      } else {
        NULL
      },
      seed = control$seed,
      parameterization = fit$parameterization_id %||% vcmoe_parameterization(fit)$id,
      u_scale = fit$u_scale %||% fit$u_scaling$method %||% "unit",
      domain_length = domain_length
    )
  )
  if (!is.null(result$tests)) out$tests <- .vcmoe_confband_filter(result$tests, coefficient_set)
  class(out) <- "vcmoe_confband"
  out
}

#' @export
print.vcmoe_confband <- function(x, ...) {
  ok <- sum(x$intervals$status == "ok", na.rm = TRUE)
  total <- nrow(x$intervals)
  cat("VCMoE confidence bands\n")
  cat("  family: ", x$settings$family, "\n", sep = "")
  cat("  components: ", x$settings$k, "\n", sep = "")
  cat("  engine: ", x$settings$engine_id, "\n", sep = "")
  cat("  type: ", x$settings$type, "\n", sep = "")
  cat("  covariance: ", paste(x$settings$covariance_adjustment, collapse = ","), "\n", sep = "")
  cat("  level: ", x$settings$level, "\n", sep = "")
  cat("  coverage theory: ", x$settings$coverage_theory, "\n", sep = "")
  cat("  interval rows ok: ", ok, "/", total, "\n", sep = "")
  warning_rows <- sum(x$intervals$status == "warning", na.rm = TRUE)
  if (warning_rows) {
    cat("  interval rows warning: ", warning_rows, "/", total, "\n", sep = "")
  }
  invisible(x)
}
