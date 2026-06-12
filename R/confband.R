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

#' Analytic-style confidence bands for a VCMoE fit
#'
#' @param fit A `vcmoe` fit with `k = 2:10`.
#' @param data Optional original data frame. The current implementation uses
#'   the data stored in `fit$fitted`; refit with `keep_data = TRUE` if needed.
#' @param level Confidence level.
#' @param type Interval columns to expose as `lower` and `upper`.
#' @param coefficient_set Coefficient blocks to return.
#' @param strict Whether weak local fits should return blocked intervals.
#' @param control Optional development inference controls. HC0 is the only
#'   active covariance adjustment.
#' @return A `vcmoe_confband` object with interval and diagnostic data frames.
#' @export
vcmoe_confband <- function(fit, data = NULL, level = 0.95,
                           type = c("pointwise", "simultaneous"),
                           coefficient_set = c("expert", "gating", "sigma", "theta"),
                           strict = TRUE,
                           control = list()) {
  if (!inherits(fit, "vcmoe")) {
    stop("`fit` must be a VCMoE fit.", call. = FALSE)
  }
  type <- match.arg(type)
  coefficient_set <- match.arg(coefficient_set, c("expert", "gating", "sigma", "theta"), several.ok = TRUE)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be a number between 0 and 1.", call. = FALSE)
  }
  domain_length <- fit$u_scaling$domain_length %||%
    vcmoe_parameterization(fit)$u_scaling$domain_length %||% 1
  inference_control <- utils::modifyList(
    list(
      simultaneous_method = "analytic_epanechnikov_path",
      covariance_adjustment = "HC0",
      scb_domain_length = domain_length
    ),
    control %||% list()
  )
  result <- .vcmoe_dev_intervals(
    fit = fit,
    data = data,
    level = level,
    strict = strict,
    control = inference_control
  )
  intervals <- .vcmoe_confband_filter(result$intervals, coefficient_set)
  if (identical(type, "pointwise")) {
    intervals$lower <- intervals$pointwise_lower
    intervals$upper <- intervals$pointwise_upper
  } else {
    intervals$lower <- intervals$simultaneous_lower
    intervals$upper <- intervals$simultaneous_upper
  }
  intervals$type <- type
  rownames(intervals) <- NULL
  out <- list(
    fit = fit,
    intervals = intervals,
    diagnostics = result$diagnostics,
    settings = list(
      family = fit$family,
      k = fit$k,
      level = level,
      type = type,
      coefficient_set = coefficient_set,
      covariance_adjustment = "HC0",
      simultaneous_method = "analytic_epanechnikov_path",
      parameterization = fit$parameterization_id %||% vcmoe_parameterization(fit)$id,
      u_scale = fit$u_scale %||% fit$u_scaling$method %||% "unit",
      domain_length = domain_length
    )
  )
  class(out) <- "vcmoe_confband"
  out
}

#' @export
print.vcmoe_confband <- function(x, ...) {
  ok <- sum(x$intervals$status == "ok", na.rm = TRUE)
  total <- nrow(x$intervals)
  cat("VCMoE analytic-style confidence bands\n")
  cat("  family: ", x$settings$family, "\n", sep = "")
  cat("  components: ", x$settings$k, "\n", sep = "")
  cat("  type: ", x$settings$type, "\n", sep = "")
  cat("  level: ", x$settings$level, "\n", sep = "")
  cat("  interval rows ok: ", ok, "/", total, "\n", sep = "")
  invisible(x)
}
