.vcmoe_parameterization_metadata <- function(family, k, bandwidth, label,
                                             alignment_method, control,
                                             u_scaling = NULL,
                                             parameterization = "a1_epanechnikov_scaled",
                                             estimation_spec = NULL) {
  defaults <- .vcmoe_default_control(control %||% list())
  estimation_spec <- estimation_spec %||% .vcmoe_estimation_spec(parameterization)
  dispersion <- if (identical(family, "gaussian")) {
    "component_specific_local_linear_log_sigma"
  } else if (identical(family, "negative-binomial")) {
    "component_specific_theta"
  } else {
    "none"
  }
  list(
    version = "current_v1",
    id = estimation_spec$id,
    estimator = "local_linear_em",
    family = family,
    k = k,
    bandwidth = bandwidth,
    kernel = list(
      name = estimation_spec$kernel,
      weight_expression = estimation_spec$weight_expression,
      weight_normalization = if (identical(estimation_spec$weight_normalization, "density")) {
        "density_over_bandwidth"
      } else {
        "mean_one_per_grid"
      },
      note = if (identical(estimation_spec$id, "a1_epanechnikov_scaled")) {
        "Epanechnikov weights use K((u-u0)/h)/h."
      } else {
        "Non-default estimator metadata."
      }
    ),
    local_linear_basis = list(
      intercept_column = "1",
      slope_column = estimation_spec$slope_column,
      slope_scale = estimation_spec$slope_scale,
      slope_storage = estimation_spec$slope_storage,
      scaled_basis_conversion = estimation_spec$slope_scaled_conversion
    ),
    u_scaling = list(
      method = u_scaling$method %||% "none",
      original_range = u_scaling$original_range %||% c(NA_real_, NA_real_),
      original_domain_length = u_scaling$original_domain_length %||% NA_real_,
      analysis_range = u_scaling$analysis_range %||% c(NA_real_, NA_real_),
      domain_length = u_scaling$domain_length %||% NA_real_
    ),
    gating = list(
      probability = "softmax_c(X beta_c(u))",
      storage = "centered_logits",
      k2_contrast = "component1 - component2",
      highk_default_contrast = "component_c - selected_baseline"
    ),
    component_structure = "free_component_specific_coefficients",
    dispersion = dispersion,
    label_alignment = list(
      requested = label,
      method = alignment_method
    ),
    optimization = list(
      ridge = defaults$ridge,
      binomial_ridge = defaults$binomial_ridge,
      negbin_ridge = defaults$negbin_ridge,
      negbin_theta_ridge = defaults$negbin_theta_ridge,
      negbin_theta_target = defaults$negbin_theta_target,
      negbin_theta_bounds = c(defaults$negbin_theta_min, defaults$negbin_theta_max),
      tolerance = defaults$tol,
      maxit = defaults$maxit,
      negbin_mstep_maxit = defaults$negbin_mstep_maxit,
      n_starts = defaults$n_starts
    ),
    inference_mode = estimation_spec$inference_mode
  )
}

.vcmoe_require_fit <- function(object) {
  if (!inherits(object, "vcmoe")) {
    stop("`object` must be a VCMoE fit.", call. = FALSE)
  }
  invisible(TRUE)
}

.vcmoe_component_index <- function(object, component) {
  labels <- dimnames(object$coefficients$gating)[[2L]]
  if (is.null(component)) {
    component <- if (object$k == 2L) 2L else 1L
  }
  if (is.character(component)) {
    match_id <- match(component, labels)
    if (is.na(match_id)) {
      stop("Unknown component label in `baseline`.", call. = FALSE)
    }
    return(match_id)
  }
  component <- as.integer(component)
  if (length(component) != 1L || is.na(component) ||
      component < 1L || component > object$k) {
    stop("`baseline` must identify one fitted component.", call. = FALSE)
  }
  component
}

#' Inspect VCMoE parameterization metadata
#'
#' @param object A `vcmoe` object.
#' @return A named list describing the estimator convention used by the fit.
#' @export
vcmoe_parameterization <- function(object) {
  .vcmoe_require_fit(object)
  object$parameterization %||%
    .vcmoe_parameterization_metadata(
      family = object$family,
      k = object$k,
      bandwidth = object$bandwidth,
      label = object$label %||% object$diagnostics$alignment_method %||% "align",
      alignment_method = object$diagnostics$alignment_method %||% object$label %||% "align",
      control = object$control %||% list(),
      u_scaling = object$u_scaling,
      parameterization = object$parameterization_id %||% "a1_epanechnikov_scaled"
    )
}

#' Inspect VCMoE local-linear slopes on the scaled basis
#'
#' VCMoE stores slopes on the scaled local-linear basis `(u - u0) / h`.
#' This helper returns the stored scaled-basis slope block.
#'
#' @param object A `vcmoe` object.
#' @param type Coefficient block, either `"expert"` or `"gating"`.
#' @param bandwidth Optional bandwidth recorded in the returned attributes.
#'   Defaults to the fitted bandwidth.
#' @return An array with the same dimensions as the stored slope block.
#' @export
vcmoe_scaled_slopes <- function(object, type = c("expert", "gating"),
                                bandwidth = NULL) {
  .vcmoe_require_fit(object)
  type <- match.arg(type)
  if (is.null(bandwidth)) {
    bandwidth <- object$bandwidth
  }
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
      !is.finite(bandwidth) || bandwidth <= 0) {
    stop("`bandwidth` must be a positive finite number.", call. = FALSE)
  }
  slopes <- object$coefficients[[paste0(type, "_slope")]]
  metadata <- vcmoe_parameterization(object)
  out <- .slope_to_scaled(slopes, bandwidth, .vcmoe_estimation_spec(metadata$id %||% "a1_epanechnikov_scaled"))
  attr(out, "basis") <- "(u - u0) / bandwidth"
  attr(out, "source_basis") <- metadata$local_linear_basis$slope_column
  attr(out, "bandwidth") <- bandwidth
  out
}

#' Report identifiable gating contrasts
#'
#' @param object A `vcmoe` object.
#' @param baseline Component used as the contrast baseline. By default,
#'   `k = 2` uses component 2 so the contrast is component 1 versus component
#'   2; `k > 2` uses component 1.
#' @param scaled If `TRUE`, slope contrasts are converted to the scaled
#'   local-linear basis `(u - u0) / h`.
#' @return A data frame with one row per grid point, contrast, term, and block.
#' @export
vcmoe_gating_contrasts <- function(object, baseline = NULL, scaled = FALSE) {
  .vcmoe_require_fit(object)
  baseline_id <- .vcmoe_component_index(object, baseline)
  labels <- dimnames(object$coefficients$gating)[[2L]]
  terms <- dimnames(object$coefficients$gating)[[3L]]
  component_ids <- setdiff(seq_len(object$k), baseline_id)
  slope <- object$coefficients$gating_slope
  metadata <- vcmoe_parameterization(object)
  slope_scale <- if (isTRUE(scaled)) "scaled_(u-u0)/bandwidth" else metadata$local_linear_basis$slope_column
  if (isTRUE(scaled)) {
    slope <- vcmoe_scaled_slopes(object, "gating")
  }

  rows <- vector(
    "list",
    length(object$u_grid) * length(component_ids) * length(terms) * 2L
  )
  row_id <- 1L
  for (grid_id in seq_along(object$u_grid)) {
    for (component_id in component_ids) {
      for (term_id in seq_along(terms)) {
        rows[[row_id]] <- data.frame(
          u = object$u_grid[[grid_id]],
          component = labels[[component_id]],
          baseline = labels[[baseline_id]],
          contrast = paste(labels[[component_id]], "vs", labels[[baseline_id]], sep = "_"),
          term = terms[[term_id]],
          block = "intercept",
          estimate = object$coefficients$gating[grid_id, component_id, term_id] -
            object$coefficients$gating[grid_id, baseline_id, term_id],
          slope_scale = NA_character_,
          stringsAsFactors = FALSE
        )
        row_id <- row_id + 1L
        rows[[row_id]] <- data.frame(
          u = object$u_grid[[grid_id]],
          component = labels[[component_id]],
          baseline = labels[[baseline_id]],
          contrast = paste(labels[[component_id]], "vs", labels[[baseline_id]], sep = "_"),
          term = terms[[term_id]],
          block = "slope",
          estimate = slope[grid_id, component_id, term_id] -
            slope[grid_id, baseline_id, term_id],
          slope_scale = slope_scale,
          stringsAsFactors = FALSE
        )
        row_id <- row_id + 1L
      }
    }
  }
  do.call(rbind, rows)
}
