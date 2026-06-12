.vcmoe_glrt_epanechnikov_kernel <- function(t) {
  0.75 * pmax(0, 1 - t^2)
}

.vcmoe_glrt_epanechnikov_convolution <- function(t) {
  f <- function(s) .vcmoe_glrt_epanechnikov_kernel(s) *
    .vcmoe_glrt_epanechnikov_kernel(t - s)
  lower <- max(-1, t - 1)
  upper <- min(1, t + 1)
  if (lower >= upper) {
    return(0)
  }
  stats::integrate(f, lower, upper, rel.tol = 1e-8)$value
}

.vcmoe_glrt_rk_delta_epanechnikov <- function(bandwidth, p_test = 1,
                                              C_const = 2, U_length = 1) {
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
      !is.finite(bandwidth) || bandwidth <= 0) {
    stop("`bandwidth` must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(p_test) || length(p_test) != 1L ||
      !is.finite(p_test) || p_test <= 0) {
    stop("`p_test` must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(C_const) || length(C_const) != 1L ||
      !is.finite(C_const) || C_const <= 0) {
    stop("`C_const` must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(U_length) || length(U_length) != 1L ||
      !is.finite(U_length) || U_length <= 0) {
    stop("`U_length` must be a positive finite number.", call. = FALSE)
  }
  int_K2 <- stats::integrate(
    function(t) .vcmoe_glrt_epanechnikov_kernel(t)^2,
    -1, 1, rel.tol = 1e-8
  )$value
  K0 <- .vcmoe_glrt_epanechnikov_kernel(0)
  g <- function(t) {
    .vcmoe_glrt_epanechnikov_kernel(t) -
      0.5 * vapply(t, .vcmoe_glrt_epanechnikov_convolution, numeric(1L))
  }
  denom <- stats::integrate(function(t) g(t)^2, -2, 2, rel.tol = 1e-8)$value
  num <- K0 - 0.5 * int_K2
  rK <- num / denom
  delta <- rK * p_test * C_const * U_length * num / bandwidth
  list(
    rK = as.numeric(rK),
    delta = as.numeric(delta),
    p_test = as.numeric(p_test),
    C_const = as.numeric(C_const),
    U_length = as.numeric(U_length),
    int_K2 = as.numeric(int_K2),
    K0 = as.numeric(K0)
  )
}

.vcmoe_glrt_control <- function(control) {
  defaults <- list(
    maxit = 200L,
    reltol = 1e-7,
    lambda_tol = 1e-7,
    strict = TRUE,
    min_converged_fraction = 0.80,
    min_component_proportion = 0.05,
    C_const = 2
  )
  utils::modifyList(defaults, control %||% list())
}

.vcmoe_glrt_normalize_calibration <- function(calibration) {
  calibration <- match.arg(
    calibration,
    c("analytic_epanechnikov", "bootstrap", "both", "none", "parametric_bootstrap")
  )
  if (identical(calibration, "parametric_bootstrap")) {
    return("bootstrap")
  }
  calibration
}

.vcmoe_glrt_public_set <- function(coefficient_set) {
  coefficient_set <- match.arg(coefficient_set, c("expert", "gating", "sigma", "theta"))
  coefficient_set
}

.vcmoe_glrt_select_term <- function(terms, term, label) {
  if (is.null(term)) {
    if (length(terms) == 1L) {
      return(terms[[1L]])
    }
    stop("`term` must be provided for coefficient-specific ", label,
         " tests. Available terms: ", paste(terms, collapse = ", "),
         call. = FALSE)
  }
  if (!is.character(term) || length(term) != 1L || !term %in% terms) {
    stop("Unknown `term` for ", label, " test. Available terms: ",
         paste(terms, collapse = ", "), call. = FALSE)
  }
  term
}

.vcmoe_glrt_select_component <- function(labels, component, label) {
  if (is.null(component)) {
    return(labels[[1L]])
  }
  if (is.character(component) && length(component) == 1L) {
    if (!component %in% labels) {
      stop("Unknown `component` for ", label, " test. Available components: ",
           paste(labels, collapse = ", "), call. = FALSE)
    }
    return(component)
  }
  component <- as.integer(component)
  if (length(component) != 1L || is.na(component) ||
      component < 1L || component > length(labels)) {
    stop("`component` must identify one ", label, " component.", call. = FALSE)
  }
  labels[[component]]
}

.vcmoe_glrt_select_gating_component <- function(fit, gating_labels, component) {
  if (is.null(component)) {
    return(gating_labels[[1L]])
  }
  component_labels <- dimnames(fit$coefficients$gating)[[2L]]
  baseline <- if (fit$k == 2L) component_labels[[2L]] else component_labels[[1L]]
  if (is.character(component) && length(component) == 1L) {
    if (component %in% gating_labels) {
      return(component)
    }
    if (component %in% component_labels) {
      if (identical(component, baseline) && fit$k > 2L) {
        stop("For `k > 2` gating tests, component 1 is the baseline; choose ",
             "a non-baseline component or a contrast label. Available contrasts: ",
             paste(gating_labels, collapse = ", "), call. = FALSE)
      }
      if (fit$k == 2L) {
        return(gating_labels[[1L]])
      }
      contrast <- paste(component, baseline, sep = "_vs_")
      if (contrast %in% gating_labels) {
        return(contrast)
      }
    }
    stop("Unknown `component` for gating test. Available components: ",
         paste(component_labels, collapse = ", "), "; available contrasts: ",
         paste(gating_labels, collapse = ", "), call. = FALSE)
  }
  component_id <- as.integer(component)
  if (length(component_id) != 1L || is.na(component_id) ||
      component_id < 1L || component_id > length(component_labels)) {
    stop("`component` must identify one gating component or contrast.",
         call. = FALSE)
  }
  if (fit$k == 2L) {
    return(gating_labels[[1L]])
  }
  if (component_id == 1L) {
    stop("For `k > 2` gating tests, component 1 is the baseline; choose ",
         "a non-baseline component or a contrast label. Available contrasts: ",
         paste(gating_labels, collapse = ", "), call. = FALSE)
  }
  contrast <- paste(component_labels[[component_id]], baseline, sep = "_vs_")
  if (!contrast %in% gating_labels) {
    stop("Selected gating contrast was not found in the fitted parameter table.",
         call. = FALSE)
  }
  contrast
}

.vcmoe_glrt_coefficient_selection <- function(fit, coefficient_set, component, term) {
  coefficient_set <- .vcmoe_glrt_public_set(coefficient_set)
  spec <- .vcmoe_dev_parameter_spec(fit, 1L)
  component_labels <- dimnames(fit$coefficients$expert)[[2L]]
  expert_terms <- dimnames(fit$coefficients$expert)[[3L]]
  gating_terms <- dimnames(fit$coefficients$gating)[[3L]]
  gating_labels <- unique(spec$parameter_table$component[
    spec$parameter_table$coefficient_set == "gating_contrast"
  ])

  if (identical(coefficient_set, "expert")) {
    return(list(
      public_set = "expert",
      internal_set = "expert",
      component = .vcmoe_glrt_select_component(component_labels, component, "expert"),
      term = .vcmoe_glrt_select_term(expert_terms, term, "expert")
    ))
  }
  if (identical(coefficient_set, "gating")) {
    selected_component <- .vcmoe_glrt_select_gating_component(
      fit,
      gating_labels,
      component
    )
    return(list(
      public_set = "gating",
      internal_set = "gating_contrast",
      component = selected_component,
      term = .vcmoe_glrt_select_term(gating_terms, term, "gating")
    ))
  }
  if (identical(coefficient_set, "sigma")) {
    if (!identical(fit$family, "gaussian")) {
      stop("`coefficient_set = \"sigma\"` is only available for Gaussian fits.",
           call. = FALSE)
    }
    return(list(
      public_set = "sigma",
      internal_set = "nuisance",
      component = .vcmoe_glrt_select_component(component_labels, component, "sigma"),
      term = .vcmoe_glrt_select_term("log_sigma", term, "sigma")
    ))
  }
  if (!identical(fit$family, "negative-binomial")) {
    stop("`coefficient_set = \"theta\"` is only available for Negative-Binomial fits.",
         call. = FALSE)
  }
  list(
    public_set = "theta",
    internal_set = "nuisance",
    component = .vcmoe_glrt_select_component(component_labels, component, "theta"),
    term = .vcmoe_glrt_select_term("log_theta", term, "theta")
  )
}

.vcmoe_glrt_constraint <- function(fit, test, coefficient_set, component, term) {
  test <- match.arg(test, c("coefficient", "constant_all"))
  table <- .vcmoe_dev_parameter_spec(fit, 1L)$parameter_table
  if (identical(test, "constant_all")) {
    shared <- table$parameter[table$block %in% c("intercept", "log_theta")]
    fixed <- table$parameter[table$block == "slope"]
    return(list(
      test = "constant_all",
      coefficient_set = "all",
      component = NA_character_,
      term = NA_character_,
      shared_parameters = shared,
      fixed_zero_parameters = fixed,
      p_test = length(shared),
      label = "all_coefficient_functions"
    ))
  }

  selected <- .vcmoe_glrt_coefficient_selection(fit, coefficient_set, component, term)
  mask <- table$coefficient_set == selected$internal_set &
    table$component == selected$component &
    table$term == selected$term
  if (!any(mask)) {
    stop("Selected coefficient path was not found in the fitted parameter table.",
         call. = FALSE)
  }
  shared <- table$parameter[mask & table$block %in% c("intercept", "log_theta")]
  fixed <- table$parameter[mask & table$block == "slope"]
  if (!length(shared)) {
    stop("Selected coefficient path has no intercept/log-theta parameter to constrain.",
         call. = FALSE)
  }
  list(
    test = "coefficient",
    coefficient_set = selected$public_set,
    internal_set = selected$internal_set,
    component = selected$component,
    term = selected$term,
    shared_parameters = shared,
    fixed_zero_parameters = fixed,
    p_test = 1,
    label = paste(selected$public_set, selected$component, selected$term, sep = ":")
  )
}

.vcmoe_glrt_packing <- function(fit, constraint) {
  grid_pars <- lapply(seq_along(fit$u_grid), function(grid_id) {
    .vcmoe_dev_parameter_vector(fit, grid_id)
  })
  shared <- unique(constraint$shared_parameters)
  fixed <- unique(constraint$fixed_zero_parameters)
  shared_start <- vapply(shared, function(name) {
    mean(vapply(grid_pars, function(par) unname(par[[name]]), numeric(1L)), na.rm = TRUE)
  }, numeric(1L))
  names(shared_start) <- paste0("shared:", shared)

  free_entries <- data.frame()
  free_start <- numeric(0L)
  for (grid_id in seq_along(grid_pars)) {
    par <- grid_pars[[grid_id]]
    free_names <- setdiff(names(par), c(shared, fixed))
    if (length(free_names)) {
      free_entries <- rbind(
        free_entries,
        data.frame(
          grid_id = grid_id,
          parameter = free_names,
          stringsAsFactors = FALSE
        )
      )
      values <- par[free_names]
      names(values) <- paste0("grid", grid_id, ":", free_names)
      free_start <- c(free_start, values)
    }
  }

  list(
    fit = fit,
    constraint = constraint,
    grid_pars = grid_pars,
    shared = shared,
    fixed = fixed,
    free_entries = free_entries,
    start = c(shared_start, free_start)
  )
}

.vcmoe_glrt_unpack_grid_parameters <- function(par, packing) {
  n_shared <- length(packing$shared)
  shared_values <- if (n_shared) par[seq_len(n_shared)] else numeric(0L)
  names(shared_values) <- packing$shared
  free_values <- if (length(par) > n_shared) par[-seq_len(n_shared)] else numeric(0L)

  out <- packing$grid_pars
  for (grid_id in seq_along(out)) {
    grid_par <- out[[grid_id]]
    if (length(shared_values)) {
      grid_par[names(shared_values)] <- shared_values
    }
    if (length(packing$fixed)) {
      grid_par[packing$fixed] <- 0
    }
    out[[grid_id]] <- grid_par
  }
  if (length(free_values)) {
    for (row_id in seq_len(nrow(packing$free_entries))) {
      entry <- packing$free_entries[row_id, , drop = FALSE]
      out[[entry$grid_id]][entry$parameter] <- free_values[[row_id]]
    }
  }
  out
}

.vcmoe_glrt_apply_grid_parameters <- function(fit, grid_parameters) {
  null_fit <- fit
  for (grid_id in seq_along(grid_parameters)) {
    spec <- .vcmoe_dev_parameter_spec(null_fit, grid_id)
    params <- .vcmoe_dev_unpack_parameter(grid_parameters[[grid_id]], spec)
    null_fit$coefficients$expert[grid_id, , ] <- params$expert_coef[, seq_len(spec$q), drop = FALSE]
    null_fit$coefficients$expert_slope[grid_id, , ] <- params$expert_coef[, spec$q + seq_len(spec$q), drop = FALSE]
    null_fit$coefficients$gating[grid_id, , ] <- params$gating_coef[, seq_len(spec$p), drop = FALSE]
    null_fit$coefficients$gating_slope[grid_id, , ] <- params$gating_coef[, spec$p + seq_len(spec$p), drop = FALSE]
    if (!is.null(null_fit$coefficients$sigma)) {
      null_fit$coefficients$sigma[grid_id, ] <- params$sigma
    }
    if (!is.null(null_fit$coefficients$sigma_slope)) {
      null_fit$coefficients$sigma_slope[grid_id, ] <- params$sigma_slope
    }
    if (!is.null(null_fit$coefficients$theta)) {
      null_fit$coefficients$theta[grid_id, ] <- params$theta
    }
  }
  null_fit
}

.vcmoe_glrt_null_fit <- function(fit, constraint, control = list()) {
  .vcmoe_dev_require_fit(fit)
  control <- .vcmoe_glrt_control(control)
  packing <- .vcmoe_glrt_packing(fit, constraint)
  objective <- function(par) {
    grid_parameters <- .vcmoe_glrt_unpack_grid_parameters(par, packing)
    value <- tryCatch(sum(vapply(seq_along(grid_parameters), function(grid_id) {
      .vcmoe_dev_local_loglik(
        fit,
        grid_id = grid_id,
        par = grid_parameters[[grid_id]],
        penalized = FALSE
      )
    }, numeric(1L))), error = function(e) NA_real_)
    if (!is.finite(value)) {
      return(1e100)
    }
    -value
  }
  opt <- tryCatch(
    stats::optim(
      packing$start,
      objective,
      method = "BFGS",
      control = list(maxit = as.integer(control$maxit), reltol = control$reltol)
    ),
    error = function(e) e
  )
  if (inherits(opt, "error")) {
    null_fit <- fit
    null_fit$diagnostics$glrt_null_optim_convergence <- NA_integer_
    null_fit$diagnostics$glrt_null_optim_value <- NA_real_
    null_fit$diagnostics$glrt_null_optim_error <- conditionMessage(opt)
    return(null_fit)
  }
  grid_parameters <- .vcmoe_glrt_unpack_grid_parameters(opt$par, packing)
  null_fit <- .vcmoe_glrt_apply_grid_parameters(fit, grid_parameters)
  null_fit$diagnostics$glrt_null_optim_convergence <- opt$convergence
  null_fit$diagnostics$glrt_null_optim_value <- opt$value
  null_fit$diagnostics$glrt_null_optim_error <- NA_character_
  null_fit$diagnostics$glrt_constraint <- constraint$label
  null_fit
}

.vcmoe_glrt_global_pseudologlik <- function(fit) {
  .vcmoe_dev_global_pseudologlik(fit)
}

.vcmoe_glrt_statistic <- function(full_fit, null_fit, lambda_tol = 1e-7) {
  full_loglik <- .vcmoe_glrt_global_pseudologlik(full_fit)
  null_loglik <- .vcmoe_glrt_global_pseudologlik(null_fit)
  lambda <- full_loglik - null_loglik
  block_reason <- NA_character_
  if (!is.finite(lambda)) {
    block_reason <- "nonfinite_lambda"
  } else if (lambda < -lambda_tol) {
    block_reason <- "null_loglik_exceeds_full_loglik"
  }
  lambda_positive <- max(0, lambda)
  list(
    lambda = lambda,
    lambda_positive = lambda_positive,
    lrt_statistic = 2 * lambda_positive,
    full_loglik = full_loglik,
    null_loglik = null_loglik,
    block_reason = block_reason
  )
}

.vcmoe_glrt_fit_warnings <- function(fit, control) {
  diagnostics <- tryCatch(vcmoe_diagnostics(fit), error = function(e) NULL)
  reasons <- character(0L)
  if (mean(fit$diagnostics$converged %in% TRUE) < control$min_converged_fraction) {
    reasons <- c(reasons, "low_converged_fraction")
  }
  if (any(fit$diagnostics$ambiguous %||% FALSE, na.rm = TRUE)) {
    reasons <- c(reasons, "ambiguous_label")
  }
  if (!is.null(diagnostics) &&
      any(diagnostics$min_posterior_mean < control$min_component_proportion, na.rm = TRUE)) {
    reasons <- c(reasons, "component_collapse")
  }
  unique(reasons)
}

.vcmoe_glrt_require_analytic_fit <- function(fit) {
  .vcmoe_dev_require_fit(fit)
  metadata <- vcmoe_parameterization(fit)
  if (!identical(metadata$id, "a1_epanechnikov_scaled")) {
    stop("Analytic Epanechnikov GLRT requires `parameterization = \"a1_epanechnikov_scaled\"`.",
         call. = FALSE)
  }
  if (!identical(metadata$kernel$name, "epanechnikov") ||
      !identical(metadata$kernel$weight_normalization, "density_over_bandwidth") ||
      !identical(metadata$local_linear_basis$slope_storage, "scaled")) {
    stop("Analytic Epanechnikov GLRT requires Epanechnikov density weights and scaled local-linear slopes.",
         call. = FALSE)
  }
  method <- fit$u_scale %||% fit$u_scaling$method %||% "unit"
  if (!identical(method, "unit")) {
    stop("Analytic Epanechnikov GLRT requires `u_scale = \"unit\"`.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.vcmoe_glrt_analytic_calibration <- function(fit, constraint, observed, control) {
  .vcmoe_glrt_require_analytic_fit(fit)
  domain_length <- fit$u_scaling$domain_length %||%
    vcmoe_parameterization(fit)$u_scaling$domain_length %||% 1
  constants <- .vcmoe_glrt_rk_delta_epanechnikov(
    bandwidth = fit$bandwidth,
    p_test = constraint$p_test,
    C_const = control$C_const,
    U_length = domain_length
  )
  analytic_statistic <- constants$rK * observed$lambda_positive
  list(
    analytic_statistic = analytic_statistic,
    analytic_p_value = stats::pchisq(analytic_statistic, df = constants$delta, lower.tail = FALSE),
    rK = constants$rK,
    delta = constants$delta,
    p_test = constants$p_test,
    C_const = constants$C_const,
    U_length = constants$U_length,
    int_K2 = constants$int_K2,
    K0 = constants$K0
  )
}

.vcmoe_glrt_bootstrap <- function(fit, data, null_fit, constraint, B, seed,
                                  control, refit_control, verbose) {
  base_data <- .bootstrap_base_data(fit, data)
  u_info <- .bootstrap_u_info(fit, data, base_data, u = NULL)
  replicate_rows <- vector("list", B)
  successful_stats <- numeric(0L)

  for (replicate_id in seq_len(B)) {
    if (!is.null(seed)) {
      set.seed(as.integer(seed + replicate_id))
    }
    if (isTRUE(verbose)) {
      message("GLRT bootstrap replicate ", replicate_id, " / ", B)
    }
    start_time <- proc.time()[["elapsed"]]
    result <- tryCatch({
      boot_data <- .simulate_bootstrap_response(null_fit, base_data, u_info$values)
      boot_fit <- suppressWarnings(vcmoe_fit(
        formula = fit$formula,
        data = boot_data,
        u = u_info$refit,
        k = fit$k,
        family = fit$family,
        bandwidth = fit$bandwidth,
        u_grid = fit$u_grid,
        control = .bootstrap_refit_control(fit$control, refit_control, seed, replicate_id),
        label = .vcmoe_refit_label(fit),
        u_scale = fit$u_scale %||% fit$u_scaling$method %||% "unit",
        parameterization = fit$parameterization_id %||% vcmoe_parameterization(fit)$id %||% "a1_epanechnikov_scaled"
      ))
      boot_null <- .vcmoe_glrt_null_fit(boot_fit, constraint, control = control)
      .vcmoe_glrt_statistic(boot_fit, boot_null, .vcmoe_glrt_control(control)$lambda_tol)
    }, error = function(e) e)
    runtime <- proc.time()[["elapsed"]] - start_time
    if (inherits(result, "error")) {
      replicate_rows[[replicate_id]] <- data.frame(
        replicate = replicate_id,
        status = "failed",
        lambda = NA_real_,
        lrt_statistic = NA_real_,
        runtime_seconds = runtime,
        error_message = conditionMessage(result),
        stringsAsFactors = FALSE
      )
    } else {
      successful_stats <- c(successful_stats, result$lambda_positive)
      replicate_rows[[replicate_id]] <- data.frame(
        replicate = replicate_id,
        status = "ok",
        lambda = result$lambda,
        lrt_statistic = result$lrt_statistic,
        runtime_seconds = runtime,
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    replicate_summary = do.call(rbind, replicate_rows),
    successful_lambdas = successful_stats
  )
}

#' Generalized likelihood-ratio test for VCMoE coefficient variation
#'
#' @param fit A `vcmoe` fit.
#' @param data Original data frame used to fit `fit`.
#' @param test Test type. `"coefficient"` tests one coefficient function;
#'   `"constant_all"` tests all fitted coefficient functions jointly.
#' @param coefficient_set Coefficient block for coefficient-specific tests.
#' @param component Component label or index for coefficient-specific tests.
#' @param term Term name for coefficient-specific tests.
#' @param calibration Calibration method. `"analytic_epanechnikov"` uses the
#'   Epanechnikov modified chi-square calibration; `"bootstrap"` uses
#'   parametric bootstrap calibration; `"both"` reports both.
#' @param B Number of bootstrap calibration replicates.
#' @param seed Optional random seed.
#' @param control Controls for constrained null optimization and diagnostics.
#' @param refit_control Controls overriding bootstrap full-model refits.
#' @param verbose Whether to message bootstrap progress.
#' @return A `vcmoe_glrt` object.
#' @export
vcmoe_glrt <- function(fit, data, test = c("coefficient", "constant_all"),
                       coefficient_set = c("expert", "gating", "sigma", "theta"),
                       component = NULL, term = NULL,
                       calibration = c("analytic_epanechnikov", "bootstrap", "both", "none", "parametric_bootstrap"),
                       B = 200L, seed = NULL, control = list(),
                       refit_control = list(), verbose = FALSE) {
  if (!inherits(fit, "vcmoe")) {
    stop("`fit` must be a VCMoE fit.", call. = FALSE)
  }
  test <- match.arg(test)
  calibration <- .vcmoe_glrt_normalize_calibration(calibration)
  glrt_control <- .vcmoe_glrt_control(control)
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 1L || abs(B - round(B)) > 1e-8) {
    stop("`B` must be a whole number of at least 1.", call. = FALSE)
  }
  B <- as.integer(B)
  if (calibration %in% c("analytic_epanechnikov", "both")) {
    .vcmoe_glrt_require_analytic_fit(fit)
  }

  constraint <- .vcmoe_glrt_constraint(fit, test, coefficient_set[[1L]], component, term)
  null_fit <- .vcmoe_glrt_null_fit(fit, constraint, control = control)
  observed <- .vcmoe_glrt_statistic(fit, null_fit, glrt_control$lambda_tol)

  warnings <- .vcmoe_glrt_fit_warnings(fit, glrt_control)
  if (any(c(.vcmoe_default_control(fit$control %||% list())$ridge,
            .vcmoe_default_control(fit$control %||% list())$binomial_ridge,
            .vcmoe_default_control(fit$control %||% list())$negbin_ridge,
            .vcmoe_default_control(fit$control %||% list())$negbin_theta_ridge) > 0)) {
    warnings <- c(warnings, "ridge_used_for_fitting_excluded_from_glrt_statistic")
  }
  null_convergence <- null_fit$diagnostics$glrt_null_optim_convergence
  if (!is.finite(null_convergence) || null_convergence != 0L) {
    warnings <- c(warnings, "null_optimizer_not_fully_converged")
  }

  block_reasons <- character(0L)
  if (!is.na(observed$block_reason)) {
    block_reasons <- c(block_reasons, observed$block_reason)
  }
  if (isTRUE(glrt_control$strict)) {
    block_reasons <- c(block_reasons, setdiff(warnings, "ridge_used_for_fitting_excluded_from_glrt_statistic"))
  }
  status <- if (length(block_reasons)) "blocked" else "ok"
  block_reason <- if (length(block_reasons)) paste(unique(block_reasons), collapse = ";") else NA_character_

  analytic <- list(
    analytic_statistic = NA_real_,
    analytic_p_value = NA_real_,
    rK = NA_real_,
    delta = NA_real_,
    p_test = constraint$p_test,
    C_const = glrt_control$C_const,
    U_length = fit$u_scaling$domain_length %||% 1,
    int_K2 = NA_real_,
    K0 = NA_real_
  )
  if (calibration %in% c("analytic_epanechnikov", "both") && identical(status, "ok")) {
    analytic <- .vcmoe_glrt_analytic_calibration(fit, constraint, observed, glrt_control)
  }

  bootstrap_p_value <- NA_real_
  replicate_summary <- data.frame()
  if (calibration %in% c("bootstrap", "both") && identical(status, "ok")) {
    boot <- .vcmoe_glrt_bootstrap(
      fit, data, null_fit, constraint, B, seed,
      control = control,
      refit_control = refit_control,
      verbose = verbose
    )
    replicate_summary <- boot$replicate_summary
    if (length(boot$successful_lambdas)) {
      bootstrap_p_value <- (1 + sum(boot$successful_lambdas >= observed$lambda_positive)) /
        (1 + length(boot$successful_lambdas))
    }
  }

  p_value <- if (identical(calibration, "analytic_epanechnikov")) {
    analytic$analytic_p_value
  } else if (identical(calibration, "bootstrap")) {
    bootstrap_p_value
  } else if (identical(calibration, "both")) {
    analytic$analytic_p_value
  } else {
    NA_real_
  }

  fit_defaults <- .vcmoe_default_control(fit$control %||% list())
  out <- list(
    fit = fit,
    null_fit = null_fit,
    status = status,
    block_reason = block_reason,
    warnings = unique(warnings),
    statistic = analytic$analytic_statistic,
    lambda = observed$lambda,
    lrt_statistic = observed$lrt_statistic,
    analytic_statistic = analytic$analytic_statistic,
    full_loglik = observed$full_loglik,
    null_loglik = observed$null_loglik,
    p_value = p_value,
    analytic_p_value = analytic$analytic_p_value,
    bootstrap_p_value = bootstrap_p_value,
    rK = analytic$rK,
    delta = analytic$delta,
    p_test = analytic$p_test,
    C_const = analytic$C_const,
    U_length = analytic$U_length,
    int_K2 = analytic$int_K2,
    K0 = analytic$K0,
    null_convergence = null_convergence,
    replicate_summary = replicate_summary,
    settings = list(
      test = test,
      coefficient_set = constraint$coefficient_set,
      component = constraint$component,
      term = constraint$term,
      constraint = constraint$label,
      calibration = calibration,
      B = B,
      seed = seed,
      family = fit$family,
      k = fit$k,
      bandwidth = fit$bandwidth,
      parameterization = fit$parameterization_id %||% vcmoe_parameterization(fit)$id,
      u_scale = fit$u_scale %||% fit$u_scaling$method %||% "unit",
      ridge = fit_defaults$ridge,
      binomial_ridge = fit_defaults$binomial_ridge,
      negbin_ridge = fit_defaults$negbin_ridge,
      negbin_theta_ridge = fit_defaults$negbin_theta_ridge,
      negbin_theta_target = fit_defaults$negbin_theta_target,
      statistic_note = "Analytic calibration uses rK * (full_loglik - null_loglik); ridge penalties are excluded."
    )
  )
  class(out) <- "vcmoe_glrt"
  out
}

#' @export
print.vcmoe_glrt <- function(x, ...) {
  cat("VCMoE GLRT\n")
  cat("  family: ", x$settings$family, "\n", sep = "")
  cat("  test: ", x$settings$test, "\n", sep = "")
  if (!identical(x$settings$test, "constant_all")) {
    cat("  coefficient: ", x$settings$coefficient_set, " / ",
        x$settings$component, " / ", x$settings$term, "\n", sep = "")
  }
  cat("  status: ", x$status, "\n", sep = "")
  if (!identical(x$status, "ok")) {
    cat("  block reason: ", x$block_reason, "\n", sep = "")
  }
  cat("  lambda: ", signif(x$lambda, 5), "\n", sep = "")
  cat("  analytic statistic: ", signif(x$analytic_statistic, 5), "\n", sep = "")
  cat("  calibration: ", x$settings$calibration, "\n", sep = "")
  if (is.finite(x$analytic_p_value)) {
    cat("  analytic p-value: ", signif(x$analytic_p_value, 4), "\n", sep = "")
  }
  if (is.finite(x$bootstrap_p_value)) {
    cat("  bootstrap p-value: ", signif(x$bootstrap_p_value, 4), "\n", sep = "")
  }
  invisible(x)
}
