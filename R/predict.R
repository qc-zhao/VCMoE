.coef_at_u <- function(array3, u_grid, u_values) {
  n <- length(u_values)
  k <- dim(array3)[2L]
  p <- dim(array3)[3L]
  out <- array(NA_real_, dim = c(n, k, p),
               dimnames = list(NULL, dimnames(array3)[[2L]], dimnames(array3)[[3L]]))
  if (length(u_grid) == 1L) {
    for (component in seq_len(k)) {
      for (term in seq_len(p)) {
        out[, component, term] <- array3[1L, component, term]
      }
    }
    return(out)
  }
  for (component in seq_len(k)) {
    for (term in seq_len(p)) {
      out[, component, term] <- stats::approx(
        x = u_grid,
        y = array3[, component, term],
        xout = u_values,
        rule = 2
      )$y
    }
  }
  out
}

.sigma_at_u <- function(sigma, u_grid, u_values) {
  if (is.null(sigma)) {
    return(NULL)
  }
  n <- length(u_values)
  k <- ncol(sigma)
  out <- matrix(NA_real_, nrow = n, ncol = k,
                dimnames = list(NULL, colnames(sigma)))
  if (length(u_grid) == 1L) {
    for (component in seq_len(k)) {
      out[, component] <- sigma[1L, component]
    }
    return(out)
  }
  for (component in seq_len(k)) {
    out[, component] <- exp(stats::approx(
      x = u_grid,
      y = log(pmax(sigma[, component], 1e-12)),
      xout = u_values,
      rule = 2
    )$y)
  }
  out
}

.theta_at_u <- function(theta, u_grid, u_values) {
  .sigma_at_u(theta, u_grid, u_values)
}

.prediction_response <- function(object, newdata) {
  if (is.null(newdata)) {
    return(object$fitted$response %||% list(y = object$fitted$y))
  }

  if (identical(object$family, "gaussian")) {
    if (!is.null(object$response) && object$response %in% names(newdata)) {
      return(.parse_vcmoe_response("gaussian", newdata[[object$response]], object$response))
    }
    return(NULL)
  }
  if (identical(object$family, "negative-binomial")) {
    if (!is.null(object$response) && object$response %in% names(newdata)) {
      return(.parse_vcmoe_response("negative-binomial", newdata[[object$response]], object$response))
    }
    return(NULL)
  }

  info <- object$response_info
  if (is.null(info)) {
    return(NULL)
  }
  if (identical(info$type, "bernoulli")) {
    if (!is.null(info$response) && info$response %in% names(newdata)) {
      return(.parse_vcmoe_response("binomial", newdata[[info$response]], info$response))
    }
    return(NULL)
  }
  if (identical(info$type, "grouped")) {
    if (!is.null(info$success) && !is.null(info$failure) &&
        info$success %in% names(newdata) && info$failure %in% names(newdata)) {
      response_matrix <- cbind(newdata[[info$success]], newdata[[info$failure]])
      colnames(response_matrix) <- c(info$success, info$failure)
      return(.parse_vcmoe_response("binomial", response_matrix, info$response))
    }
  }
  NULL
}

.prediction_design <- function(object, newdata, u) {
  if (is.null(newdata)) {
    if (is.null(object$fitted)) {
      stop("The fit does not contain training data; provide `newdata` and `u`.", call. = FALSE)
    }
    return(list(
      y = object$fitted$y,
      response = .prediction_response(object, NULL),
      z_design = object$fitted$z_design,
      x_design = object$fitted$x_design,
      expert_offset = object$fitted$expert_offset %||% rep(0, length(object$fitted$u)),
      u = object$fitted$u
    ))
  }

  expert_terms <- stats::delete.response(object$terms$expert)
  expert_frame <- stats::model.frame(expert_terms, data = newdata, na.action = stats::na.pass)
  gating_frame <- stats::model.frame(object$terms$gating, data = newdata, na.action = stats::na.pass)
  z_design <- stats::model.matrix(expert_terms, data = expert_frame)
  x_design <- stats::model.matrix(object$terms$gating, data = gating_frame)
  expert_offset <- .extract_model_offset(expert_frame)
  gating_offset <- .extract_model_offset(gating_frame)
  if (any(abs(gating_offset) > 0)) {
    stop("Gating-side `offset()` terms are not supported in v0.", call. = FALSE)
  }

  if (is.null(u)) {
    if (is.null(object$u_name) || !object$u_name %in% names(newdata)) {
      stop("Provide `u` for `newdata`, or include the original `u` column.", call. = FALSE)
    }
    u_values <- newdata[[object$u_name]]
  } else if (is.character(u) && length(u) == 1L) {
    u_values <- newdata[[u]]
  } else {
    u_values <- u
  }

  response <- .prediction_response(object, newdata)
  y <- if (is.null(response)) NULL else response$observed
  list(
    y = y,
    response = response,
    z_design = z_design,
    x_design = x_design,
    expert_offset = expert_offset,
    u = .vcmoe_scale_u(as.numeric(u_values), object$u_scaling)
  )
}

.predict_components <- function(object, newdata = NULL, u = NULL) {
  design <- .prediction_design(object, newdata, u)
  expert <- .coef_at_u(object$coefficients$expert, object$u_grid, design$u)
  gating <- .coef_at_u(object$coefficients$gating, object$u_grid, design$u)
  sigma <- .sigma_at_u(object$coefficients$sigma, object$u_grid, design$u)
  theta <- .theta_at_u(object$coefficients$theta, object$u_grid, design$u)
  n <- length(design$u)
  k <- object$k

  component_mean <- matrix(NA_real_, nrow = n, ncol = k,
                           dimnames = list(NULL, paste0("component", seq_len(k))))
  logits <- component_mean
  for (component in seq_len(k)) {
    eta <- rowSums(design$z_design * expert[, component, ]) + design$expert_offset
    if (identical(object$family, "binomial")) {
      component_mean[, component] <- stats::plogis(eta)
    } else if (identical(object$family, "negative-binomial")) {
      component_mean[, component] <- exp(pmin(pmax(eta, -20), 20))
    } else {
      component_mean[, component] <- eta
    }
    logits[, component] <- rowSums(design$x_design * gating[, component, ])
  }
  prior <- .row_softmax(logits)

  posterior <- prior
  posterior_source <- "gating"
  if (!is.null(design$response)) {
    log_terms <- matrix(0, nrow = n, ncol = k)
    for (component in seq_len(k)) {
      if (identical(object$family, "binomial")) {
        expert_loglik <- stats::dbinom(
          design$response$success,
          size = design$response$trials,
          prob = component_mean[, component],
          log = TRUE
        )
      } else if (identical(object$family, "negative-binomial")) {
        expert_loglik <- stats::dnbinom(
          design$response$y,
          size = theta[, component],
          mu = component_mean[, component],
          log = TRUE
        )
      } else {
        expert_loglik <- stats::dnorm(
          design$response$y,
          component_mean[, component],
          sigma[, component],
          log = TRUE
        )
      }
      log_terms[, component] <- log(pmax(prior[, component], 1e-12)) + expert_loglik
    }
    posterior <- exp(log_terms - .row_log_sum_exp(log_terms))
    posterior <- posterior / rowSums(posterior)
    posterior_source <- "posterior"
  }

  list(
    mean = component_mean,
    sigma = sigma,
    theta = theta,
    prior = prior,
    posterior = posterior,
    posterior_source = posterior_source
  )
}

#' Extract VCMoE coefficients
#'
#' @param object A `vcmoe` object.
#' @param type Coefficient block to return.
#' @param ... Unused.
#' @return A list or array of fitted coefficient functions.
#' @export
coef.vcmoe <- function(object, type = c("all", "expert", "gating", "sigma", "sigma_slope", "theta"), ...) {
  type <- match.arg(type)
  if (type == "all") {
    return(object$coefficients)
  }
  object$coefficients[[type]]
}

#' Predict from a VCMoE fit
#'
#' @param object A `vcmoe` object.
#' @param newdata Optional data frame.
#' @param u Optional index values for `newdata`.
#' @param type Prediction type: fitted mean, posterior probabilities,
#'   component-specific means, or gating probabilities.
#' @param ... Unused.
#' @return A vector or matrix depending on `type`.
#' @export
predict.vcmoe <- function(object, newdata = NULL, u = NULL,
                          type = c("mean", "posterior", "component", "prior"), ...) {
  type <- match.arg(type)
  pieces <- .predict_components(object, newdata = newdata, u = u)
  if (type == "component") {
    return(pieces$mean)
  }
  if (type == "posterior") {
    out <- pieces$posterior
    attr(out, "source") <- pieces$posterior_source
    return(out)
  }
  if (type == "prior") {
    return(pieces$prior)
  }
  if (identical(object$family, "binomial") || identical(object$family, "negative-binomial")) {
    return(rowSums(pieces$prior * pieces$mean))
  }
  rowSums(pieces$posterior * pieces$mean)
}
