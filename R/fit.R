.make_expert_design <- function(z_design, du) {
  cbind(z_design, z_design * as.numeric(du))
}

.fit_binomial_logistic <- function(success, trials, design, weights, start, control) {
  weights <- pmax(weights, 0)
  if (sum(weights) <= 1e-10) {
    return(list(
      coef = start,
      convergence = NA_integer_,
      gradient_norm = NA_real_,
      value = NA_real_
    ))
  }
  start <- as.numeric(start)
  ridge <- control$binomial_ridge %||% control$ridge
  objective <- function(par) {
    eta <- as.numeric(design %*% par)
    prob <- stats::plogis(eta)
    prob <- pmin(pmax(prob, 1e-10), 1 - 1e-10)
    value <- -sum(weights * (
      success * log(prob) + (trials - success) * log1p(-prob)
    ))
    value + 0.5 * ridge * sum(par^2)
  }
  gradient <- function(par) {
    eta <- as.numeric(design %*% par)
    prob <- stats::plogis(eta)
    as.numeric(-crossprod(design, weights * (success - trials * prob)) + ridge * par)
  }
  opt <- tryCatch(
    stats::optim(start, objective, gradient, method = "BFGS",
                 control = list(maxit = control$gating_maxit)),
    error = function(e) NULL
  )
  if (is.null(opt) || any(!is.finite(opt$par))) {
    return(list(
      coef = start,
      convergence = NA_integer_,
      gradient_norm = NA_real_,
      value = NA_real_
    ))
  }
  coef <- as.numeric(opt$par)
  list(
    coef = coef,
    convergence = opt$convergence,
    gradient_norm = sqrt(sum(gradient(coef)^2)),
    value = opt$value
  )
}

.bounded_exp <- function(eta, control) {
  exp(pmin(pmax(eta, control$negbin_eta_min), control$negbin_eta_max))
}

.fit_negbin_loglink <- function(y, design, offset, weights, start, theta, control) {
  weights <- pmax(weights, 0)
  if (sum(weights) <= control$min_component_weight) {
    return(list(
      coef = start,
      convergence = NA_integer_,
      gradient_norm = NA_real_,
      value = NA_real_
    ))
  }
  start <- as.numeric(start)
  theta <- pmin(pmax(theta, control$negbin_theta_min), control$negbin_theta_max)
  ridge <- control$negbin_ridge %||% control$ridge
  objective <- function(par) {
    eta <- offset + as.numeric(design %*% par)
    mu <- .bounded_exp(eta, control)
    value <- -sum(weights * stats::dnbinom(y, size = theta, mu = mu, log = TRUE))
    value + 0.5 * ridge * sum(par^2)
  }
  gradient <- function(par) {
    eta_raw <- offset + as.numeric(design %*% par)
    eta <- pmin(pmax(eta_raw, control$negbin_eta_min), control$negbin_eta_max)
    mu <- exp(eta)
    score_eta <- y - mu * (y + theta) / (mu + theta)
    score_eta <- score_eta * as.numeric(eta_raw > control$negbin_eta_min & eta_raw < control$negbin_eta_max)
    as.numeric(-crossprod(design, weights * score_eta) + ridge * par)
  }
  opt <- tryCatch(
    stats::optim(start, objective, gradient, method = "BFGS",
                 control = list(maxit = control$negbin_mstep_maxit)),
    error = function(e) NULL
  )
  if (is.null(opt) || any(!is.finite(opt$par))) {
    return(list(
      coef = start,
      convergence = NA_integer_,
      gradient_norm = NA_real_,
      value = NA_real_
    ))
  }
  coef <- as.numeric(opt$par)
  list(
    coef = coef,
    convergence = opt$convergence,
    gradient_norm = sqrt(sum(gradient(coef)^2)),
    value = opt$value
  )
}

.fit_negbin_theta <- function(y, mu, weights, previous_theta, control,
                              penalty_mass = sum(pmax(weights, 0))) {
  weights <- pmax(weights, 0)
  if (sum(weights) <= control$min_component_weight) {
    return(pmin(pmax(previous_theta, control$negbin_theta_min), control$negbin_theta_max))
  }
  lower <- log(control$negbin_theta_min)
  upper <- log(control$negbin_theta_max)
  theta_target <- control$negbin_theta_target %||% previous_theta
  theta_target <- as.numeric(theta_target[[1L]])
  if (!is.finite(theta_target) || theta_target <= 0) {
    theta_target <- previous_theta
  }
  theta_target <- pmin(pmax(theta_target, control$negbin_theta_min), control$negbin_theta_max)
  log_theta_target <- log(theta_target)
  theta_ridge <- as.numeric(control$negbin_theta_ridge %||% 0)
  if (!is.finite(theta_ridge) || theta_ridge < 0) {
    theta_ridge <- 0
  }
  penalty_scale <- penalty_mass
  objective <- function(log_theta) {
    theta <- exp(log_theta)
    value <- -sum(weights * stats::dnbinom(y, size = theta, mu = mu, log = TRUE))
    value + 0.5 * theta_ridge * penalty_scale * (log_theta - log_theta_target)^2
  }
  opt <- tryCatch(
    stats::optimize(objective, interval = c(lower, upper)),
    error = function(e) NULL
  )
  if (is.null(opt) || !is.finite(opt$minimum)) {
    return(pmin(pmax(previous_theta, control$negbin_theta_min), control$negbin_theta_max))
  }
  exp(opt$minimum)
}

.gaussian_sigma_from_eta <- function(log_sigma_eta, control) {
  exp(pmin(pmax(log_sigma_eta, log(control$min_sigma)), log(control$max_sigma)))
}

.gaussian_component_sigma <- function(params, component, du, control) {
  slope <- if (!is.null(params$sigma_slope)) params$sigma_slope[[component]] else 0
  .gaussian_sigma_from_eta(log(pmax(params$sigma[[component]], control$min_sigma)) + slope * du, control)
}

.fit_gaussian_locallinear <- function(y, design, du, weights, start_coef,
                                      start_sigma, start_sigma_slope, control) {
  weights <- pmax(weights, 0)
  if (sum(weights) <= control$min_component_weight) {
    return(list(
      coef = start_coef,
      sigma = pmin(pmax(start_sigma, control$min_sigma), control$max_sigma),
      sigma_slope = start_sigma_slope,
      convergence = NA_integer_,
      gradient_norm = NA_real_,
      value = NA_real_
    ))
  }
  start_coef <- as.numeric(start_coef)
  start_sigma <- pmin(pmax(start_sigma, control$min_sigma), control$max_sigma)
  start_sigma_slope <- as.numeric(start_sigma_slope %||% 0)
  if (!is.finite(start_sigma_slope)) {
    start_sigma_slope <- 0
  }
  start <- c(start_coef, log(start_sigma), start_sigma_slope)
  coef_len <- length(start_coef)
  ridge <- control$ridge
  sigma_ridge <- control$gaussian_sigma_ridge %||% control$ridge
  objective <- function(par) {
    if (any(!is.finite(par))) {
      return(1e100)
    }
    beta <- par[seq_len(coef_len)]
    log_sigma_eta <- par[[coef_len + 1L]] + par[[coef_len + 2L]] * du
    sigma <- .gaussian_sigma_from_eta(log_sigma_eta, control)
    mu <- as.numeric(design %*% beta)
    value <- -sum(weights * stats::dnorm(y, mean = mu, sd = sigma, log = TRUE))
    value + 0.5 * ridge * sum(beta^2) + 0.5 * sigma_ridge * par[[coef_len + 2L]]^2
  }
  opt <- tryCatch(
    stats::optim(
      start,
      objective,
      method = "BFGS",
      control = list(maxit = control$gaussian_mstep_maxit)
    ),
    error = function(e) NULL
  )
  if (is.null(opt) || any(!is.finite(opt$par))) {
    return(list(
      coef = start_coef,
      sigma = start_sigma,
      sigma_slope = start_sigma_slope,
      convergence = NA_integer_,
      gradient_norm = NA_real_,
      value = NA_real_
    ))
  }
  coef <- as.numeric(opt$par[seq_len(coef_len)])
  sigma <- .gaussian_sigma_from_eta(opt$par[[coef_len + 1L]], control)
  sigma_slope <- as.numeric(opt$par[[coef_len + 2L]])
  list(
    coef = coef,
    sigma = sigma,
    sigma_slope = sigma_slope,
    convergence = opt$convergence,
    gradient_norm = NA_real_,
    value = opt$value
  )
}

.mstep_gating <- function(gating_design, posterior, weights, beta_full, control) {
  k <- ncol(posterior)
  n_terms <- ncol(gating_design)
  weights <- pmax(weights, 0)
  ridge <- control$ridge

  if (k == 2L) {
    start <- beta_full[1L, ] - beta_full[2L, ]
    objective <- function(par) {
      eta <- as.numeric(gating_design %*% par)
      p1 <- stats::plogis(eta)
      p1 <- pmin(pmax(p1, 1e-10), 1 - 1e-10)
      value <- -sum(weights * (
        posterior[, 1L] * log(p1) + posterior[, 2L] * log1p(-p1)
      ))
      value + 0.5 * ridge * sum(par^2)
    }
    gradient <- function(par) {
      eta <- as.numeric(gating_design %*% par)
      p1 <- stats::plogis(eta)
      as.numeric(-crossprod(gating_design, weights * (posterior[, 1L] - p1)) + ridge * par)
    }
    opt <- tryCatch(
      stats::optim(start, objective, gradient, method = "BFGS",
                   control = list(maxit = control$gating_maxit)),
      error = function(e) NULL
    )
    if (is.null(opt) || any(!is.finite(opt$par))) {
      return(beta_full)
    }
    beta <- opt$par
    return(rbind(beta / 2, -beta / 2))
  }

  start_matrix <- beta_full[-1L, , drop = FALSE] -
    matrix(beta_full[1L, ], nrow = k - 1L, ncol = n_terms, byrow = TRUE)
  start <- as.vector(t(start_matrix))

  objective <- function(par) {
    theta <- matrix(par, nrow = k - 1L, byrow = TRUE)
    eta <- cbind(0, gating_design %*% t(theta))
    prob <- .row_softmax(eta)
    -sum(weights * rowSums(posterior * log(pmax(prob, 1e-12)))) +
      0.5 * ridge * sum(par^2)
  }
  gradient <- function(par) {
    theta <- matrix(par, nrow = k - 1L, byrow = TRUE)
    eta <- cbind(0, gating_design %*% t(theta))
    prob <- .row_softmax(eta)
    grad <- matrix(0, nrow = k - 1L, ncol = n_terms)
    for (component in seq_len(k - 1L)) {
      grad[component, ] <- -crossprod(
        gating_design,
        weights * (posterior[, component + 1L] - prob[, component + 1L])
      )
    }
    as.vector(t(grad)) + ridge * par
  }
  opt <- tryCatch(
    stats::optim(start, objective, gradient, method = "BFGS",
                 control = list(maxit = control$gating_maxit)),
    error = function(e) NULL
  )
  if (is.null(opt) || any(!is.finite(opt$par))) {
    return(beta_full)
  }

  theta <- matrix(opt$par, nrow = k - 1L, byrow = TRUE)
  beta_baseline <- rbind(rep(0, n_terms), theta)
  .center_logits(beta_baseline)
}

.m_step <- function(y, trials, family, z_design, x_design, du, offset, weights, posterior, previous, control) {
  k <- ncol(posterior)
  q <- ncol(z_design)
  p <- ncol(x_design)
  expert_design <- .make_expert_design(z_design, du)
  gating_design <- .make_expert_design(x_design, du)

  expert_coef <- matrix(0, nrow = k, ncol = 2L * q)
  sigma <- if (identical(family, "gaussian")) numeric(k) else NULL
  sigma_slope <- if (identical(family, "gaussian")) numeric(k) else NULL
  theta <- if (identical(family, "negative-binomial")) numeric(k) else NULL
  expert_convergence <- rep(NA_integer_, k)
  expert_gradient_norm <- rep(NA_real_, k)
  expert_coef_norm <- rep(NA_real_, k)
  for (component in seq_len(k)) {
    wc <- weights * posterior[, component]
    if (identical(family, "gaussian")) {
      gaussian_fit <- .fit_gaussian_locallinear(
        y,
        expert_design,
        du,
        wc,
        previous$expert_coef[component, ],
        previous$sigma[component],
        previous$sigma_slope[component] %||% 0,
        control
      )
      coef_c <- gaussian_fit$coef
      sigma[component] <- gaussian_fit$sigma
      sigma_slope[component] <- gaussian_fit$sigma_slope
      expert_convergence[component] <- gaussian_fit$convergence
      expert_gradient_norm[component] <- gaussian_fit$gradient_norm
    } else if (identical(family, "binomial")) {
      expert_fit <- .fit_binomial_logistic(
        y,
        trials,
        expert_design,
        wc,
        previous$expert_coef[component, ],
        control
      )
      coef_c <- expert_fit$coef
      expert_convergence[component] <- expert_fit$convergence
      expert_gradient_norm[component] <- expert_fit$gradient_norm
    } else {
      expert_fit <- .fit_negbin_loglink(
        y,
        expert_design,
        offset,
        wc,
        previous$expert_coef[component, ],
        previous$theta[component],
        control
      )
      coef_c <- expert_fit$coef
      eta_c <- offset + as.numeric(expert_design %*% coef_c)
      mu_c <- .bounded_exp(eta_c, control)
      theta[component] <- .fit_negbin_theta(
        y,
        mu_c,
        wc,
        previous$theta[component],
        control,
        penalty_mass = sum(if (identical(.vcmoe_theta_penalty_mass(control), "kernel")) weights else wc)
      )
      expert_convergence[component] <- expert_fit$convergence
      expert_gradient_norm[component] <- expert_fit$gradient_norm
    }
    expert_coef[component, ] <- coef_c
    expert_coef_norm[component] <- sqrt(sum(coef_c^2))
  }

  beta_full <- .mstep_gating(gating_design, posterior, weights, previous$gating_coef, control)
  list(
    expert_coef = expert_coef,
    gating_coef = beta_full,
    sigma = sigma,
    sigma_slope = sigma_slope,
    theta = theta,
    expert_convergence = expert_convergence,
    expert_gradient_norm = expert_gradient_norm,
    expert_coef_norm = expert_coef_norm
  )
}

# Known-label fits: `control$fixed_membership` (an n x k 0/1 matrix on the
# fitted rows) drops every other component from each row's mixture, so the
# E-step returns the labels and each local fit is the kernel-weighted
# complete-data fit (component-wise experts, multinomial-logit gate).
.vcmoe_mask_log_terms <- function(log_terms, membership) {
  if (!is.null(membership)) {
    log_terms[membership == 0] <- -Inf
  }
  log_terms
}

# Accepts component indices or a factor (levels in component order) of length
# n, or an n x k 0/1 matrix with one 1 per row.
.vcmoe_fixed_membership <- function(value, n, k) {
  if (is.data.frame(value)) {
    value <- as.matrix(value)
  }
  if (is.matrix(value)) {
    if (!identical(dim(value), c(as.integer(n), as.integer(k))) || anyNA(value) ||
        !all(value %in% c(0, 1)) || any(rowSums(value) != 1)) {
      stop("A `fixed_membership` matrix must be nrow(data) x k with one 1 per row.", call. = FALSE)
    }
    out <- matrix(as.numeric(value), nrow = n, ncol = k)
  } else {
    index <- if (is.factor(value)) as.integer(value) else suppressWarnings(as.integer(value))
    if (length(index) != n || anyNA(index) || any(index < 1L | index > k) ||
        (!is.factor(value) && any(index != value))) {
      stop("`fixed_membership` must give a component in 1..k for every row of `data`.", call. = FALSE)
    }
    out <- matrix(0, nrow = n, ncol = k)
    out[cbind(seq_len(n), index)] <- 1
  }
  if (any(colSums(out) == 0)) {
    stop("Every component must have at least one row in `fixed_membership`.", call. = FALSE)
  }
  out
}

.local_expert_penalty <- function(family, expert_coef, sigma_slope, theta, weights, posterior, control) {
  if (identical(family, "gaussian")) {
    sigma_ridge <- control$gaussian_sigma_ridge %||% control$ridge
    return(0.5 * control$ridge * sum(expert_coef^2) + 0.5 * sigma_ridge * sum(sigma_slope^2))
  }
  if (identical(family, "binomial")) {
    return(0.5 * (control$binomial_ridge %||% control$ridge) * sum(expert_coef^2))
  }
  out <- 0.5 * (control$negbin_ridge %||% control$ridge) * sum(expert_coef^2)
  if (.vcmoe_theta_ridge_active(control)) {
    target <- as.numeric(control$negbin_theta_target[1L])
    target <- min(max(target, control$negbin_theta_min), control$negbin_theta_max)
    mass <- if (identical(.vcmoe_theta_penalty_mass(control), "kernel")) {
      rep(sum(pmax(weights, 0)), length(theta))
    } else {
      colSums(pmax(weights, 0) * posterior)
    }
    out <- out + 0.5 * as.numeric(control$negbin_theta_ridge) * sum(mass * (log(theta) - log(target))^2)
  }
  out
}

.local_gating_penalty <- function(gating_coef, control) {
  contrast <- gating_coef[-1L, , drop = FALSE] -
    matrix(gating_coef[1L, ], nrow = nrow(gating_coef) - 1L, ncol = ncol(gating_coef), byrow = TRUE)
  0.5 * control$ridge * sum(contrast^2)
}

# Ridge terms of the local objective that EM ascends. Each M-step maximises the
# kernel-weighted complete-data log-likelihood minus these penalties, so EM is
# monotone in `loglik - penalty`, not in the log-likelihood; the stop rule and
# start selection use the penalized value. The gating ridge acts on contrasts
# with component 1 (the M-step parameterization). A theta ridge without a fixed
# target is centred on the previous theta (proximal), vanishes at a fixed point
# and is not part of the objective.
.local_ridge_penalty <- function(family, params, weights, posterior, control) {
  .local_gating_penalty(params$gating_coef, control) +
    .local_expert_penalty(family, params$expert_coef, params$sigma_slope, params$theta,
                          weights, posterior, control)
}

.local_loglik <- function(y, trials, family, z_design, x_design, du, offset, weights, params, control) {
  expert_design <- .make_expert_design(z_design, du)
  gating_design <- .make_expert_design(x_design, du)
  pi <- .gating_prob(gating_design, params$gating_coef)
  k <- nrow(params$expert_coef)
  log_terms <- matrix(0, nrow = length(y), ncol = k)
  for (component in seq_len(k)) {
    eta <- offset + as.numeric(expert_design %*% params$expert_coef[component, ])
    if (identical(family, "gaussian")) {
      expert_loglik <- stats::dnorm(
        y,
        mean = eta,
        sd = .gaussian_component_sigma(params, component, du, control),
        log = TRUE
      )
    } else if (identical(family, "binomial")) {
      expert_loglik <- stats::dbinom(y, size = trials, prob = stats::plogis(eta), log = TRUE)
    } else {
      expert_loglik <- stats::dnbinom(
        y,
        size = params$theta[component],
        mu = .bounded_exp(eta, control),
        log = TRUE
      )
    }
    log_terms[, component] <- log(pmax(pi[, component], 1e-12)) + expert_loglik
  }
  log_terms <- .vcmoe_mask_log_terms(log_terms, control$fixed_membership)
  sum(weights * .row_log_sum_exp(log_terms))
}

.e_step <- function(y, trials, family, z_design, x_design, du, offset, params, control) {
  expert_design <- .make_expert_design(z_design, du)
  gating_design <- .make_expert_design(x_design, du)
  pi <- .gating_prob(gating_design, params$gating_coef)
  k <- nrow(params$expert_coef)
  log_terms <- matrix(0, nrow = length(y), ncol = k)
  for (component in seq_len(k)) {
    eta <- offset + as.numeric(expert_design %*% params$expert_coef[component, ])
    if (identical(family, "gaussian")) {
      expert_loglik <- stats::dnorm(
        y,
        mean = eta,
        sd = .gaussian_component_sigma(params, component, du, control),
        log = TRUE
      )
    } else if (identical(family, "binomial")) {
      expert_loglik <- stats::dbinom(y, size = trials, prob = stats::plogis(eta), log = TRUE)
    } else {
      expert_loglik <- stats::dnbinom(
        y,
        size = params$theta[component],
        mu = .bounded_exp(eta, control),
        log = TRUE
      )
    }
    log_terms[, component] <- log(pmax(pi[, component], 1e-12)) + expert_loglik
  }
  log_terms <- .vcmoe_mask_log_terms(log_terms, control$fixed_membership)
  log_norm <- .row_log_sum_exp(log_terms)
  posterior <- exp(log_terms - log_norm)
  posterior / rowSums(posterior)
}

.params_to_local <- function(params, posterior, z_design, x_design, du, loglik, iterations, converged) {
  q <- ncol(z_design)
  p <- ncol(x_design)
  gating_design <- .make_expert_design(x_design, du)
  list(
    expert_intercept = params$expert_coef[, seq_len(q), drop = FALSE],
    expert_slope = params$expert_coef[, q + seq_len(q), drop = FALSE],
    gating_intercept = params$gating_coef[, seq_len(p), drop = FALSE],
    gating_slope = params$gating_coef[, p + seq_len(p), drop = FALSE],
    sigma = params$sigma,
    sigma_slope = params$sigma_slope,
    theta = params$theta,
    posterior = posterior,
    pi = .gating_prob(gating_design, params$gating_coef),
    loglik = loglik,
    iterations = iterations,
    converged = converged,
    expert_convergence = params$expert_convergence,
    expert_gradient_norm = params$expert_gradient_norm,
    expert_coef_norm = params$expert_coef_norm
  )
}

.local_to_params <- function(local) {
  list(
    expert_coef = cbind(local$expert_intercept, local$expert_slope),
    gating_coef = cbind(local$gating_intercept, local$gating_slope),
    sigma = local$sigma,
    sigma_slope = local$sigma_slope,
    theta = local$theta
  )
}

.initial_tau <- function(initial_y, k, random = FALSE) {
  n <- length(initial_y)
  if (random) {
    raw <- matrix(stats::runif(n * k), nrow = n, ncol = k)
    return(raw / rowSums(raw))
  }

  groups <- tryCatch(stats::kmeans(initial_y, centers = k, nstart = 10)$cluster,
                     error = function(e) NULL)
  if (is.null(groups) || any(is.na(groups))) {
    groups <- cut(rank(initial_y, ties.method = "first"),
                  breaks = k, labels = FALSE, include.lowest = TRUE)
  }
  tau <- matrix(0.05 / max(k - 1L, 1L), nrow = n, ncol = k)
  for (component in seq_len(k)) {
    tau[groups == component, component] <- 0.95
  }
  tau / rowSums(tau)
}

.initial_params_from_tau <- function(y, trials, family, z_design, x_design, du, offset, weights, posterior, control) {
  q <- ncol(z_design)
  p <- ncol(x_design)
  k <- ncol(posterior)
  previous <- list(
    expert_coef = matrix(0, nrow = k, ncol = 2L * q),
    gating_coef = matrix(0, nrow = k, ncol = 2L * p),
    sigma = if (identical(family, "gaussian")) rep(max(stats::sd(y), control$min_sigma), k) else NULL,
    sigma_slope = if (identical(family, "gaussian")) rep(0, k) else NULL,
    theta = if (identical(family, "negative-binomial")) rep(control$negbin_init_theta, k) else NULL
  )
  .m_step(y, trials, family, z_design, x_design, du, offset, weights, posterior, previous, control)
}

.initial_binomial_structured_params <- function(y, trials, z_design, x_design, du, weights, k, control, start_id) {
  q <- ncol(z_design)
  p <- ncol(x_design)
  expert_design <- .make_expert_design(z_design, du)
  base_fit <- .fit_binomial_logistic(
    y,
    trials,
    expert_design,
    weights,
    rep(0, 2L * q),
    control
  )
  base_coef <- base_fit$coef
  centers <- seq(-(k - 1) / 2, (k - 1) / 2, length.out = k)
  direction <- if (start_id %% 2L == 0L) 1 else -1
  scale <- control$binomial_start_separation * (1 + 0.20 * max(start_id - 1L, 0L))
  expert_coef <- matrix(base_coef, nrow = k, ncol = length(base_coef), byrow = TRUE)
  jitter <- control$binomial_start_jitter
  for (component in seq_len(k)) {
    component_shift <- direction * centers[component] * scale
    expert_coef[component, 1L] <- expert_coef[component, 1L] + component_shift
    if (q >= 2L) {
      expert_coef[component, 2L] <- expert_coef[component, 2L] - 0.50 * component_shift
    }
    if (2L * q >= q + 1L) {
      expert_coef[component, q + 1L] <- expert_coef[component, q + 1L] + 0.20 * component_shift
    }
    if (jitter > 0 && start_id > 1L) {
      expert_coef[component, ] <- expert_coef[component, ] +
        stats::rnorm(ncol(expert_coef), sd = jitter / sqrt(start_id))
    }
  }

  gating_coef <- matrix(0, nrow = k, ncol = 2L * p)
  for (component in seq_len(k)) {
    gating_coef[component, 1L] <- 0.20 * direction * centers[component]
    if (p >= 2L) {
      gating_coef[component, 2L] <- -0.10 * direction * centers[component]
    }
  }
  gating_coef <- .center_logits(gating_coef)

  list(
    expert_coef = expert_coef,
    gating_coef = gating_coef,
    sigma = NULL,
    sigma_slope = NULL,
    expert_convergence = rep(NA_integer_, k),
    expert_gradient_norm = rep(NA_real_, k),
    expert_coef_norm = apply(expert_coef, 1L, function(x) sqrt(sum(x^2)))
  )
}

.initial_negbin_structured_params <- function(y, z_design, x_design, du, offset, weights, k, control, start_id) {
  q <- ncol(z_design)
  p <- ncol(x_design)
  expert_design <- .make_expert_design(z_design, du)
  base_y <- log(pmax(y, 0) + 0.5) - offset
  base_coef <- .fit_wls(base_y, expert_design, weights, control$negbin_ridge %||% control$ridge)
  centers <- seq(-(k - 1) / 2, (k - 1) / 2, length.out = k)
  direction <- if (start_id %% 2L == 0L) 1 else -1
  scale <- control$negbin_start_separation * (1 + 0.15 * max(start_id - 1L, 0L))
  expert_coef <- matrix(base_coef, nrow = k, ncol = length(base_coef), byrow = TRUE)
  jitter <- control$negbin_start_jitter
  for (component in seq_len(k)) {
    component_shift <- direction * centers[component] * scale
    expert_coef[component, 1L] <- expert_coef[component, 1L] + component_shift
    if (q >= 2L) {
      expert_coef[component, 2L] <- expert_coef[component, 2L] - 0.25 * component_shift
    }
    if (jitter > 0 && start_id > 1L) {
      expert_coef[component, ] <- expert_coef[component, ] +
        stats::rnorm(ncol(expert_coef), sd = jitter / sqrt(start_id))
    }
  }

  gating_coef <- matrix(0, nrow = k, ncol = 2L * p)
  for (component in seq_len(k)) {
    gating_coef[component, 1L] <- 0.15 * direction * centers[component]
    if (p >= 2L) {
      gating_coef[component, 2L] <- -0.08 * direction * centers[component]
    }
  }
  gating_coef <- .center_logits(gating_coef)

  list(
    expert_coef = expert_coef,
    gating_coef = gating_coef,
    sigma = NULL,
    sigma_slope = NULL,
    theta = rep(control$negbin_init_theta, k),
    expert_convergence = rep(NA_integer_, k),
    expert_gradient_norm = rep(NA_real_, k),
    expert_coef_norm = apply(expert_coef, 1L, function(x) sqrt(sum(x^2)))
  )
}

.local_em <- function(y, trials, family, initial_y, z_design, x_design, offset, u, u0, bandwidth, k, control,
                      estimation_spec, init_local = NULL, init_tau = NULL, init_params = NULL) {
  weights <- .local_weights(u, u0, bandwidth, estimation_spec)
  du <- .local_du(u, u0, bandwidth, estimation_spec)
  if (!is.null(init_params)) {
    params <- init_params
  } else if (!is.null(init_local)) {
    params <- .local_to_params(init_local)
  } else {
    posterior0 <- init_tau %||% .initial_tau(initial_y, k)
    params <- .initial_params_from_tau(y, trials, family, z_design, x_design, du, offset, weights, posterior0, control)
  }

  loglik <- -Inf
  objective <- -Inf
  objective_delta <- NA_real_
  converged <- FALSE
  posterior <- NULL
  iter <- 0L
  for (iter in seq_len(control$maxit)) {
    posterior <- .e_step(y, trials, family, z_design, x_design, du, offset, params, control)
    new_params <- .m_step(y, trials, family, z_design, x_design, du, offset, weights, posterior, params, control)
    new_loglik <- .local_loglik(y, trials, family, z_design, x_design, du, offset, weights, new_params, control)
    new_objective <- new_loglik - .local_ridge_penalty(family, new_params, weights, posterior, control)
    objective_delta <- if (is.finite(objective)) new_objective - objective else NA_real_
    # Stop on the penalized objective EM ascends; the unpenalized log-likelihood
    # is not monotone and its change can pass through zero while EM still climbs.
    # `legacy_em_rule` (internal) restores the VCMoE 0.2.0 rule for regression tests.
    did_converge <- if (isTRUE(control$legacy_em_rule)) {
      is.finite(loglik) && abs(new_loglik - loglik) <= control$tol * (abs(loglik) + 1)
    } else {
      is.finite(objective) && abs(new_objective - objective) <= control$tol * (abs(objective) + 1)
    }
    params <- new_params
    loglik <- new_loglik
    objective <- new_objective
    if (did_converge) {
      converged <- TRUE
      break
    }
  }

  posterior <- .e_step(y, trials, family, z_design, x_design, du, offset, params, control)
  out <- .params_to_local(params, posterior, z_design, x_design, du, loglik, iter, converged)
  out$objective <- objective
  out$objective_delta <- objective_delta
  out$stop_reason <- if (converged) "tolerance" else if (is.finite(objective)) "maxit" else "nonfinite_objective"
  out
}

.fit_local_best <- function(y, trials, family, initial_y, z_design, x_design, offset, u, u0, bandwidth, k, control,
                            previous, estimation_spec, init_locals = NULL) {
  # `init_locals` (a refit or a grid-parallel second pass) replaces the usual
  # starts with EM runs initialised from the given local fits.
  # Known labels give one complete-data fit: a single start from the labels.
  fixed <- control$fixed_membership
  starts <- if (!is.null(fixed)) {
    1L
  } else if (is.null(init_locals)) {
    max(1L, as.integer(control$n_starts))
  } else {
    length(init_locals)
  }
  candidates <- vector("list", starts)
  loglik <- rep(-Inf, starts)
  objective <- rep(-Inf, starts)
  start_converged <- rep(FALSE, starts)

  for (start_id in seq_len(starts)) {
    init_local <- NULL
    init_tau <- NULL
    init_params <- NULL
    if (!is.null(fixed)) {
      init_tau <- fixed
    } else if (!is.null(init_locals)) {
      init_local <- init_locals[[start_id]]
    } else if (start_id == 1L && !is.null(previous)) {
      init_local <- previous
    } else if (identical(family, "binomial") && isTRUE(control$binomial_structured_starts)) {
      weights <- .local_weights(u, u0, bandwidth, estimation_spec)
      du <- .local_du(u, u0, bandwidth, estimation_spec)
      init_params <- .initial_binomial_structured_params(
        y,
        trials,
        z_design,
        x_design,
        du,
        weights,
        k,
        control,
        start_id
      )
    } else if (identical(family, "negative-binomial") && isTRUE(control$negbin_structured_starts)) {
      weights <- .local_weights(u, u0, bandwidth, estimation_spec)
      du <- .local_du(u, u0, bandwidth, estimation_spec)
      init_params <- .initial_negbin_structured_params(
        y,
        z_design,
        x_design,
        du,
        offset,
        weights,
        k,
        control,
        start_id
      )
    } else {
      init_tau <- .initial_tau(initial_y, k, random = start_id > 1L)
    }
    fit <- tryCatch(
      .local_em(y, trials, family, initial_y, z_design, x_design, offset, u, u0, bandwidth, k, control,
                estimation_spec = estimation_spec,
                init_local = init_local, init_tau = init_tau, init_params = init_params),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      candidates[[start_id]] <- fit
      loglik[start_id] <- fit$loglik
      objective[start_id] <- fit$objective %||% fit$loglik
      start_converged[start_id] <- isTRUE(fit$converged)
    }
  }

  .select_local_start(candidates, loglik, objective, start_converged, u0,
                      legacy = isTRUE(control$legacy_em_rule))
}

# Choose among the EM starts at one grid point: the highest penalized objective
# among converged starts; when no start converged, the highest finite objective,
# which keeps `converged = FALSE`. `better_unconverged_start` flags a grid where
# an unconverged start reached a higher objective than the selected one.
.select_local_start <- function(candidates, loglik, objective, start_converged, u0, legacy = FALSE) {
  finite <- is.finite(objective)
  if (!any(finite)) {
    stop("All EM starts failed at u = ", signif(u0, 4), ".", call. = FALSE)
  }
  pool <- if (any(finite & start_converged)) which(finite & start_converged) else which(finite)
  best <- pool[which.max(objective[pool])]
  if (legacy) {
    best <- which.max(loglik)
  }
  out <- candidates[[best]]
  out$start_loglik <- loglik
  out$start_objective <- objective
  out$start_converged <- start_converged
  out$selected_start <- best
  out$better_unconverged_start <- any(finite & !start_converged & objective > objective[best] + 1e-8)
  out
}

.nearest_grid_assignment <- function(u, u_grid) {
  u <- as.numeric(u)
  u_grid <- sort(unique(as.numeric(u_grid)))
  if (length(u_grid) == 1L) {
    return(rep(1L, length(u)))
  }
  mids <- (u_grid[-1L] + u_grid[-length(u_grid)]) / 2
  findInterval(u, mids) + 1L
}

.joint_path_cost_guard <- function(n, n_grid, control) {
  n <- as.integer(n)
  n_grid <- as.integer(n_grid)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a positive finite observation count.", call. = FALSE)
  }
  if (length(n_grid) != 1L || !is.finite(n_grid) || n_grid < 1L) {
    stop("`n_grid` must be a positive finite grid-point count.", call. = FALSE)
  }
  dense_limit <- suppressWarnings(as.integer(
    control$joint_path_dense_grid_limit %||% 100L
  ))
  if (length(dense_limit) != 1L || !is.finite(dense_limit) || dense_limit < 1L) {
    dense_limit <- 100L
  }
  dense_ratio <- suppressWarnings(as.numeric(
    control$joint_path_dense_grid_ratio %||% 0.2
  ))
  if (length(dense_ratio) != 1L || !is.finite(dense_ratio) || dense_ratio <= 0) {
    dense_ratio <- 0.2
  }

  too_dense <- n_grid > dense_limit || n_grid / max(n, 1L) > dense_ratio
  if (too_dense && !isTRUE(control$allow_dense_u_grid)) {
    stop(
      sprintf(
        paste0(
          "`u_grid` has %d point(s) for %d observation(s), exceeding the ",
          "joint-path dense-grid guard. Set `control$allow_dense_u_grid = TRUE` ",
          "only when this computational cost is intended."
        ),
        n_grid,
        n
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.joint_path_unique_grid_warning <- function(u_values, u_grid) {
  u_grid <- sort(unique(as.numeric(u_grid)))
  unique_u <- sort(unique(as.numeric(u_values)))
  if (length(unique_u) == length(u_grid) && length(u_grid) > 20L &&
      max(abs(unique_u - sort(u_grid))) <= sqrt(.Machine$double.eps)) {
    warning(
      "`u_grid` is effectively `unique(u)`; this fits almost one local model per observation.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.joint_path_dense_grid_guard <- function(u_values, u_grid, control) {
  .joint_path_cost_guard(length(u_values), length(u_grid), control)
  .joint_path_unique_grid_warning(u_values, u_grid)
  invisible(TRUE)
}

.joint_path_param_vector <- function(params) {
  as.numeric(unlist(
    list(
      params$expert_coef,
      params$gating_coef,
      params$sigma,
      params$sigma_slope,
      params$theta
    ),
    use.names = FALSE
  ))
}

.joint_path_param_delta <- function(new_params, old_params) {
  if (is.null(old_params)) {
    return(NA_real_)
  }
  new_vec <- .joint_path_param_vector(new_params)
  old_vec <- .joint_path_param_vector(old_params)
  if (length(new_vec) != length(old_vec)) {
    return(NA_real_)
  }
  delta <- abs(new_vec - old_vec)
  if (!any(is.finite(delta))) {
    return(NA_real_)
  }
  max(delta[is.finite(delta)])
}

.joint_path_sample_loglik <- function(y, trials, family, z_design, x_design,
                                      offset, u, u_grid, bandwidth, assignment,
                                      params_list, control, estimation_spec) {
  total <- 0
  for (grid_id in seq_along(u_grid)) {
    rows <- assignment == grid_id
    if (!any(rows)) {
      next
    }
    du <- .local_du(u, u_grid[[grid_id]], bandwidth, estimation_spec)
    total <- total + .local_loglik(
      y,
      trials,
      family,
      z_design,
      x_design,
      du,
      offset,
      as.numeric(rows),
      params_list[[grid_id]],
      control
    )
  }
  as.numeric(total)
}

# Posterior of a u-free mixture of experts (constant kernel weights, local slope
# columns zero). As the first joint-path start it gives component labels that
# agree along u, which random and k-means starts often do not; those can
# converge to paths whose components swap partway along u.
.joint_path_global_tau <- function(y, trials, family, initial_y, z_design, x_design,
                                   offset, k, control, estimation_spec) {
  n <- length(y)
  global <- tryCatch(
    .fit_local_best(y, trials, family, initial_y, z_design, x_design, offset,
                    rep(0, n), 0, 1, k, control, previous = NULL,
                    estimation_spec = estimation_spec),
    error = function(e) NULL
  )
  posterior <- global$posterior
  if (is.null(posterior) || !all(is.finite(posterior))) {
    return(NULL)
  }
  posterior <- pmax(posterior, control$min_component_weight %||% 1e-8)
  posterior / rowSums(posterior)
}

.joint_path_em_one_start <- function(y, trials, family, initial_y, z_design, x_design,
                                     offset, u, u_grid, bandwidth, k, control,
                                     estimation_spec, assignment, progress = NULL,
                                     start_id = 1L, n_starts = 1L) {
  n <- length(y)
  n_grid <- length(u_grid)
  path_posterior <- NULL
  if (start_id == 1L && isTRUE(control$joint_path_global_start)) {
    path_posterior <- .joint_path_global_tau(y, trials, family, initial_y, z_design, x_design,
                                             offset, k, control, estimation_spec)
  }
  path_posterior <- path_posterior %||% .initial_tau(initial_y, k, random = start_id > 1L)
  params_list <- vector("list", n_grid)
  loglik <- rep(NA_real_, n_grid)
  trace <- data.frame(
    iteration = integer(0L),
    objective = numeric(0L),
    objective_delta = numeric(0L),
    posterior_delta = numeric(0L),
    parameter_delta = numeric(0L),
    converged = logical(0L)
  )
  objective <- -Inf
  converged <- FALSE
  posterior_by_grid <- array(NA_real_, dim = c(n, k, n_grid))
  start_start_time <- proc.time()[["elapsed"]]

  for (iter in seq_len(control$maxit)) {
    previous_posterior <- path_posterior
    previous_params <- params_list
    max_param_delta <- NA_real_
    objective_old <- objective

    for (grid_id in seq_along(u_grid)) {
      weights <- .local_weights(u, u_grid[[grid_id]], bandwidth, estimation_spec)
      du <- .local_du(u, u_grid[[grid_id]], bandwidth, estimation_spec)
      previous <- params_list[[grid_id]]
      if (is.null(previous)) {
        previous <- .initial_params_from_tau(
          y,
          trials,
          family,
          z_design,
          x_design,
          du,
          offset,
          weights,
          path_posterior,
          control
        )
      }
      params <- .m_step(
        y,
        trials,
        family,
        z_design,
        x_design,
        du,
        offset,
        weights,
        path_posterior,
        previous,
        control
      )
      params_list[[grid_id]] <- params
      posterior_by_grid[, , grid_id] <- .e_step(
        y,
        trials,
        family,
        z_design,
        x_design,
        du,
        offset,
        params,
        control
      )
      loglik[[grid_id]] <- .local_loglik(
        y,
        trials,
        family,
        z_design,
        x_design,
        du,
        offset,
        weights,
        params,
        control
      )
      grid_delta <- .joint_path_param_delta(params, previous_params[[grid_id]])
      if (is.finite(grid_delta)) {
        max_param_delta <- if (is.finite(max_param_delta)) {
          max(max_param_delta, grid_delta)
        } else {
          grid_delta
        }
      }
    }

    for (grid_id in seq_along(u_grid)) {
      rows <- assignment == grid_id
      if (any(rows)) {
        path_posterior[rows, ] <- matrix(
          posterior_by_grid[rows, , grid_id],
          ncol = k,
          byrow = FALSE
        )
      }
    }
    path_posterior <- pmax(path_posterior, 1e-12)
    path_posterior <- path_posterior / rowSums(path_posterior)

    objective <- .joint_path_sample_loglik(
      y,
      trials,
      family,
      z_design,
      x_design,
      offset,
      u,
      u_grid,
      bandwidth,
      assignment,
      params_list,
      control,
      estimation_spec
    )
    objective_delta <- if (is.finite(objective_old)) objective - objective_old else NA_real_
    posterior_delta <- max(abs(path_posterior - previous_posterior), na.rm = TRUE)
    did_converge <- is.finite(posterior_delta) &&
      posterior_delta <= control$tol &&
      (!is.finite(max_param_delta) || max_param_delta <= sqrt(control$tol))
    trace <- rbind(
      trace,
      data.frame(
        iteration = iter,
        objective = objective,
        objective_delta = objective_delta,
        posterior_delta = posterior_delta,
        parameter_delta = max_param_delta,
        converged = did_converge
      )
    )
    if (.vcmoe_progress_due(control, iter, did_converge)) {
      .vcmoe_progress_log(
        progress,
        event = "joint_path_iteration",
        status = if (did_converge) "converged" else "running",
        family = family,
        k = k,
        n = n,
        bandwidth = bandwidth,
        n_grid = n_grid,
        start_id = start_id,
        n_starts = n_starts,
        iteration = iter,
        maxit = control$maxit,
        start_start_time = start_start_time,
        loglik = objective,
        loglik_delta = objective_delta,
        converged = did_converge
      )
    }
    if (did_converge) {
      converged <- TRUE
      break
    }
  }

  iterations <- if (nrow(trace)) utils::tail(trace$iteration, 1L) else 0L
  locals <- vector("list", n_grid)
  for (grid_id in seq_along(u_grid)) {
    du <- .local_du(u, u_grid[[grid_id]], bandwidth, estimation_spec)
    weights <- .local_weights(u, u_grid[[grid_id]], bandwidth, estimation_spec)
    posterior <- .e_step(
      y,
      trials,
      family,
      z_design,
      x_design,
      du,
      offset,
      params_list[[grid_id]],
      control
    )
    locals[[grid_id]] <- .params_to_local(
      params_list[[grid_id]],
      posterior,
      z_design,
      x_design,
      du,
      .local_loglik(
        y,
        trials,
        family,
        z_design,
        x_design,
        du,
        offset,
        weights,
        params_list[[grid_id]],
        control
      ),
      iterations,
      converged
    )
    locals[[grid_id]]$selected_start <- start_id
    locals[[grid_id]]$start_loglik <- objective
  }

  list(
    locals = locals,
    path_posterior = path_posterior,
    assignment = assignment,
    loglik = objective,
    loglik_by_grid = loglik,
    trace = trace,
    converged = converged,
    iterations = iterations,
    selected_start = start_id
  )
}

.joint_path_fit_best <- function(y, trials, family, initial_y, z_design, x_design,
                                 offset, u, u_grid, bandwidth, k, control,
                                 estimation_spec, assignment, progress = NULL) {
  starts <- max(1L, as.integer(control$n_starts))
  candidates <- vector("list", starts)
  loglik <- rep(-Inf, starts)
  errors <- rep("", starts)

  for (start_id in seq_len(starts)) {
    start_start_time <- proc.time()[["elapsed"]]
    error_message <- NULL
    fit <- tryCatch(
      .joint_path_em_one_start(
        y,
        trials,
        family,
        initial_y,
        z_design,
        x_design,
        offset,
        u,
        u_grid,
        bandwidth,
        k,
        control,
        estimation_spec,
        assignment,
        progress = progress,
        start_id = start_id,
        n_starts = starts
      ),
      error = function(e) {
        error_message <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(fit) && is.finite(fit$loglik)) {
      candidates[[start_id]] <- fit
      loglik[[start_id]] <- fit$loglik
      .vcmoe_progress_log(
        progress,
        event = "joint_path_start_finished",
        status = "completed",
        family = family,
        k = k,
        n = length(y),
        bandwidth = bandwidth,
        n_grid = length(u_grid),
        start_id = start_id,
        n_starts = starts,
        maxit = control$maxit,
        start_start_time = start_start_time,
        loglik = fit$loglik,
        converged = fit$converged
      )
    } else {
      errors[[start_id]] <- error_message %||% "non-finite joint-path objective"
      .vcmoe_progress_log(
        progress,
        event = "joint_path_start_failed",
        status = "failed",
        family = family,
        k = k,
        n = length(y),
        bandwidth = bandwidth,
        n_grid = length(u_grid),
        start_id = start_id,
        n_starts = starts,
        maxit = control$maxit,
        start_start_time = start_start_time,
        error_message = errors[[start_id]]
      )
    }
  }

  if (!any(is.finite(loglik))) {
    detail <- unique(errors[nzchar(errors)])
    suffix <- if (length(detail)) paste0(" ", paste(detail, collapse = "; ")) else ""
    stop("All joint-path EM starts failed.", suffix, call. = FALSE)
  }
  selected <- which.max(loglik)
  out <- candidates[[selected]]
  out$selected_start <- selected
  out$start_loglik <- loglik
  out
}

.vcmoe_fit_joint_path <- function(formula, data, u, k, family, bandwidth, u_grid,
                                  control, label, parameterization, u_scale,
                                  progress, fit_call) {
  if (is.list(formula)) {
    stop("Joint-path EM supports single-response formulas only.", call. = FALSE)
  }
  family <- match.arg(family, c("gaussian", "binomial", "negative-binomial"))
  label <- match.arg(label, c("align", "global", "greedy"))
  parameterization <- match.arg(parameterization, "a1_epanechnikov_scaled")
  u_scale <- .vcmoe_validate_u_scale(u_scale)
  estimation_spec <- .vcmoe_estimation_spec(parameterization)
  k <- as.integer(k)
  if (k < 2L || k > 10L) {
    stop("`k` must be between 2 and 10.", call. = FALSE)
  }
  alignment_method <- if (identical(label, "greedy")) {
    "greedy"
  } else if (identical(label, "global") && k <= 6L) {
    "global"
  } else if (k <= 6L) {
    "global"
  } else {
    "sequential"
  }

  user_control <- control %||% list()
  ridge_user_supplied <- "ridge" %in% names(user_control)
  control <- .vcmoe_joint_path_control(user_control)
  if (!is.null(control$seed)) {
    set.seed(control$seed)
  }

  pieces <- .split_vcmoe_formula(formula)
  u_info <- .extract_u(u, data)
  u_values_all <- as.numeric(u_info$values)
  expert_frame_all <- stats::model.frame(
    pieces$expert_formula,
    data = data,
    na.action = stats::na.pass
  )
  gating_frame_all <- stats::model.frame(
    pieces$gating_formula,
    data = data,
    na.action = stats::na.pass
  )
  expert_offset_all <- .extract_model_offset(expert_frame_all)
  gating_offset_all <- .extract_model_offset(gating_frame_all)
  if (any(abs(gating_offset_all) > 0)) {
    stop("Gating-side `offset()` terms are not supported in v0.", call. = FALSE)
  }
  if (!identical(family, "negative-binomial") && any(abs(expert_offset_all) > 0)) {
    stop("Expert-side `offset()` terms are currently supported only for `family = \"negative-binomial\"`.", call. = FALSE)
  }
  complete_rows <- .complete_cases_frame(expert_frame_all, nrow(data)) &
    .complete_cases_frame(gating_frame_all, nrow(data)) &
    is.finite(u_values_all) &
    is.finite(expert_offset_all)
  if (!any(complete_rows)) {
    stop("No complete rows remain after checking response, covariates, and `u`.", call. = FALSE)
  }
  if (!all(complete_rows)) {
    warning(
      sprintf("Removed %d row(s) with missing or non-finite response, covariates, or `u`.",
              sum(!complete_rows)),
      call. = FALSE
    )
  }

  expert_frame <- expert_frame_all[complete_rows, , drop = FALSE]
  gating_frame <- gating_frame_all[complete_rows, , drop = FALSE]
  response <- .parse_vcmoe_response(family, stats::model.response(expert_frame), pieces$response)
  if (identical(family, "binomial") &&
      identical(response$type, "bernoulli") &&
      !ridge_user_supplied) {
    control$ridge <- 1
  }
  y <- response$y
  trials <- response$trials
  initial_y <- response$observed
  z_design <- stats::model.matrix(stats::terms(pieces$expert_formula), expert_frame)
  x_design <- stats::model.matrix(stats::terms(pieces$gating_formula), data = gating_frame)
  expert_offset <- expert_offset_all[complete_rows]
  u_original <- u_values_all[complete_rows]
  u_scaling <- .vcmoe_make_u_scaling(u_original, u_scale)
  u_values <- .vcmoe_scale_u(u_original, u_scaling)

  if (is.null(bandwidth)) {
    bandwidth <- .default_bandwidth(u_values)
  }
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0) {
    stop("`bandwidth` must be a positive finite number.", call. = FALSE)
  }
  if (is.null(u_grid)) {
    u_grid <- .default_u_grid(u_values)
  }
  u_grid <- sort(unique(as.numeric(u_grid)))
  if (!length(u_grid)) {
    stop("`u_grid` must contain at least one value.", call. = FALSE)
  }
  if (any(!is.finite(u_grid))) {
    stop("`u_grid` must contain finite values on the analysis scale.", call. = FALSE)
  }
  if (identical(u_scale, "unit") && (min(u_grid) < -1e-8 || max(u_grid) > 1 + 1e-8)) {
    warning(
      "`u_scale = \"unit\"` interprets `u_grid` on the scaled [0, 1] analysis scale.",
      call. = FALSE
    )
  }
  .joint_path_dense_grid_guard(u_values, u_grid, control)
  u_grid_original <- .vcmoe_unscale_u(u_grid, u_scaling)

  n_grid <- length(u_grid)
  q <- ncol(z_design)
  p <- ncol(x_design)
  assignment <- .nearest_grid_assignment(u_values, u_grid)
  assignment_counts <- tabulate(assignment, nbins = n_grid)
  progress_state <- .vcmoe_progress_setup(progress)
  warning(
    paste0(
      "Joint-path EM fits one local model per grid point in each EM iteration; ",
      "runtime scales approximately with `control$n_starts * control$maxit * ",
      "length(u_grid) * n`."
    ),
    call. = FALSE
  )
  .vcmoe_progress_log(
    progress_state,
    event = "fit_started",
    status = "running",
    family = family,
    k = k,
    n = length(y),
    bandwidth = bandwidth,
    n_grid = n_grid,
    n_starts = control$n_starts,
    maxit = control$maxit
  )
  joint <- .joint_path_fit_best(
    y,
    trials,
    family,
    initial_y,
    z_design,
    x_design,
    expert_offset,
    u_values,
    u_grid,
    bandwidth,
    k,
    control,
    estimation_spec,
    assignment,
    progress = progress_state
  )
  if (!isTRUE(joint$converged)) {
    warning(
      sprintf(
        paste0(
          "Joint-path EM did not converge within %d outer iteration(s); ",
          "inspect `fit$diagnostics$joint_path_trace` before interpretation."
        ),
        control$maxit
      ),
      call. = FALSE
    )
  }
  locals <- joint$locals

  previous <- NULL
  for (grid_id in seq_along(u_grid)) {
    locals[[grid_id]] <- .align_local_fit(
      locals[[grid_id]],
      previous,
      z_design,
      x_design,
      initial_y,
      control
    )
    previous <- locals[[grid_id]]
  }
  diagnostics <- list(
    engine = .vcmoe_joint_path_engine_id,
    alignment_method = alignment_method,
    loglik = numeric(n_grid),
    converged = logical(n_grid),
    iterations = integer(n_grid),
    selected_start = integer(n_grid),
    permutations = matrix(NA_integer_, nrow = n_grid, ncol = k),
    greedy_permutations = matrix(NA_integer_, nrow = n_grid, ncol = k),
    global_permutations = matrix(NA_integer_, nrow = n_grid, ncol = k),
    expert_optimizer_convergence = matrix(NA_integer_, nrow = n_grid, ncol = k),
    expert_gradient_norm = matrix(NA_real_, nrow = n_grid, ncol = k),
    expert_coef_norm = matrix(NA_real_, nrow = n_grid, ncol = k),
    global_path_cost = NA_real_,
    global_transition_cost = rep(NA_real_, n_grid),
    alignment_margin = numeric(n_grid),
    ambiguous = logical(n_grid),
    posterior_entropy = numeric(n_grid),
    warnings = character(0L),
    start_loglik = vector("list", n_grid),
    joint_path_iterations = joint$iterations,
    joint_path_converged = joint$converged,
    joint_path_objective = "nearest_grid_sample_loglik",
    joint_path_trace = joint$trace,
    joint_path_assignment = data.frame(
      grid_id = seq_len(n_grid),
      u = u_grid,
      u_original = u_grid_original,
      n_assigned = assignment_counts
    ),
    joint_path_selected_start = joint$selected_start,
    joint_path_start_loglik = joint$start_loglik
  )

  if (identical(alignment_method, "global")) {
    global_alignment <- .global_align_local_fits(
      locals,
      u_grid,
      z_design,
      x_design,
      initial_y,
      control,
      bandwidth = bandwidth,
      estimation_spec = estimation_spec
    )
    locals <- global_alignment$locals
    diagnostics$global_path_cost <- global_alignment$path_cost
    diagnostics$global_transition_cost <- global_alignment$transition_cost
  }

  expert <- array(NA_real_, dim = c(n_grid, k, q),
                  dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k)),
                                  term = colnames(z_design)))
  expert_slope <- expert
  gating <- array(NA_real_, dim = c(n_grid, k, p),
                  dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k)),
                                  term = colnames(x_design)))
  gating_slope <- gating
  sigma <- if (identical(family, "gaussian")) {
    matrix(NA_real_, nrow = n_grid, ncol = k,
           dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k))))
  } else {
    NULL
  }
  sigma_slope <- if (identical(family, "gaussian")) {
    matrix(NA_real_, nrow = n_grid, ncol = k,
           dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k))))
  } else {
    NULL
  }
  theta <- if (identical(family, "negative-binomial")) {
    matrix(NA_real_, nrow = n_grid, ncol = k,
           dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k))))
  } else {
    NULL
  }
  posterior <- array(NA_real_, dim = c(length(y), k, n_grid),
                     dimnames = list(observation = NULL, component = paste0("component", seq_len(k)),
                                     u = signif(u_grid, 8)))
  path_posterior <- matrix(
    NA_real_,
    nrow = length(y),
    ncol = k,
    dimnames = list(observation = NULL, component = paste0("component", seq_len(k)))
  )

  for (grid_id in seq_along(u_grid)) {
    local <- locals[[grid_id]]
    diagnostics$ambiguous[grid_id] <- isTRUE(local$ambiguous)
    if (isTRUE(local$ambiguous)) {
      diagnostics$warnings <- c(
        diagnostics$warnings,
        sprintf("Ambiguous label alignment near u = %.4f.", u_grid[grid_id])
      )
    }

    expert[grid_id, , ] <- local$expert_intercept
    expert_slope[grid_id, , ] <- local$expert_slope
    gating[grid_id, , ] <- local$gating_intercept
    gating_slope[grid_id, , ] <- local$gating_slope
    if (!is.null(sigma)) {
      sigma[grid_id, ] <- local$sigma
    }
    if (!is.null(sigma_slope)) {
      sigma_slope[grid_id, ] <- local$sigma_slope
    }
    if (!is.null(theta)) {
      theta[grid_id, ] <- local$theta
    }
    posterior[, , grid_id] <- local$posterior
    rows <- assignment == grid_id
    if (any(rows)) {
      path_posterior[rows, ] <- matrix(
        local$posterior[rows, ],
        ncol = k,
        byrow = FALSE
      )
    }

    diagnostics$loglik[grid_id] <- local$loglik
    diagnostics$converged[grid_id] <- local$converged
    diagnostics$iterations[grid_id] <- local$iterations
    diagnostics$selected_start[grid_id] <- joint$selected_start
    diagnostics$permutations[grid_id, ] <- local$permutation
    diagnostics$greedy_permutations[grid_id, ] <- local$greedy_permutation %||% local$permutation
    diagnostics$global_permutations[grid_id, ] <- local$global_permutation %||% seq_len(k)
    diagnostics$expert_optimizer_convergence[grid_id, ] <- local$expert_convergence %||% rep(NA_integer_, k)
    diagnostics$expert_gradient_norm[grid_id, ] <- local$expert_gradient_norm %||% rep(NA_real_, k)
    diagnostics$expert_coef_norm[grid_id, ] <- local$expert_coef_norm %||% rep(NA_real_, k)
    diagnostics$alignment_margin[grid_id] <- local$alignment_margin
    diagnostics$posterior_entropy[grid_id] <- .weighted_entropy(
      local$posterior,
      .local_weights(u_values, u_grid[grid_id], bandwidth, estimation_spec)
    )
    diagnostics$start_loglik[[grid_id]] <- joint$start_loglik
  }

  if (length(diagnostics$warnings) && isTRUE(control$warn_ambiguous)) {
    warning(
      sprintf("Label alignment was ambiguous at %d grid point(s); inspect `fit$diagnostics$warnings`.",
              length(diagnostics$warnings)),
      call. = FALSE
    )
  }
  .vcmoe_progress_log(
    progress_state,
    event = "fit_finished",
    status = "completed",
    family = family,
    k = k,
    n = length(y),
    bandwidth = bandwidth,
    n_grid = n_grid,
    n_starts = control$n_starts,
    maxit = control$maxit,
    selected_start = joint$selected_start,
    converged = isTRUE(joint$converged)
  )

  fit_call$progress <- NULL
  parameterization_meta <- .vcmoe_parameterization_metadata(
    family = family,
    k = k,
    bandwidth = bandwidth,
    label = label,
    alignment_method = alignment_method,
    control = control,
    u_scaling = u_scaling,
    parameterization = parameterization,
    estimation_spec = estimation_spec
  )
  parameterization_meta$estimator <- .vcmoe_joint_path_engine_id
  parameterization_meta$engine <- list(
    id = .vcmoe_joint_path_engine_id,
    responsibility_update = "nearest_grid_observation_path",
    diagnostic_objective = "nearest_grid_sample_loglik",
    runtime = "Fits one local model per grid point in each EM iteration."
  )

  object <- list(
    call = fit_call,
    formula = formula,
    family = family,
    k = k,
    label = alignment_method,
    engine_id = .vcmoe_joint_path_engine_id,
    parameterization_id = parameterization,
    u_scale = u_scale,
    u_scaling = u_scaling,
    bandwidth = bandwidth,
    u_grid = u_grid,
    u_grid_original = u_grid_original,
    control = control,
    parameterization = parameterization_meta,
    coefficients = list(
      expert = expert,
      expert_slope = expert_slope,
      gating = gating,
      gating_slope = gating_slope,
      sigma = sigma,
      sigma_slope = sigma_slope,
      theta = theta
    ),
    diagnostics = diagnostics,
    posterior = posterior,
    joint_path_posterior = path_posterior,
    terms = list(
      expert = stats::terms(pieces$expert_formula),
      gating = stats::terms(pieces$gating_formula)
    ),
    response = pieces$response,
    response_info = response$info,
    u_name = u_info$name,
    rows_used = which(complete_rows),
    rows_removed = sum(!complete_rows),
    fitted = if (isTRUE(control$keep_data)) {
      list(
        y = response$observed,
        success = response$success,
        trials = response$trials,
        response = response,
        z_design = z_design,
        x_design = x_design,
        expert_offset = expert_offset,
        u = u_values,
        u_original = u_original
      )
    } else {
      NULL
    }
  )
  class(object) <- "vcmoe"
  object
}

#' Fit a varying-coefficient mixture-of-experts model
#'
#' @param formula A formula of the form `y ~ expert_terms | gating_terms`.
#'   For grouped binomial data, use `cbind(success, failure) ~ expert_terms |
#'   gating_terms`. For Negative-Binomial count data, expert-side
#'   `offset(log_size_factor)` terms are supported.
#' @param data A data frame.
#' @param u Continuous index column name or numeric vector.
#' @param k Number of mixture components. `k = 2` is the primary v0 target;
#'   `k = 3:10` are high-k candidate support and require diagnostics before
#'   being treated as stable.
#' @param family Model family. `"gaussian"`, `"binomial"`, and
#'   `"negative-binomial"` are implemented.
#' @param bandwidth Kernel bandwidth. If `NULL`, a Silverman-style default is
#'   used for `u`.
#' @param u_grid Grid where coefficient functions are estimated.
#' @param control Named list overriding EM and label-alignment settings.
#'   Common entries: `maxit` (default 100), `tol` (default `1e-5`), `n_starts`
#'   (default 3) and `seed`. The local EM stops when the relative change of
#'   the penalized local objective, `|delta F| / (|F| + 1)`, falls below `tol`;
#'   with large samples this can stop well before the local maximum, so check
#'   `vcmoe_diagnostics()` and use [vcmoe_refine()] before inference.
#'   Negative-Binomial: `negbin_theta_ridge` and `negbin_theta_target` add a
#'   ridge on `log(theta)` towards the target, scaled by the kernel mass of
#'   the grid point (`negbin_theta_penalty_mass = "kernel"`, the default;
#'   `"component"` is the rule of VCMoE 0.2.0). Grid-parallel fitting
#'   (local-grid engine): `grid_subset` fits only the given grid indices and
#'   returns a `vcmoe_grid_part` to merge with [vcmoe_merge_grid_fits()];
#'   `grid_seed = "per_grid"` reseeds each grid point with `seed + index`;
#'   `grid_warm_start = FALSE` drops the warm start from the previous grid
#'   point; `grid_init_locals` replaces the starts at given grid points (see
#'   [vcmoe_grid_neighbor_inits()]). The defaults reproduce the serial fit.
#'   Known labels (local-grid engine): `fixed_membership` gives each row of
#'   `data` a component, as indices in `1..k`, a factor (levels in component
#'   order) or an `nrow(data) x k` 0/1 matrix. Each local fit is then the
#'   kernel-weighted complete-data fit (component-wise experts and a
#'   multinomial-logit gate on the labels), with one start and no label
#'   alignment, and inference uses the complete-data local likelihood. Joint
#'   path: `joint_path_global_start` (default `TRUE`) makes the first start
#'   the posterior of a mixture of experts without `u` dependence, which
#'   keeps component labels consistent along `u`.
#' @param label Label strategy. `"align"` uses exact global alignment for
#'   `k <= 6` and sequential assignment with ambiguity margins for `k >= 7`;
#'   `"global"` requests exact global alignment when feasible and falls back
#'   to the same sequential assignment path for `k >= 7`;
#'   `"greedy"` keeps the older one-step alignment.
#' @param u_scale How to transform `u` before fitting. The default `"unit"`
#'   maps complete-row `u` values to `[0, 1]`; `"none"` leaves `u` on the
#'   supplied scale. `bandwidth` and `u_grid` are interpreted on the transformed
#'   analysis scale.
#' @param parameterization Estimator convention. The public package uses
#'   `"a1_epanechnikov_scaled"`: Epanechnikov density weights `K((u-u0)/h)/h`
#'   and the scaled local-linear basis `(u-u0)/h`.
#' @param engine Fitting engine. The default `"local_grid_em"` preserves the
#'   original independent local-grid EM path. `"joint_path_em"` updates a
#'   shared observation-level responsibility path across grid points.
#' @param progress Joint-path progress reporting. The default `NULL` is silent;
#'   use `TRUE` for messages or a single CSV file path for structured logging.
#'   Iteration frequency is controlled by `control$progress_every`.
#' @return An object of class `vcmoe`, or a `vcmoe_grid_part` when
#'   `control$grid_subset` is set.
#' @details At each grid point, starts that met the stopping rule are
#'   preferred and the one with the highest penalized local objective is kept;
#'   when no start met it, the best finite start is kept with
#'   `converged = FALSE`. `fit$diagnostics` records, per grid point, the stop
#'   reason, the last objective change, the selected start and whether an
#'   unconverged start reached a higher objective. With global label
#'   alignment, a grid point relabelled against a confident match with its
#'   neighbour is flagged in `alignment_override` with a warning; compare with
#'   `label = "greedy"`.
#'
#'   Rows with missing or non-finite response, covariates, or `u` are
#'   removed consistently before fitting, with a warning. For single-trial
#'   Bernoulli responses, the default gating ridge is strengthened to
#'   `control$ridge = 1` unless the user explicitly supplies `control$ridge`;
#'   grouped Binomial and other families keep the global default. Joint-path
#'   traces record the sample-level nearest-grid log-likelihood as a diagnostic
#'   criterion. The label-consistent updates do not guarantee that this
#'   diagnostic is monotone; convergence is based on posterior and parameter
#'   deltas instead.
#' @export
vcmoe_fit <- function(formula, data, u, k = 2L, family = "gaussian",
                      bandwidth = NULL, u_grid = NULL, control = list(),
                      label = "align",
                      parameterization = "a1_epanechnikov_scaled",
                      u_scale = c("unit", "none"),
                      engine = c("local_grid_em", "joint_path_em"),
                      progress = NULL) {
  engine <- match.arg(engine)
  if (identical(engine, .vcmoe_joint_path_engine_id)) {
    if (!is.null(control$fixed_membership)) {
      stop("`fixed_membership` is supported by `engine = \"local_grid_em\"` only.", call. = FALSE)
    }
    return(.vcmoe_fit_joint_path(
      formula = formula,
      data = data,
      u = u,
      k = k,
      family = family,
      bandwidth = bandwidth,
      u_grid = u_grid,
      control = control,
      label = label,
      parameterization = parameterization,
      u_scale = u_scale,
      progress = progress,
      fit_call = match.call()
    ))
  }
  if (!is.null(progress) && !identical(progress, FALSE)) {
    stop("`progress` is available only for `engine = \"joint_path_em\"`.", call. = FALSE)
  }
  family <- match.arg(family, c("gaussian", "binomial", "negative-binomial"))
  label <- match.arg(label, c("align", "global", "greedy"))
  parameterization <- match.arg(parameterization, "a1_epanechnikov_scaled")
  u_scale <- .vcmoe_validate_u_scale(u_scale)
  estimation_spec <- .vcmoe_estimation_spec(parameterization)
  k <- as.integer(k)
  if (k < 2L || k > 10L) {
    stop("`k` must be between 2 and 10.", call. = FALSE)
  }
  alignment_method <- .vcmoe_alignment_method(label, k)
  if (identical(family, "binomial")) {
    if (k >= 3L) {
      warning(
        "Binomial VCMoE with k >= 3 is experimental high-k candidate support; inspect diagnostics before interpretation.",
        call. = FALSE
      )
    }
  } else if (identical(family, "negative-binomial") && k >= 3L) {
    warning(
      "Negative-Binomial VCMoE with k >= 3 is high-k candidate support; stable status depends on validation diagnostics.",
      call. = FALSE
    )
  }

  user_control <- control %||% list()
  ridge_user_supplied <- "ridge" %in% names(user_control)
  control <- .vcmoe_default_control(user_control)
  if (!is.null(control$seed)) {
    set.seed(control$seed)
  }

  pieces <- .split_vcmoe_formula(formula)
  u_info <- .extract_u(u, data)
  u_values_all <- as.numeric(u_info$values)
  expert_frame_all <- stats::model.frame(
    pieces$expert_formula,
    data = data,
    na.action = stats::na.pass
  )
  gating_frame_all <- stats::model.frame(
    pieces$gating_formula,
    data = data,
    na.action = stats::na.pass
  )
  expert_offset_all <- .extract_model_offset(expert_frame_all)
  gating_offset_all <- .extract_model_offset(gating_frame_all)
  if (any(abs(gating_offset_all) > 0)) {
    stop("Gating-side `offset()` terms are not supported in v0.", call. = FALSE)
  }
  if (!identical(family, "negative-binomial") && any(abs(expert_offset_all) > 0)) {
    stop("Expert-side `offset()` terms are currently supported only for `family = \"negative-binomial\"`.", call. = FALSE)
  }
  complete_rows <- .complete_cases_frame(expert_frame_all, nrow(data)) &
    .complete_cases_frame(gating_frame_all, nrow(data)) &
    is.finite(u_values_all) &
    is.finite(expert_offset_all)
  if (!any(complete_rows)) {
    stop("No complete rows remain after checking response, covariates, and `u`.", call. = FALSE)
  }
  if (!all(complete_rows)) {
    warning(
      sprintf("Removed %d row(s) with missing or non-finite response, covariates, or `u`.",
              sum(!complete_rows)),
      call. = FALSE
    )
  }

  expert_frame <- expert_frame_all[complete_rows, , drop = FALSE]
  gating_frame <- gating_frame_all[complete_rows, , drop = FALSE]
  response <- .parse_vcmoe_response(family, stats::model.response(expert_frame), pieces$response)
  if (identical(family, "binomial") &&
      identical(response$type, "bernoulli") &&
      !ridge_user_supplied) {
    control$ridge <- 1
  }
  y <- response$y
  trials <- response$trials
  initial_y <- response$observed
  z_design <- stats::model.matrix(stats::terms(pieces$expert_formula), expert_frame)
  x_design <- stats::model.matrix(stats::terms(pieces$gating_formula), data = gating_frame)
  expert_offset <- expert_offset_all[complete_rows]
  u_original <- u_values_all[complete_rows]
  u_scaling <- .vcmoe_make_u_scaling(u_original, u_scale)
  u_values <- .vcmoe_scale_u(u_original, u_scaling)
  # Stored on the fitted rows; components keep the label order, so no alignment.
  if (!is.null(control$fixed_membership)) {
    control$fixed_membership <- .vcmoe_fixed_membership(control$fixed_membership, nrow(data), k)[
      complete_rows, , drop = FALSE]
    alignment_method <- "none"
  }

  if (is.null(bandwidth)) {
    bandwidth <- .default_bandwidth(u_values)
  }
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0) {
    stop("`bandwidth` must be a positive finite number.", call. = FALSE)
  }
  if (is.null(u_grid)) {
    u_grid <- .default_u_grid(u_values)
  }
  u_grid <- sort(unique(as.numeric(u_grid)))
  if (!length(u_grid)) {
    stop("`u_grid` must contain at least one value.", call. = FALSE)
  }
  if (any(!is.finite(u_grid))) {
    stop("`u_grid` must contain finite values on the analysis scale.", call. = FALSE)
  }
  if (identical(u_scale, "unit") && (min(u_grid) < -1e-8 || max(u_grid) > 1 + 1e-8)) {
    warning(
      "`u_scale = \"unit\"` interprets `u_grid` on the scaled [0, 1] analysis scale.",
      call. = FALSE
    )
  }
  u_grid_original <- .vcmoe_unscale_u(u_grid, u_scaling)

  n_grid <- length(u_grid)
  grid_options <- .vcmoe_grid_options(control, n_grid)
  grid_ids <- grid_options$subset %||% seq_len(n_grid)
  fit_call <- match.call()
  # Stored controls never carry the per-task subset or the warm-start fits.
  control$grid_subset <- NULL
  control$grid_init_locals <- NULL

  previous <- NULL
  previous_grid_id <- NA_integer_
  locals <- vector("list", n_grid)
  raw_locals <- vector("list", n_grid)
  for (grid_id in grid_ids) {
    if (identical(grid_options$seed, "per_grid")) {
      set.seed(control$seed + grid_id)
    }
    # Warm start only from the immediately preceding grid point fitted in this run.
    warm <- if (isTRUE(grid_options$warm_start) && identical(previous_grid_id, grid_id - 1L)) previous else NULL
    local <- .fit_local_best(y, trials, family, initial_y, z_design, x_design, expert_offset,
                             u_values, u_grid[grid_id],
                             bandwidth, k, control, warm, estimation_spec,
                             init_locals = grid_options$init_locals[[as.character(grid_id)]])
    raw_locals[[grid_id]] <- local
    local <- .align_local_fit(local, previous, z_design, x_design, initial_y, control)
    locals[[grid_id]] <- local
    previous <- local
    previous_grid_id <- grid_id
  }

  spec <- list(
    formula = formula, pieces = pieces, family = family, k = k, label = label,
    alignment_method = alignment_method, parameterization = parameterization,
    estimation_spec = estimation_spec, u_scale = u_scale, u_scaling = u_scaling,
    bandwidth = bandwidth, u_grid = u_grid, u_grid_original = u_grid_original, control = control,
    response = response, y = y, trials = trials, initial_y = initial_y, z_design = z_design,
    x_design = x_design, expert_offset = expert_offset, u_values = u_values, u_original = u_original,
    u_name = u_info$name, complete_rows = complete_rows, fit_call = fit_call
  )
  if (!is.null(grid_options$subset)) {
    part_locals <- raw_locals[grid_ids]
    names(part_locals) <- as.character(grid_ids)
    return(structure(
      list(spec = spec, grid_ids = grid_ids, locals = part_locals,
           grid_seed = grid_options$seed, grid_warm_start = grid_options$warm_start,
           second_pass = !is.null(grid_options$init_locals)),
      class = "vcmoe_grid_part"
    ))
  }
  .vcmoe_finalize_fit(spec, locals)
}


.vcmoe_alignment_method <- function(label, k) {
  if (identical(label, "greedy")) {
    "greedy"
  } else if (k <= 6L) {
    "global"
  } else {
    "sequential"
  }
}

.vcmoe_grid_control_names <- c("grid_subset", "grid_seed", "grid_warm_start", "grid_init_locals")

# Grid-parallel fitting. Defaults reproduce the serial fit: one global seed and
# a warm start from the previous grid point. `grid_subset` fits only those grid
# indices and returns a "vcmoe_grid_part"; `grid_init_locals` (named by grid
# index) replaces the starts at those grid points with the given local fits.
.vcmoe_grid_options <- function(control, n_grid) {
  seed_mode <- control$grid_seed %||% "global"
  if (!is.character(seed_mode) || length(seed_mode) != 1L || !seed_mode %in% c("global", "per_grid")) {
    stop("`grid_seed` must be \"global\" or \"per_grid\".", call. = FALSE)
  }
  if (identical(seed_mode, "per_grid") && is.null(control$seed)) {
    stop("`grid_seed = \"per_grid\"` requires `control$seed`.", call. = FALSE)
  }
  warm_start <- control$grid_warm_start %||% TRUE
  if (!is.logical(warm_start) || length(warm_start) != 1L || is.na(warm_start)) {
    stop("`grid_warm_start` must be TRUE or FALSE.", call. = FALSE)
  }
  subset <- control$grid_subset
  if (!is.null(subset)) {
    subset <- suppressWarnings(as.integer(subset))
    if (!length(subset) || anyNA(subset) || anyDuplicated(subset) || any(subset < 1L | subset > n_grid)) {
      stop(sprintf("`grid_subset` must be distinct grid indices in 1..%d.", n_grid), call. = FALSE)
    }
    subset <- sort(subset)
  }
  init_locals <- control$grid_init_locals
  if (!is.null(init_locals)) {
    if (is.null(subset)) {
      stop("`grid_init_locals` requires `grid_subset`.", call. = FALSE)
    }
    if (!is.list(init_locals) || is.null(names(init_locals)) ||
        !all(names(init_locals) %in% as.character(subset)) ||
        !all(vapply(init_locals, function(x) is.list(x) && length(x) > 0L, logical(1L)))) {
      stop("`grid_init_locals` must be a non-empty list of local fits per grid index in `grid_subset`.", call. = FALSE)
    }
  }
  list(seed = seed_mode, warm_start = warm_start, subset = subset, init_locals = init_locals)
}

# Post-loop half of vcmoe_fit(): global label alignment over the grid and
# assembly of the fit object from aligned local fits.
.vcmoe_finalize_fit <- function(spec, locals) {
  family <- spec$family
  k <- spec$k
  label <- spec$label
  alignment_method <- spec$alignment_method
  parameterization <- spec$parameterization
  estimation_spec <- spec$estimation_spec
  u_scale <- spec$u_scale
  u_scaling <- spec$u_scaling
  bandwidth <- spec$bandwidth
  u_grid <- spec$u_grid
  u_grid_original <- spec$u_grid_original
  control <- spec$control
  pieces <- spec$pieces
  formula <- spec$formula
  response <- spec$response
  y <- spec$y
  initial_y <- spec$initial_y
  z_design <- spec$z_design
  x_design <- spec$x_design
  expert_offset <- spec$expert_offset
  u_values <- spec$u_values
  u_original <- spec$u_original
  complete_rows <- spec$complete_rows
  u_info <- list(name = spec$u_name)
  n_grid <- length(u_grid)
  q <- ncol(z_design)
  p <- ncol(x_design)
  expert <- array(NA_real_, dim = c(n_grid, k, q),
                  dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k)),
                                  term = colnames(z_design)))
  expert_slope <- expert
  gating <- array(NA_real_, dim = c(n_grid, k, p),
                  dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k)),
                                  term = colnames(x_design)))
  gating_slope <- gating
  sigma <- if (identical(family, "gaussian")) {
    matrix(NA_real_, nrow = n_grid, ncol = k,
           dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k))))
  } else {
    NULL
  }
  sigma_slope <- if (identical(family, "gaussian")) {
    matrix(NA_real_, nrow = n_grid, ncol = k,
           dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k))))
  } else {
    NULL
  }
  theta <- if (identical(family, "negative-binomial")) {
    matrix(NA_real_, nrow = n_grid, ncol = k,
           dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k))))
  } else {
    NULL
  }
  posterior <- array(NA_real_, dim = c(length(y), k, n_grid),
                     dimnames = list(observation = NULL, component = paste0("component", seq_len(k)),
                                     u = signif(u_grid, 8)))

  diagnostics <- list(
    engine = "local_grid_em",
    alignment_method = alignment_method,
    loglik = numeric(n_grid),
    converged = logical(n_grid),
    iterations = integer(n_grid),
    selected_start = integer(n_grid),
    permutations = matrix(NA_integer_, nrow = n_grid, ncol = k),
    greedy_permutations = matrix(NA_integer_, nrow = n_grid, ncol = k),
    global_permutations = matrix(NA_integer_, nrow = n_grid, ncol = k),
    expert_optimizer_convergence = matrix(NA_integer_, nrow = n_grid, ncol = k),
    expert_gradient_norm = matrix(NA_real_, nrow = n_grid, ncol = k),
    expert_coef_norm = matrix(NA_real_, nrow = n_grid, ncol = k),
    global_path_cost = NA_real_,
    global_transition_cost = rep(NA_real_, n_grid),
    alignment_margin = numeric(n_grid),
    ambiguous = logical(n_grid),
    posterior_entropy = numeric(n_grid),
    warnings = character(0L),
    start_loglik = vector("list", n_grid),
    objective = rep(NA_real_, n_grid),
    objective_delta = rep(NA_real_, n_grid),
    stop_reason = rep(NA_character_, n_grid),
    start_objective = vector("list", n_grid),
    start_converged = vector("list", n_grid),
    better_unconverged_start = logical(n_grid),
    alignment_override = logical(n_grid)
  )

  if (identical(alignment_method, "global")) {
    global_alignment <- .global_align_local_fits(
      locals,
      u_grid,
      z_design,
      x_design,
      initial_y,
      control,
      bandwidth = bandwidth,
      estimation_spec = estimation_spec
    )
    locals <- global_alignment$locals
    diagnostics$global_path_cost <- global_alignment$path_cost
    diagnostics$global_transition_cost <- global_alignment$transition_cost
  }

  for (grid_id in seq_along(u_grid)) {
    local <- locals[[grid_id]]
    diagnostics$ambiguous[grid_id] <- isTRUE(local$ambiguous)
    if (isTRUE(local$ambiguous)) {
      diagnostics$warnings <- c(
        diagnostics$warnings,
        sprintf("Ambiguous label alignment near u = %.4f.", u_grid[grid_id])
      )
    }

    expert[grid_id, , ] <- local$expert_intercept
    expert_slope[grid_id, , ] <- local$expert_slope
    gating[grid_id, , ] <- local$gating_intercept
    gating_slope[grid_id, , ] <- local$gating_slope
    if (!is.null(sigma)) {
      sigma[grid_id, ] <- local$sigma
    }
    if (!is.null(sigma_slope)) {
      sigma_slope[grid_id, ] <- local$sigma_slope
    }
    if (!is.null(theta)) {
      theta[grid_id, ] <- local$theta
    }
    posterior[, , grid_id] <- local$posterior

    diagnostics$loglik[grid_id] <- local$loglik
    diagnostics$converged[grid_id] <- local$converged
    diagnostics$iterations[grid_id] <- local$iterations
    diagnostics$selected_start[grid_id] <- local$selected_start
    diagnostics$permutations[grid_id, ] <- local$permutation
    diagnostics$greedy_permutations[grid_id, ] <- local$greedy_permutation %||% local$permutation
    diagnostics$global_permutations[grid_id, ] <- local$global_permutation %||% seq_len(k)
    diagnostics$expert_optimizer_convergence[grid_id, ] <- local$expert_convergence %||% rep(NA_integer_, k)
    diagnostics$expert_gradient_norm[grid_id, ] <- local$expert_gradient_norm %||% rep(NA_real_, k)
    diagnostics$expert_coef_norm[grid_id, ] <- local$expert_coef_norm %||% rep(NA_real_, k)
    diagnostics$alignment_margin[grid_id] <- local$alignment_margin
    diagnostics$posterior_entropy[grid_id] <- .weighted_entropy(
      local$posterior,
      .local_weights(u_values, u_grid[grid_id], bandwidth, estimation_spec)
    )
    diagnostics$start_loglik[[grid_id]] <- local$start_loglik
    diagnostics$objective[grid_id] <- local$objective %||% NA_real_
    diagnostics$objective_delta[grid_id] <- local$objective_delta %||% NA_real_
    diagnostics$stop_reason[grid_id] <- local$stop_reason %||% NA_character_
    diagnostics$start_objective[[grid_id]] <- local$start_objective
    diagnostics$start_converged[[grid_id]] <- local$start_converged
    diagnostics$better_unconverged_start[grid_id] <- isTRUE(local$better_unconverged_start)
  }
  # Global alignment can relabel a grid point against a confident greedy match
  # with its neighbour (an isolated swap). Flag every such link.
  if (identical(alignment_method, "global") && n_grid > 1L) {
    for (grid_id in 2:n_grid) {
      changed <- any(diagnostics$global_permutations[grid_id, ] != diagnostics$global_permutations[grid_id - 1L, ])
      diagnostics$alignment_override[grid_id] <- changed && !isTRUE(diagnostics$ambiguous[grid_id])
    }
    if (any(diagnostics$alignment_override) && isTRUE(control$warn_ambiguous)) {
      warning(
        sprintf(paste0("Global label alignment overrode a confident neighbour match at grid point(s) %s; ",
                       "compare with `label = \"greedy\"`."),
                paste(which(diagnostics$alignment_override), collapse = ", ")),
        call. = FALSE
      )
    }
  }

  if (length(diagnostics$warnings) && isTRUE(control$warn_ambiguous)) {
    warning(
      sprintf("Label alignment was ambiguous at %d grid point(s); inspect `fit$diagnostics$warnings`.",
              length(diagnostics$warnings)),
      call. = FALSE
    )
  }

  object <- list(
    call = spec$fit_call,
    formula = formula,
    family = family,
    k = k,
    label = alignment_method,
    engine_id = "local_grid_em",
    parameterization_id = parameterization,
    u_scale = u_scale,
    u_scaling = u_scaling,
    bandwidth = bandwidth,
    u_grid = u_grid,
    u_grid_original = u_grid_original,
    control = control,
    parameterization = .vcmoe_parameterization_metadata(
      family = family,
      k = k,
      bandwidth = bandwidth,
      label = label,
      alignment_method = alignment_method,
      control = control,
      u_scaling = u_scaling,
      parameterization = parameterization,
      estimation_spec = estimation_spec
    ),
    coefficients = list(
      expert = expert,
      expert_slope = expert_slope,
      gating = gating,
      gating_slope = gating_slope,
      sigma = sigma,
      sigma_slope = sigma_slope,
      theta = theta
    ),
    diagnostics = diagnostics,
    posterior = posterior,
    terms = list(
      expert = stats::terms(pieces$expert_formula),
      gating = stats::terms(pieces$gating_formula)
    ),
    response = pieces$response,
    response_info = response$info,
    u_name = u_info$name,
    rows_used = which(complete_rows),
    rows_removed = sum(!complete_rows),
    fitted = if (isTRUE(control$keep_data)) {
      list(
        y = response$observed,
        success = response$success,
        trials = response$trials,
        response = response,
        z_design = z_design,
        x_design = x_design,
        expert_offset = expert_offset,
        u = u_values,
        u_original = u_original
      )
    } else {
      NULL
    }
  )
  class(object) <- "vcmoe"
  object
}

# Ranking of local fits used when grid-parallel parts or refits compete for a
# grid point: converged before unconverged, then the higher penalized objective.
.vcmoe_local_better <- function(candidate, incumbent) {
  if (is.null(incumbent)) return(TRUE)
  value <- function(x) {
    v <- x$objective %||% x$loglik
    if (is.finite(v)) v else -Inf
  }
  cand_conv <- isTRUE(candidate$converged)
  inc_conv <- isTRUE(incumbent$converged)
  if (cand_conv != inc_conv) return(cand_conv)
  value(candidate) > value(incumbent)
}

.vcmoe_best_part_locals <- function(parts) {
  best <- list()
  source <- list()
  for (part_id in seq_along(parts)) {
    for (key in names(parts[[part_id]]$locals)) {
      local <- parts[[part_id]]$locals[[key]]
      if (.vcmoe_local_better(local, best[[key]])) {
        best[[key]] <- local
        source[[key]] <- part_id
      }
    }
  }
  list(best = best, source = source)
}

#' Merge grid-subset VCMoE fits
#'
#' Combines `vcmoe_grid_part` objects returned by [vcmoe_fit()] with
#' `control$grid_subset` (for example one per cluster job) into an ordinary
#' `vcmoe` fit. Every grid point must be covered. When several parts cover the
#' same grid point (for example a cold-start pass and a neighbour warm-start
#' pass), a converged local fit is preferred, then the one with the higher
#' penalized local objective. Label alignment is then run over the whole grid
#' as in a serial fit.
#'
#' @param parts A list of `vcmoe_grid_part` objects from the same data and
#'   settings.
#' @param label Optional label-alignment override: `"align"`, `"global"` or
#'   `"greedy"`. The default keeps the alignment of the original call.
#' @return A `vcmoe` object; `diagnostics$grid_source` records which part
#'   supplied each grid point.
#' @export
vcmoe_merge_grid_fits <- function(parts, label = NULL) {
  if (inherits(parts, "vcmoe_grid_part")) {
    parts <- list(parts)
  }
  if (!is.list(parts) || !length(parts) ||
      !all(vapply(parts, inherits, logical(1L), what = "vcmoe_grid_part"))) {
    stop("`parts` must be a list of `vcmoe_grid_part` objects.", call. = FALSE)
  }
  spec <- parts[[1L]]$spec
  # Data and settings only: formulas carry environments that differ between
  # the R processes that produced the parts.
  core <- function(s) {
    control <- s$control
    control[c("grid_seed", "grid_warm_start")] <- NULL
    c(s[c("family", "k", "alignment_method", "parameterization", "u_scale", "u_scaling", "bandwidth",
          "u_grid", "y", "trials", "initial_y", "z_design", "x_design", "expert_offset",
          "u_values", "complete_rows")],
      list(control = control, formula = paste(deparse(s$formula), collapse = " ")))
  }
  same <- vapply(parts, function(p) identical(core(p$spec), core(spec)), logical(1L))
  if (!all(same)) {
    stop("All parts must come from the same data, formula and settings (except the grid subset).", call. = FALSE)
  }
  if (!is.null(label)) {
    label <- match.arg(label, c("align", "global", "greedy"))
    spec$label <- label
    spec$alignment_method <- .vcmoe_alignment_method(label, spec$k)
  }
  n_grid <- length(spec$u_grid)
  picked <- .vcmoe_best_part_locals(parts)
  covered <- as.integer(names(picked$best))
  missing <- setdiff(seq_len(n_grid), covered)
  if (length(missing)) {
    stop(sprintf("Grid points not covered by any part: %s.", paste(missing, collapse = ", ")), call. = FALSE)
  }
  candidates <- integer(n_grid)
  for (p in parts) {
    ids <- as.integer(names(p$locals))
    candidates[ids] <- candidates[ids] + 1L
  }
  locals <- vector("list", n_grid)
  previous <- NULL
  for (grid_id in seq_len(n_grid)) {
    locals[[grid_id]] <- .align_local_fit(picked$best[[as.character(grid_id)]], previous, spec$z_design,
                                          spec$x_design, spec$initial_y, spec$control)
    previous <- locals[[grid_id]]
  }
  object <- .vcmoe_finalize_fit(spec, locals)
  object$diagnostics$grid_source <- vapply(as.character(seq_len(n_grid)),
                                           function(key) picked$source[[key]], integer(1L), USE.NAMES = FALSE)
  object$diagnostics$grid_candidates <- candidates
  object$diagnostics$grid_seed <- unique(vapply(parts, `[[`, character(1L), "grid_seed"))
  object
}

#' Neighbour warm starts for a second grid-parallel pass
#'
#' For each requested grid point, returns the best local fits across `parts`
#' (ranked as in [vcmoe_merge_grid_fits()]) at the neighbouring grid points,
#' ready to pass as `control$grid_init_locals` in a second [vcmoe_fit()] call
#' with `control$grid_subset`. This restores the serial fit's warm start in
#' both directions without the serial dependency.
#'
#' @param parts A list of `vcmoe_grid_part` objects.
#' @param grid_ids Grid indices to build starts for; default all covered ones.
#' @return A named list (by grid index) of lists of local fits.
#' @export
vcmoe_grid_neighbor_inits <- function(parts, grid_ids = NULL) {
  if (inherits(parts, "vcmoe_grid_part")) {
    parts <- list(parts)
  }
  best <- .vcmoe_best_part_locals(parts)$best
  covered <- sort(as.integer(names(best)))
  grid_ids <- as.integer(grid_ids %||% covered)
  out <- lapply(grid_ids, function(g) {
    Filter(Negate(is.null), list(best[[as.character(g - 1L)]], best[[as.character(g + 1L)]]))
  })
  names(out) <- as.character(grid_ids)
  out[vapply(out, length, integer(1L)) > 0L]
}

#' @export
print.vcmoe_grid_part <- function(x, ...) {
  cat("VCMoE grid-subset fit (merge with vcmoe_merge_grid_fits())\n")
  cat("  grid points: ", paste(x$grid_ids, collapse = ", "), " of ", length(x$spec$u_grid), "\n", sep = "")
  cat("  seed mode: ", x$grid_seed, "; warm start: ", x$grid_warm_start,
      "; neighbour pass: ", x$second_pass, "\n", sep = "")
  invisible(x)
}

# Local fit at one grid point, rebuilt from a fitted object (aligned labels).
.vcmoe_local_from_fit <- function(fit, grid_id) {
  k <- fit$k
  take <- function(arr) matrix(arr[grid_id, , ], nrow = k)
  coefs <- fit$coefficients
  list(
    expert_intercept = take(coefs$expert),
    expert_slope = take(coefs$expert_slope),
    gating_intercept = take(coefs$gating),
    gating_slope = take(coefs$gating_slope),
    sigma = if (!is.null(coefs$sigma)) as.numeric(coefs$sigma[grid_id, ]) else NULL,
    sigma_slope = if (!is.null(coefs$sigma_slope)) as.numeric(coefs$sigma_slope[grid_id, ]) else NULL,
    theta = if (!is.null(coefs$theta)) as.numeric(coefs$theta[grid_id, ]) else NULL,
    posterior = fit$posterior[, , grid_id]
  )
}

#' Refine a local-grid VCMoE fit
#'
#' Continues the local EM at selected grid points from the fitted solution,
#' with a tighter tolerance and more iterations. The relative stopping rule of
#' [vcmoe_fit()] (`tol = 1e-5` by default) corresponds to a large absolute
#' change of the local objective when the sample is large, so a fit can stop
#' before the local maximum is reached; the bread used by [vcmoe_confband()]
#' is then often not positive definite. Refining fixes this without refitting
#' from scratch. Each refined grid point keeps the labels of the original fit;
#' a refinement that would change them is reported, not applied silently.
#'
#' @param fit A local-grid `vcmoe` fit that retains its data
#'   (`control$keep_data = TRUE`).
#' @param grids Grid indices to refine. The default refines every grid point.
#' @param tol Relative tolerance on the change of the penalized local objective.
#' @param maxit Maximum number of EM iterations per grid point.
#' @return The refined `vcmoe` fit. `diagnostics$refined` marks the refined
#'   grid points, `diagnostics$refine_objective_gain` gives the change of the
#'   penalized local objective and `diagnostics$refine_relabelled` flags grid
#'   points whose refinement moved to a relabelled solution (kept with the
#'   original labels).
#' @export
vcmoe_refine <- function(fit, grids = NULL, tol = 1e-8, maxit = 1000L) {
  if (!inherits(fit, "vcmoe")) {
    stop("`fit` must be a VCMoE fit.", call. = FALSE)
  }
  if (!identical(fit$engine_id %||% "local_grid_em", "local_grid_em")) {
    stop("`vcmoe_refine()` supports local-grid fits only.", call. = FALSE)
  }
  if (is.null(fit$fitted)) {
    stop("`vcmoe_refine()` needs a fit with retained data (`control$keep_data = TRUE`).", call. = FALSE)
  }
  n_grid <- length(fit$u_grid)
  grids <- as.integer(grids %||% seq_len(n_grid))
  if (!length(grids) || anyNA(grids) || any(grids < 1L | grids > n_grid)) {
    stop(sprintf("`grids` must be grid indices in 1..%d.", n_grid), call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a positive number.", call. = FALSE)
  }
  maxit <- as.integer(maxit)
  if (length(maxit) != 1L || is.na(maxit) || maxit < 1L) {
    stop("`maxit` must be a positive integer.", call. = FALSE)
  }
  control <- .vcmoe_default_control(.vcmoe_legacy_fit_control(fit$control))
  control$tol <- tol
  control$maxit <- maxit
  fitted <- fit$fitted
  y <- fitted$response$y
  trials <- fitted$response$trials
  initial_y <- fitted$y
  estimation_spec <- .vcmoe_estimation_spec(fit$parameterization_id %||% "a1_epanechnikov_scaled")
  out <- fit
  d <- out$diagnostics
  if (is.null(d$refined)) d$refined <- logical(n_grid)
  if (is.null(d$refine_objective_gain)) d$refine_objective_gain <- rep(NA_real_, n_grid)
  if (is.null(d$refine_relabelled)) d$refine_relabelled <- logical(n_grid)
  for (field in c("objective", "objective_delta")) {
    if (is.null(d[[field]])) d[[field]] <- rep(NA_real_, n_grid)
  }
  if (is.null(d$stop_reason)) d$stop_reason <- rep(NA_character_, n_grid)
  for (grid_id in sort(unique(grids))) {
    original <- .vcmoe_local_from_fit(fit, grid_id)
    u0 <- fit$u_grid[grid_id]
    weights <- .local_weights(fitted$u, u0, fit$bandwidth, estimation_spec)
    du <- .local_du(fitted$u, u0, fit$bandwidth, estimation_spec)
    params0 <- .local_to_params(original)
    objective0 <- .local_loglik(y, trials, fit$family, fitted$z_design, fitted$x_design, du,
                                fitted$expert_offset, weights, params0, control) -
      .local_ridge_penalty(fit$family, params0, weights, original$posterior, control)
    local <- .local_em(y, trials, fit$family, initial_y, fitted$z_design, fitted$x_design, fitted$expert_offset,
                       fitted$u, u0, fit$bandwidth, fit$k, control,
                       estimation_spec = estimation_spec, init_local = original)
    aligned <- .align_local_fit(local, original, fitted$z_design, fitted$x_design, initial_y, control)
    if (any(aligned$permutation != seq_len(fit$k))) {
      d$refine_relabelled[grid_id] <- TRUE
    }
    out$coefficients$expert[grid_id, , ] <- aligned$expert_intercept
    out$coefficients$expert_slope[grid_id, , ] <- aligned$expert_slope
    out$coefficients$gating[grid_id, , ] <- aligned$gating_intercept
    out$coefficients$gating_slope[grid_id, , ] <- aligned$gating_slope
    if (!is.null(out$coefficients$sigma)) out$coefficients$sigma[grid_id, ] <- aligned$sigma
    if (!is.null(out$coefficients$sigma_slope)) out$coefficients$sigma_slope[grid_id, ] <- aligned$sigma_slope
    if (!is.null(out$coefficients$theta)) out$coefficients$theta[grid_id, ] <- aligned$theta
    out$posterior[, , grid_id] <- aligned$posterior
    d$loglik[grid_id] <- aligned$loglik
    d$converged[grid_id] <- isTRUE(aligned$converged)
    d$iterations[grid_id] <- aligned$iterations
    d$objective[grid_id] <- aligned$objective
    d$objective_delta[grid_id] <- aligned$objective_delta
    d$stop_reason[grid_id] <- aligned$stop_reason
    d$refined[grid_id] <- TRUE
    d$refine_objective_gain[grid_id] <- aligned$objective - objective0
    d$expert_optimizer_convergence[grid_id, ] <- aligned$expert_convergence %||% rep(NA_integer_, fit$k)
    d$expert_gradient_norm[grid_id, ] <- aligned$expert_gradient_norm %||% rep(NA_real_, fit$k)
    d$expert_coef_norm[grid_id, ] <- aligned$expert_coef_norm %||% rep(NA_real_, fit$k)
    d$posterior_entropy[grid_id] <- .weighted_entropy(aligned$posterior, weights)
  }
  if (any(d$refine_relabelled[grids])) {
    warning(sprintf("Refinement moved to a relabelled solution at grid point(s) %s; labels were kept.",
                    paste(which(d$refine_relabelled), collapse = ", ")), call. = FALSE)
  }
  d$refine_control <- list(tol = tol, maxit = maxit)
  out$diagnostics <- d
  out
}

#' @export
print.vcmoe <- function(x, ...) {
  cat("VCMoE fit\n")
  cat("  family: ", x$family, "\n", sep = "")
  cat("  components: ", x$k, "\n", sep = "")
  cat("  label alignment: ", x$label %||% x$diagnostics$alignment_method %||% "unknown", "\n", sep = "")
  parameterization <- vcmoe_parameterization(x)
  cat("  kernel: ", parameterization$kernel$name, " (",
      parameterization$kernel$weight_normalization, ")\n", sep = "")
  cat("  local basis: ", parameterization$local_linear_basis$slope_column, "\n", sep = "")
  cat("  u scale: ", parameterization$u_scaling$method %||% x$u_scale %||% "none", "\n", sep = "")
  cat("  grid points: ", length(x$u_grid), "\n", sep = "")
  cat("  bandwidth: ", signif(x$bandwidth, 4), "\n", sep = "")
  cat("  converged grid points: ", sum(x$diagnostics$converged), "/",
      length(x$u_grid), "\n", sep = "")
  if (length(x$diagnostics$warnings)) {
    cat("  label warnings: ", length(x$diagnostics$warnings), "\n", sep = "")
  }
  invisible(x)
}
