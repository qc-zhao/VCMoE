.match_stress_scenario <- function(scenario) {
  match.arg(
    scenario,
    c("well_separated", "moderate", "near_overlap", "crossing", "imbalanced_gating")
  )
}

.truth_functions_gaussian_generic <- function(u, k, separation, scenario) {
  n <- length(u)
  component_center <- seq(-(k - 1) / 2, (k - 1) / 2, length.out = k)
  phase <- seq_len(k) * 0.65
  effective_sep <- separation
  if (identical(scenario, "moderate")) {
    effective_sep <- 0.55 * separation
  }
  if (identical(scenario, "near_overlap")) {
    effective_sep <- min(separation, 0.12)
  }

  expert <- array(NA_real_, dim = c(n, k, 2L),
                  dimnames = list(NULL, paste0("component", seq_len(k)),
                                  c("(Intercept)", "z1")))
  gating <- array(NA_real_, dim = c(n, k, 2L),
                  dimnames = list(NULL, paste0("component", seq_len(k)),
                                  c("(Intercept)", "x1")))
  sigma <- matrix(NA_real_, nrow = n, ncol = k,
                  dimnames = list(NULL, paste0("component", seq_len(k))))

  base_intercept <- 0.08 * sin(pi * u)
  base_slope <- 0.12 + 0.08 * u
  for (component in seq_len(k)) {
    center <- component_center[component]
    if (identical(scenario, "crossing")) {
      intercept_pattern <- center * (2 * u - 1) +
        0.10 * sin(2 * pi * u + phase[component])
      slope_pattern <- -0.16 * center +
        0.08 * cos(2 * pi * u + phase[component] / 2)
    } else {
      intercept_pattern <- 0.55 * center +
        0.14 * sin(2 * pi * u + phase[component])
      slope_pattern <- -0.18 * center +
        0.08 * cos(2 * pi * u + phase[component] / 2)
    }
    expert[, component, 1L] <- base_intercept + effective_sep * intercept_pattern
    expert[, component, 2L] <- base_slope + effective_sep * slope_pattern
    sigma[, component] <- 0.28 + 0.025 * component + 0.015 * sin(2 * pi * u + phase[component])

    gating[, component, 1L] <- 0.18 * center + 0.10 * sin(pi * u + phase[component])
    gating[, component, 2L] <- 0.12 * cos(pi * u + phase[component] / 3)
  }

  if (identical(scenario, "imbalanced_gating")) {
    gating[, , 1L] <- sweep(
      gating[, , 1L, drop = FALSE][, , 1L],
      2L,
      seq(0.90, -0.90, length.out = k),
      FUN = "+"
    )
  }

  for (idx in seq_along(u)) {
    gating[idx, , ] <- .center_logits(gating[idx, , , drop = FALSE][1L, , ])
  }

  list(expert = expert, gating = gating, sigma = sigma)
}

.truth_functions_gaussian <- function(u, k, separation, scenario = "well_separated") {
  u <- as.numeric(u)
  scenario <- .match_stress_scenario(scenario)
  n <- length(u)
  if (k > 3L || !identical(scenario, "well_separated")) {
    return(.truth_functions_gaussian_generic(u, k, separation, scenario))
  }

  expert <- array(NA_real_, dim = c(n, k, 2L),
                  dimnames = list(NULL, paste0("component", seq_len(k)),
                                  c("(Intercept)", "z1")))
  gating <- array(NA_real_, dim = c(n, k, 2L),
                  dimnames = list(NULL, paste0("component", seq_len(k)),
                                  c("(Intercept)", "x1")))
  sigma <- matrix(NA_real_, nrow = n, ncol = k,
                  dimnames = list(NULL, paste0("component", seq_len(k))))

  base_intercept <- 0.10 * sin(pi * u)
  base_slope <- 0.15 + 0.10 * u
  expert[, 1L, 1L] <- base_intercept + separation * (-0.45 + 0.35 * sin(2 * pi * u))
  expert[, 1L, 2L] <- base_slope + separation * (0.35 + 0.25 * u)
  expert[, 2L, 1L] <- base_intercept + separation * (0.55 + 0.25 * cos(2 * pi * u))
  expert[, 2L, 2L] <- base_slope + separation * (-0.35 + 0.30 * (u - 0.5)^2)
  sigma[, 1L] <- 0.35 + 0.05 * u
  sigma[, 2L] <- 0.45 - 0.03 * u

  eta <- -0.2 + 0.7 * (u - 0.5)
  eta_x <- 0.45 - 0.55 * u
  gating[, 1L, 1L] <- eta / 2
  gating[, 2L, 1L] <- -eta / 2
  gating[, 1L, 2L] <- eta_x / 2
  gating[, 2L, 2L] <- -eta_x / 2

  if (k >= 3L) {
    expert[, 3L, 1L] <- 0.15 * separation - 0.30 * sin(pi * u)
    expert[, 3L, 2L] <- 0.10 + 0.45 * cos(pi * u)
    sigma[, 3L] <- 0.40 + 0.02 * sin(2 * pi * u)
    gating[, 3L, 1L] <- 0.20 - 0.25 * u
    gating[, 3L, 2L] <- -0.15 + 0.35 * u

    for (idx in seq_along(u)) {
      gating[idx, , ] <- .center_logits(gating[idx, , , drop = FALSE][1L, , ])
    }
  }

  list(expert = expert, gating = gating, sigma = sigma)
}

.truth_functions_binomial <- function(u, k, separation, scenario = "well_separated") {
  truth <- .truth_functions_gaussian(u, k, separation, scenario)
  truth$expert <- 0.85 * truth$expert
  truth$sigma <- NULL
  truth
}

.truth_functions_negbin <- function(u, k, separation, scenario = "well_separated",
                                    mean_count = 5) {
  truth <- .truth_functions_gaussian(u, k, separation, scenario)
  truth$expert <- 0.35 * truth$expert
  truth$expert[, , 1L] <- truth$expert[, , 1L] + log(mean_count)
  theta <- matrix(NA_real_, nrow = length(u), ncol = k,
                  dimnames = list(NULL, paste0("component", seq_len(k))))
  phase <- seq_len(k) * 0.41
  for (component in seq_len(k)) {
    theta[, component] <- 7 + 0.55 * component + 0.8 * sin(2 * pi * u + phase[component])
  }
  truth$theta <- pmax(theta, 0.5)
  truth$sigma <- NULL
  truth
}

.validate_simulation_trials <- function(trials, n) {
  if (length(trials) == 1L) {
    trials <- rep(trials, n)
  }
  if (length(trials) != n) {
    stop("`trials` must be length 1 or length `n`.", call. = FALSE)
  }
  trials <- as.numeric(trials)
  if (any(!.is_whole_count(trials)) || any(trials <= 0)) {
    stop("`trials` must contain positive finite whole numbers.", call. = FALSE)
  }
  trials
}

#' Simulate Gaussian VCMoE data
#'
#' @param n Number of observations.
#' @param k Number of components. Values 2 through 10 are supported.
#' @param seed Optional random seed.
#' @param separation Controls expert separation. Values near zero create a
#'   deliberately weakly identified scenario.
#' @param u Optional numeric vector of index values.
#' @param scenario Simulation scenario: `"well_separated"`, `"moderate"`,
#'   `"near_overlap"`, `"crossing"`, or `"imbalanced_gating"`.
#' @return A list with `data` and `truth`. The `truth` entry includes
#'   coefficients, gating logits, probabilities, means, standard deviations,
#'   and sampled class labels.
#' @export
simulate_vcmoe_gaussian <- function(n = 300L, k = 2L, seed = NULL,
                                    separation = 1, u = NULL,
                                    scenario = "well_separated") {
  k <- as.integer(k)
  scenario <- .match_stress_scenario(scenario)
  if (!k %in% 2L:10L) {
    stop("`simulate_vcmoe_gaussian()` currently supports k = 2 through 10.", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(u)) {
    u <- sort(stats::runif(n))
  } else {
    if (length(u) != n) {
      stop("`u` must have length `n`.", call. = FALSE)
    }
    u <- as.numeric(u)
  }

  z1 <- stats::rnorm(n)
  x1 <- stats::rnorm(n)
  truth <- .truth_functions_gaussian(u, k, separation, scenario)
  z_design <- cbind("(Intercept)" = 1, z1 = z1)
  x_design <- cbind("(Intercept)" = 1, x1 = x1)

  logits <- matrix(0, nrow = n, ncol = k)
  means <- matrix(0, nrow = n, ncol = k)
  for (component in seq_len(k)) {
    logits[, component] <- rowSums(x_design * truth$gating[, component, ])
    means[, component] <- rowSums(z_design * truth$expert[, component, ])
  }
  probs <- .row_softmax(logits)
  component <- vapply(seq_len(n), function(i) {
    sample.int(k, size = 1L, prob = probs[i, ])
  }, integer(1L))
  y <- stats::rnorm(n, mean = means[cbind(seq_len(n), component)],
                    sd = truth$sigma[cbind(seq_len(n), component)])

  data <- data.frame(
    y = y,
    u = u,
    z1 = z1,
    x1 = x1,
    component = component
  )

  truth$probability <- probs
  truth$logits <- logits
  truth$mean <- means
  truth$component <- component
  truth$scenario <- scenario
  truth$terms <- list(expert = colnames(z_design), gating = colnames(x_design))
  list(data = data, truth = truth)
}

#' Simulate Binomial VCMoE data
#'
#' @param n Number of observations.
#' @param k Number of components. Values 2 through 10 are supported.
#' @param seed Optional random seed.
#' @param separation Controls expert separation. Values near zero create a
#'   deliberately weakly identified scenario.
#' @param u Optional numeric vector of index values.
#' @param scenario Simulation scenario: `"well_separated"`, `"moderate"`,
#'   `"near_overlap"`, `"crossing"`, or `"imbalanced_gating"`.
#' @param trials Binomial trial counts. Use `1` for Bernoulli data, or a
#'   positive integer scalar/vector for grouped binomial data.
#' @return A list with `data` and `truth`. Expert truth is on the logit scale.
#' @export
simulate_vcmoe_binomial <- function(n = 300L, k = 2L, seed = NULL,
                                    separation = 1, u = NULL,
                                    scenario = "well_separated",
                                    trials = 1L) {
  k <- as.integer(k)
  scenario <- .match_stress_scenario(scenario)
  if (!k %in% 2L:10L) {
    stop("`simulate_vcmoe_binomial()` currently supports k = 2 through 10.", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(u)) {
    u <- sort(stats::runif(n))
  } else {
    if (length(u) != n) {
      stop("`u` must have length `n`.", call. = FALSE)
    }
    u <- as.numeric(u)
  }
  trials <- .validate_simulation_trials(trials, n)

  z1 <- stats::rnorm(n)
  x1 <- stats::rnorm(n)
  truth <- .truth_functions_binomial(u, k, separation, scenario)
  z_design <- cbind("(Intercept)" = 1, z1 = z1)
  x_design <- cbind("(Intercept)" = 1, x1 = x1)

  logits <- matrix(0, nrow = n, ncol = k)
  success_probability <- matrix(0, nrow = n, ncol = k)
  for (component in seq_len(k)) {
    logits[, component] <- rowSums(x_design * truth$gating[, component, ])
    success_probability[, component] <- stats::plogis(rowSums(z_design * truth$expert[, component, ]))
  }
  probs <- .row_softmax(logits)
  component <- vapply(seq_len(n), function(i) {
    sample.int(k, size = 1L, prob = probs[i, ])
  }, integer(1L))
  success <- stats::rbinom(
    n,
    size = trials,
    prob = success_probability[cbind(seq_len(n), component)]
  )
  failure <- trials - success

  if (all(trials == 1)) {
    data <- data.frame(
      y = success,
      u = u,
      z1 = z1,
      x1 = x1,
      component = component,
      success = success,
      failure = failure,
      trials = trials
    )
  } else {
    data <- data.frame(
      success = success,
      failure = failure,
      trials = trials,
      u = u,
      z1 = z1,
      x1 = x1,
      component = component
    )
  }

  truth$probability <- probs
  truth$logits <- logits
  truth$success_probability <- success_probability
  truth$component <- component
  truth$success <- success
  truth$failure <- failure
  truth$trials <- trials
  truth$scenario <- scenario
  truth$terms <- list(expert = colnames(z_design), gating = colnames(x_design))
  list(data = data, truth = truth)
}

#' Simulate Negative-Binomial VCMoE count data
#'
#' @param n Number of observations.
#' @param k Number of components. Values 2 through 10 are supported.
#' @param seed Optional random seed.
#' @param separation Controls expert separation.
#' @param u Optional numeric vector of index values.
#' @param scenario Simulation scenario: `"well_separated"`, `"moderate"`,
#'   `"near_overlap"`, `"crossing"`, or `"imbalanced_gating"`.
#' @param size_factor Optional positive size factors. If `NULL`, log-normal
#'   size factors are generated.
#' @param mean_count Baseline count scale.
#' @return A list with `data` and `truth`. Expert truth is on the log mean scale.
#' @export
simulate_vcmoe_negbin <- function(n = 300L, k = 2L, seed = NULL,
                                  separation = 1, u = NULL,
                                  scenario = "well_separated",
                                  size_factor = NULL,
                                  mean_count = 5) {
  k <- as.integer(k)
  scenario <- .match_stress_scenario(scenario)
  if (!k %in% 2L:10L) {
    stop("`simulate_vcmoe_negbin()` currently supports k = 2 through 10.", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(u)) {
    u <- sort(stats::runif(n))
  } else {
    if (length(u) != n) {
      stop("`u` must have length `n`.", call. = FALSE)
    }
    u <- as.numeric(u)
  }
  if (is.null(size_factor)) {
    size_factor <- exp(stats::rnorm(n, sd = 0.25))
  } else {
    if (length(size_factor) == 1L) {
      size_factor <- rep(size_factor, n)
    }
    if (length(size_factor) != n || any(!is.finite(size_factor)) || any(size_factor <= 0)) {
      stop("`size_factor` must contain positive finite values.", call. = FALSE)
    }
    size_factor <- as.numeric(size_factor)
  }
  log_size_factor <- log(size_factor)

  z1 <- stats::rnorm(n)
  x1 <- stats::rnorm(n)
  truth <- .truth_functions_negbin(u, k, separation, scenario, mean_count = mean_count)
  z_design <- cbind("(Intercept)" = 1, z1 = z1)
  x_design <- cbind("(Intercept)" = 1, x1 = x1)

  logits <- matrix(0, nrow = n, ncol = k)
  means <- matrix(0, nrow = n, ncol = k)
  for (component in seq_len(k)) {
    logits[, component] <- rowSums(x_design * truth$gating[, component, ])
    means[, component] <- exp(pmin(pmax(log_size_factor +
      rowSums(z_design * truth$expert[, component, ]), -20), 20))
  }
  probs <- .row_softmax(logits)
  component <- vapply(seq_len(n), function(i) {
    sample.int(k, size = 1L, prob = probs[i, ])
  }, integer(1L))
  y <- stats::rnbinom(
    n,
    size = truth$theta[cbind(seq_len(n), component)],
    mu = means[cbind(seq_len(n), component)]
  )

  data <- data.frame(
    y = y,
    u = u,
    z1 = z1,
    x1 = x1,
    size_factor = size_factor,
    log_size_factor = log_size_factor,
    component = component
  )

  truth$probability <- probs
  truth$logits <- logits
  truth$mean <- means
  truth$component <- component
  truth$offset <- log_size_factor
  truth$scenario <- scenario
  truth$terms <- list(expert = colnames(z_design), gating = colnames(x_design))
  list(data = data, truth = truth)
}
