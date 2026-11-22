refine_grid_input <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      sim <- simulate_vcmoe_negbin(n = 300, k = 2, seed = 30003, separation = 1.5)
      part <- vcmoe_fit(
        y ~ z1 + offset(log_size_factor) | x1,
        data = sim$data, u = "u", family = "negative-binomial",
        bandwidth = 0.3, u_grid = seq(0.1, 0.9, length.out = 7),
        control = list(maxit = 200, n_starts = 3, seed = 31000, warn_ambiguous = FALSE,
                       grid_seed = "per_grid", grid_warm_start = FALSE, grid_subset = 1:7)
      )
      cache <<- list(sim = sim, part = part)
    }
    cache
  }
})

refine_case <- function(label = "align", warn = FALSE) {
  input <- refine_grid_input()
  part <- input$part
  part$spec$control$warn_ambiguous <- warn
  messages <- character(0L)
  fit <- withCallingHandlers(vcmoe_merge_grid_fits(part, label = label), warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(sim = input$sim, fit = fit, warnings = messages)
}

test_that("vcmoe_refine() raises the penalized objective and keeps the labels", {
  skip_on_cran()
  case <- refine_case()
  fit <- case$fit
  refined <- vcmoe_refine(fit, tol = 1e-9, maxit = 2000)
  d <- refined$diagnostics
  expect_true(all(d$refined))
  expect_true(all(d$refine_objective_gain >= -1e-8))
  expect_true(any(d$refine_objective_gain > 1e-3))
  expect_true(all(d$objective >= fit$diagnostics$objective - 1e-8))
  expect_false(any(d$refine_relabelled))
  expect_identical(d$stop_reason, rep("tolerance", 7L))
  # Labels follow the original fit: refined coefficients stay close.
  expect_lt(max(abs(refined$coefficients$expert - fit$coefficients$expert)), 1)
  expect_true(all(vcmoe_diagnostics(refined)$refined))
  partial <- vcmoe_refine(fit, grids = 2L, tol = 1e-9)
  expect_identical(partial$diagnostics$refined, seq_len(7) == 2L)
  expect_identical(partial$coefficients$expert[-2L, , , drop = FALSE], fit$coefficients$expert[-2L, , , drop = FALSE])
  # A refined grid point is the fixed point of a further refinement.
  again <- vcmoe_refine(refined, tol = 1e-9, maxit = 2000)
  expect_lt(max(abs(again$diagnostics$refine_objective_gain)), 1e-3)
})

test_that("vcmoe_refine() validates its inputs", {
  case <- refine_case()
  expect_error(vcmoe_refine(case$fit, grids = 9L), "grid indices")
  expect_error(vcmoe_refine(case$fit, tol = 0), "tol")
  expect_error(vcmoe_refine(case$fit, maxit = 0), "maxit")
  no_data <- case$fit
  no_data$fitted <- NULL
  expect_error(vcmoe_refine(no_data), "retained data")
  expect_error(vcmoe_refine(list()), "VCMoE fit")
})

test_that("global alignment that overrides a confident neighbour match is flagged", {
  case <- refine_case()
  fit <- case$fit
  expect_identical(which(fit$diagnostics$alignment_override), 7L)
  expect_true(any(grepl("overrode a confident neighbour match at grid point\\(s\\) 7",
                         refine_case(warn = TRUE)$warnings)))
  expect_false(any(grepl("overrode", case$warnings)))
  greedy <- refine_case(label = "greedy")$fit
  expect_false(any(greedy$diagnostics$alignment_override))
  # The greedy labels agree with the simulated components at the flagged grid;
  # the global labels are swapped there.
  window <- abs(case$sim$data$u - fit$u_grid[[7L]]) < 0.15
  truth <- case$sim$data$component[window] == 1L
  expect_gt(stats::cor(greedy$posterior[window, 1L, 7L], truth), 0)
  expect_lt(stats::cor(fit$posterior[window, 1L, 7L], truth), 0)
})

test_that("Negative-Binomial eta clipping ignores observations outside the kernel support", {
  sim <- simulate_vcmoe_negbin(n = 80, k = 2, seed = 405, separation = 1.25)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 + offset(log_size_factor) | x1,
    data = sim$data, u = "u", family = "negative-binomial",
    bandwidth = 0.2, u_grid = c(0, 0.5, 1),
    control = list(maxit = 8, n_starts = 1, seed = 406, warn_ambiguous = FALSE)
  ))
  control <- VCMoE:::.vcmoe_dev_control(list())
  clipping <- function(object) {
    c(count = vcmoe_diagnostics(object)$negbin_eta_clipping_count[[3L]],
      blocked = grepl("negbin_eta_clipping", VCMoE:::.vcmoe_dev_grid_block_reason(object, 3L, control)))
  }
  # At the u = 1 boundary a steep intercept slope passes the eta bound only
  # where (u - 1) / h < -1, outside the kernel support.
  far <- fit
  far$coefficients$expert_slope[3L, 1L, "(Intercept)"] <- -5
  du <- (sim$data$u - 1) / fit$bandwidth
  expect_true(any(du < -4.5))
  expect_equal(clipping(far), c(count = 0, blocked = 0))
  inside <- far
  inside$coefficients$expert[3L, 1L, "(Intercept)"] <- 25
  expect_true(clipping(inside)[["count"]] > 0)
  expect_equal(clipping(inside)[["blocked"]], 1)
})
