test_that("Gaussian simulator returns expected data and truth", {
  sim <- simulate_vcmoe_gaussian(n = 40, k = 2, seed = 1)

  expect_s3_class(sim$data, "data.frame")
  expect_equal(nrow(sim$data), 40)
  expect_true(all(c("y", "u", "z1", "x1", "component") %in% names(sim$data)))
  expect_equal(dim(sim$truth$expert), c(40, 2, 2))
  expect_equal(dim(sim$truth$gating), c(40, 2, 2))
  expect_equal(dim(sim$truth$sigma), c(40, 2))
  expect_true(all(sim$data$component %in% c(1L, 2L)))
})

test_that("Gaussian simulator supports multi-class stress scenarios", {
  sim <- simulate_vcmoe_gaussian(
    n = 50,
    k = 5,
    seed = 11,
    separation = 1.1,
    scenario = "crossing"
  )

  expect_equal(nrow(sim$data), 50)
  expect_equal(dim(sim$truth$expert), c(50, 5, 2))
  expect_equal(dim(sim$truth$gating), c(50, 5, 2))
  expect_equal(dim(sim$truth$sigma), c(50, 5))
  expect_equal(dim(sim$truth$logits), c(50, 5))
  expect_equal(dim(sim$truth$probability), c(50, 5))
  expect_equal(as.numeric(rowSums(sim$truth$probability)), rep(1, 50), tolerance = 1e-8)
  expect_equal(sim$truth$scenario, "crossing")
  expect_true(all(sim$data$component %in% seq_len(5)))
})

test_that("Binomial simulator returns Bernoulli data and truth", {
  sim <- simulate_vcmoe_binomial(n = 40, k = 2, seed = 13)

  expect_s3_class(sim$data, "data.frame")
  expect_equal(nrow(sim$data), 40)
  expect_true(all(c("y", "u", "z1", "x1", "component", "success", "failure", "trials") %in% names(sim$data)))
  expect_true(all(sim$data$y %in% c(0, 1)))
  expect_equal(dim(sim$truth$expert), c(40, 2, 2))
  expect_equal(dim(sim$truth$gating), c(40, 2, 2))
  expect_equal(dim(sim$truth$success_probability), c(40, 2))
  expect_equal(dim(sim$truth$probability), c(40, 2))
  expect_true(all(sim$truth$success_probability >= 0 & sim$truth$success_probability <= 1))
  expect_true(all(sim$data$component %in% c(1L, 2L)))
})

test_that("Binomial simulator returns grouped data", {
  sim <- simulate_vcmoe_binomial(
    n = 35,
    k = 3,
    seed = 14,
    trials = 5,
    scenario = "crossing"
  )

  expect_equal(nrow(sim$data), 35)
  expect_false("y" %in% names(sim$data))
  expect_true(all(c("success", "failure", "trials", "u", "z1", "x1", "component") %in% names(sim$data)))
  expect_equal(sim$data$success + sim$data$failure, sim$data$trials)
  expect_equal(dim(sim$truth$expert), c(35, 3, 2))
  expect_equal(dim(sim$truth$gating), c(35, 3, 2))
  expect_equal(dim(sim$truth$success_probability), c(35, 3))
  expect_equal(sim$truth$scenario, "crossing")
})
