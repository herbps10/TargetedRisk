data <- simulate_providers(N = 5e2, providers = 5, covariates = 5, seed = 1, effect_seed = 1)
trt <- "A"
outcome <- "Y"
baseline <- paste0("W", 1:5)

test_that("direct probability weighting method works", {
  set.seed(1)
  result <- direct_pw(data, trt = trt, outcome = outcome, baseline = baseline, learners = c("mean"))

  expect_equal(unname(result$estimates[,"direct"]), c(0.618, 0.828, 0.621, 0.636, 0.758), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "direct"]), as.character(1:5))
})

test_that("direct substitution method works", {
  set.seed(1)
  result <- direct_sub(data, trt = trt, outcome = outcome, baseline = baseline, learners = c("mean"))

  expect_equal(unname(result$estimates[,"direct"]), rep(0.644, 5), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "direct"]), as.character(1:5))
})

test_that("direct one-step method works", {
  set.seed(1)
  result <- direct_onestep(data, trt = trt, outcome = outcome, baseline = baseline, learners_trt = c("mean"), learners_outcome = c("mean"))

  expect_equal(unname(result$estimates[,"direct"]), c(0.641, 0.656, 0.632, 0.638, 0.652), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "direct"]), as.character(1:5))
})

test_that("direct entropy balancing method works", {
  set.seed(1)
  result <- direct_weightit(data, trt = trt, outcome = outcome, baseline = baseline, method = "ebal")

  expect_equal(unname(result$estimates[,"direct"]), c(0.623, 0.85, 0.57, 0.628, 0.841), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "direct"]), as.character(1:5))
})

test_that("direct energy balancing method works", {
  set.seed(1)
  result <- direct_weightit(data, trt = trt, outcome = outcome, baseline = baseline, method = "energy")

  expect_equal(unname(result$estimates[,"direct"]), c(0.623, 0.85, 0.57, 0.628, 0.841), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "direct"]), as.character(1:5))
})

test_that("direct MatchIt method works", {
  set.seed(1)
  result <- direct_matchit(data, trt = trt, outcome = outcome, baseline = baseline, method = "quick", distance = "glm")

  expect_equal(unname(result$estimates[,"direct"]), c(0.623, 0.85, 0.57, 0.628, 0.841), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "direct"]), as.character(1:5))
})

test_that("direct TMLE method works", {
  set.seed(1)
  result <- direct_tmle(data, trt = trt, outcome = outcome, baseline = baseline, learners_trt = c("mean"), learners_outcome = c("mean"))

  expect_equal(unname(result$estimates[,"direct"]), c(0.641, 0.656, 0.632, 0.638, 0.652), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "direct"]), as.character(1:5))
})
