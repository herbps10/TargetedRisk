data <- simulate_providers(N = 5e2, providers = 5, covariates = 5, seed = 1, effect_seed = 1)
trt <- "A"
outcome <- "Y"
baseline <- paste0("W", 1:5)
.torch_params = list(
  hidden = 20,
  hidden2 = 20,
  learning_rate = 1e-2,
  epochs = 25,
  dropout = 0.05,
  seed = 1
)

test_that("direct probability weighting method works", {
  set.seed(1)
  result <- direct_pw(data, trt = trt, outcome = outcome, baseline = baseline, trt_method = "default", learners = c("mean"))

  expect_equal(unname(result$estimates[,"direct"]), c(0.618, 0.828, 0.621, 0.636, 0.758), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "direct"]), as.character(1:5))
})

test_that("direct probability weighting method with SuperRiesz works", {
  set.seed(1)
  result <- direct_pw(data, trt = trt, outcome = outcome, baseline = baseline, trt_method = "SuperRiesz", learners = list(list("torch", epochs = 25)), control = standardization_control(.learners_trt_folds = 1))

  expect_equal(unname(result$estimates[,"direct"]), c(0.037, 0.044, 0.171, 0.186, 0.045), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "direct"]), as.character(1:5))
})

test_that("direct probability weighting method with Torch works", {
  set.seed(1)
  result <- direct_pw(data, trt = trt, outcome = outcome, baseline = baseline, trt_method = "torch", control = standardization_control(.learners_trt_folds = 1, .torch_params = .torch_params))

  expect_equal(unname(result$estimates[,"direct"]), c(0.1, 0.135, 0.411, 0.409, 0.139), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "direct"]), as.character(1:5))
})

test_that("direct substitution method works", {
  set.seed(1)
  result <- direct_sub(data, trt = trt, outcome = outcome, baseline = baseline, learners = c("mean"))

  expect_equal(unname(result$estimates[,"direct"]), rep(0.644, 5), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "direct"]), as.character(1:5))
})


test_that("direct TMLE method works", {
  set.seed(1)
  result <- direct_tmle(data, trt = trt, outcome = outcome, baseline = baseline, trt_method = "default", learners_trt = c("mean"), learners_outcome = c("mean"))

  expect_equal(unname(result$estimates[,"direct"]), c(0.641, 0.656, 0.632, 0.638, 0.652), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "direct"]), as.character(1:5))
})

test_that("direct TMLE method with SuperRiesz works", {
  set.seed(1)
  result <- direct_tmle(data, trt = trt, outcome = outcome, baseline = baseline, trt_method = "SuperRiesz", learners_trt = list(list("torch", epochs = 25)), learners_outcome = c("mean"), control = standardization_control(.learners_trt_folds = 1))

  expect_equal(unname(result$estimates[,"direct"]), c(0.641, 0.656, 0.62, 0.637, 0.656), tolerance = 1e-2)
  expect_equal(names(result$estimates[, "direct"]), as.character(1:5))
})

test_that("direct TMLE method with Torch works", {
  set.seed(1)
  result <- direct_tmle(data, trt = trt, outcome = outcome, baseline = baseline, trt_method = "torch", control = standardization_control(.learners_trt_folds = 1, .torch_params = .torch_params))

  expect_equal(unname(result$estimates[,"direct"]), c(0.645, 0.792, 0.598, 0.627, 0.765), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "direct"]), as.character(1:5))
})
