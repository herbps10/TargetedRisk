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

test_that("indirect probability weighting method works", {
  set.seed(1)
  result <- indirect_pw(data, trt = trt, outcome = outcome, baseline = baseline, trt_method = "default", learners = c("mean"))

  expect_equal(unname(result$estimates[,"SMR"]), c(0.949, 1.261, 0.952, 0.98, 1.168), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "SMR"]), as.character(1:5))
})

test_that("indirect probability weighting method with SuperRiesz works", {
  set.seed(1)
  result <- indirect_pw(data, trt = trt, outcome = outcome, baseline = baseline, trt_method = "SuperRiesz", learners = list(list("nn", epochs = 25)), control = standardization_control(.learners_trt_folds = 1))

  expect_equal(unname(result$estimates[,"SMR"]), c(0.116, 0.135, 0.524, 0.573, 0.139), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "SMR"]), as.character(1:5))
})

test_that("indirect probability weighting method with Torch works", {
  set.seed(1)
  result <- indirect_pw(data, trt = trt, outcome = outcome, baseline = baseline, trt_method = "torch", control = standardization_control(.learners_trt_folds = 1, .torch_params = .torch_params))

  expect_equal(unname(result$estimates[,"SMR"]), c(19.49, 11.20, 0.88, 1.05, 8.98), tolerance = 1e-2)
  expect_equal(names(result$estimates[, "SMR"]), as.character(1:5))
})

test_that("indirect substitution method works", {
  set.seed(1)
  result <- indirect_sub(data, trt = trt, outcome = outcome, baseline = baseline, learners = c("mean"))

  expect_equal(unname(result$estimates[,"SMR"]), c(0.949, 1.249, 0.95, 0.976, 1.155), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "SMR"]), as.character(1:5))
})


test_that("indirect TMLE method works", {
  set.seed(1)
  result <- indirect_tmle(data, trt = trt, outcome = outcome, baseline = baseline, trt_method = "default", learners_trt = c("mean"), learners_outcome = c("mean"))

  expect_equal(unname(result$estimates[,"SMR"]), c(0.946, 1.239, 0.957, 0.977, 1.165), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "SMR"]), as.character(1:5))
})

test_that("indirect TMLE method with SuperRiesz works", {
  set.seed(1)
  result <- indirect_tmle(data, trt = trt, outcome = outcome, baseline = baseline, trt_method = "SuperRiesz", learners_trt = list(list("nn", epochs = 25)), learners_outcome = c("mean"), control = standardization_control(.learners_trt_folds = 1))

  expect_equal(unname(result$estimates[,"SMR"]), c(0.939, 1.269, 0.947, 0.976, 1.146), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "SMR"]), as.character(1:5))
})

test_that("indirect TMLE method with Torch works", {
  set.seed(1)
  result <- indirect_tmle(data, trt = trt, outcome = outcome, baseline = baseline, trt_method = "torch", control = standardization_control(.learners_trt_folds = 1, .torch_params = .torch_params))

  expect_equal(unname(result$estimates[,"SMR"]), c(0.976, 1.336, 0.919, 0.969, 1.272), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "SMR"]), as.character(1:5))
})
