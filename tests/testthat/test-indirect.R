data <- simulate_providers(N = 5e2, providers = 5, covariates = 5, seed = 1, effect_seed = 1)
trt <- "A"
outcome <- "Y"
baseline <- paste0("W", 1:5)

test_that("indirect probability weighting method works", {
  set.seed(1)
  result <- indirect_pw(data, trt = trt, outcome = outcome, baseline = baseline, learners = c("mean"))

  expect_equal(unname(result$estimates[,"SMR"]), c(0.949, 1.261, 0.952, 0.98, 1.168), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "SMR"]), as.character(1:5))
})

test_that("indirect substitution method works", {
  set.seed(1)
  result <- indirect_sub(data, trt = trt, outcome = outcome, baseline = baseline, learners = c("mean"))

  expect_equal(unname(result$estimates[,"SMR"]), c(0.949, 1.249, 0.95, 0.976, 1.155), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "SMR"]), as.character(1:5))
})

test_that("indirect one-step method works", {
  set.seed(1)
  result <- indirect_onestep(data, trt = trt, outcome = outcome, baseline = baseline, learners_trt = c("mean"), learners_outcome = c("mean"))

  expect_equal(unname(result$estimates[,"SMR"]), c(0.96, 1.24, 0.95, 0.97, 1.13), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "SMR"]), as.character(1:5))
})

test_that("indirect TMLE method works", {
  set.seed(1)
  result <- indirect_tmle(data, trt = trt, outcome = outcome, baseline = baseline, learners_trt = c("mean"), learners_outcome = c("mean"))

  expect_equal(unname(result$estimates[,"SMR"]), c(0.946, 1.239, 0.957, 0.977, 1.165), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "SMR"]), as.character(1:5))
})

test_that("indirect MatchIt method works", {
  set.seed(1)
  result <- indirect_matchit(data, trt = trt, outcome = outcome, baseline = baseline, method = "nearest", distance = "glm")

  expect_equal(unname(result$estimates[,"SMR"]), c(1, 1.4, 0.88, 0.91, 1.24), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "SMR"]), as.character(1:5))
})

test_that("indirect entropy balancing method works", {
  set.seed(1)
  result <- indirect_weightit(data, trt = trt, outcome = outcome, baseline = baseline, method = "ebal")

  expect_equal(unname(result$estimates[,"SMR"]), c(1, 1.4, 0.88, 0.91, 1.24), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "SMR"]), as.character(1:5))
})

test_that("indirect energy balancing method works", {
  set.seed(1)
  result <- indirect_weightit(data, trt = trt, outcome = outcome, baseline = baseline, method = "energy")

  expect_equal(unname(result$estimates[,"SMR"]), c(1, 1.4, 0.88, 0.91, 1.24), tolerance = 1e-1)
  expect_equal(names(result$estimates[, "SMR"]), as.character(1:5))
})
