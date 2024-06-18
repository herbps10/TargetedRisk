test_that("simulation works", {
  data <- simulate_providers(N = 1e3, providers = 5, covariates = 5, seed = 1, effect_seed = 1)

  expect_true(all(paste0("W", 1:5) %in% names(data)))
  expect_true(all(paste0("g", 1:5) %in% names(data)))
  expect_true(all(paste0("Qbar", 1:5) %in% names(data)))
  expect_true(all(paste0("Y", 1:5) %in% names(data)))

  expect_equal(sort(unique(data$A)), 1:5)
  expect_equal(sort(unique(data$Y)), 0:1)
  expect_equal(sort(unique(data$provider_effect)), 0:1)

  expect_equal(nrow(data), 1e3)
})
