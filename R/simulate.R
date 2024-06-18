#' Simulate test data for provider profiling
#'
#' @param N number of observations
#' @param providers number of providers
#' @param covariates number of baseline covariates
#' @param seed random seed for generating provider effects
#' @param effect_seed random seed for generating observations
#'
#' @importFrom purrr map_int
#' @export
simulate_providers <- function(N = 1e3, providers = 5, covariates = 5, seed = NA, effect_seed = NA) {
  if(!is.na(seed)) set.seed(effect_seed)
  provider_effects <- rbinom(providers, 1, 0.5)

  if(!is.na(effect_seed)) set.seed(seed)
  W <- matrix(runif(N * covariates, 0, 1), ncol = covariates, nrow = N)

  colnames(W) <- paste0("W", 1:covariates)
  g <- matrix(0, ncol = providers, nrow = N)
  g <- 1 + 10 * (matrix(W[, 1], ncol = providers, nrow = N, byrow = FALSE)) * matrix(provider_effects, ncol = providers, nrow = N, byrow = TRUE)
  g <- g / rowSums(g)
  colnames(g) <- paste0("g", 1:providers)

  A <- map_int(1:N, \(i) sample(1:providers, 1, prob = g[i, ]))

  trt_indicator <- matrix(0, ncol = providers, nrow = N)
  for(a in 1:providers) {
    trt_indicator[, a] <- A == a
  }

  Qbar <- plogis(2 * matrix(W[,1], ncol = providers, nrow = N, byrow = FALSE) + matrix(W[,5], ncol = providers, nrow = N, byrow = FALSE) - matrix(provider_effects, ncol = providers, nrow = N, byrow = TRUE))
  colnames(Qbar) <- paste0("Qbar", 1:providers)

  Qtilde <- rowSums(Qbar * g)

  Ya <- matrix(rbinom(N * providers, 1, as.vector(Qbar)), ncol = providers, nrow = N, byrow = FALSE)
  colnames(Ya) <- paste0("Y", 1:providers)

  Y <- rowSums(Ya * trt_indicator)
  YA <- rowSums(Ya * g)

  data <- cbind(
    as.data.frame(W),
    as.data.frame(g),
    as.data.frame(Qbar),
    as.data.frame(Ya)
  )

  data$Y <- Y
  data$A <- A
  data$YA <- YA
  data$Qtilde <- Qtilde
  data$provider_effect <- provider_effects[A]
  data
}
