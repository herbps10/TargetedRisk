eif_direct <- function(a, y, trt_levels, g, Qtilde, theta) {
  eif <- matrix(nrow = length(a), ncol = length(trt_levels))
  colnames(eif) <- trt_levels

  for(trt_level in trt_levels) {
    eif[, trt_level] <- (a == trt_level) / g$treatment_probs[, trt_level] * (y - Qtilde[, trt_level]) + Qtilde[, trt_level] - theta[trt_level]
  }
  eif
}

#' @importFrom stats pnorm qnorm sd
theta_direct_onestep <- function(task, Qtilde, g) {
  theta <- colMeans(Qtilde)
  theta_eif <- eif_direct(task$data[[task$trt]], task$data[[task$outcome]], task$trt_levels, g, Qtilde, theta)

  # Apply one-step correction
  theta <- theta + colMeans(theta_eif)
  theta_eif <- eif_direct(task$data[[task$trt]], task$data[[task$outcome]], task$trt_levels, g, Qtilde, theta)

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("direct")
  colnames(p_values) <- "direct"
  estimates[, 1] <- theta

  se[, 1] <- apply(theta_eif, 2, sd) / sqrt(task$n)

  alpha <- 0.95
  low  <- estimates + qnorm((1 - alpha) / 2) * se
  high <- estimates + qnorm(1 - (1 - alpha) / 2) * se

  p_values <- 2 * pnorm(abs((estimates[, 1] - 1) / se[, 1]), lower.tail = FALSE)

  result <- list(
    estimator = "Onestep",
    estimates = estimates,
    se = se,
    alpha = 0.95,
    low = low,
    high = high,
    p_values = p_values
  )
  class(result) <- "direct"
  result
}

#' @importFrom marginaleffects avg_predictions
theta_direct_weightit <- function(task, weights) {
  effects <- lapply(weights$outcome_fits, marginaleffects::avg_predictions, variables = "A")
  theta <- unlist(lapply(effects, \(effect) effect$estimate[2]))

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("direct")

  estimates[, 1] <- theta
  se[, 1] <- unlist(lapply(effects, \(effect) effect$std.error[2]))
  low[, 1] <- unlist(lapply(effects, \(effect) effect$conf.low[2]))
  high[, 1] <- unlist(lapply(effects, \(effect) effect$conf.high[2]))
  p_values[, 1] <- unlist(lapply(effects, \(effect) effect$p.value[2]))

  result <- list(
    estimator = "weightit",
    estimates = estimates,
    se = se,
    alpha = 0.95,
    low = low,
    high = high,
    p_values = p_values
  )
  class(result) <- "direct"
  result
}

#' @importFrom marginaleffects avg_predictions
theta_direct_matchit <- function(task, matches) {
  theta <- numeric(length(task$trt_levels))
  names(theta) <- task$trt_levels

  effects <- lapply(matches$outcome_fits, marginaleffects::avg_predictions, variables = "A")

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("direct")

  estimates[, 1] <- unlist(lapply(effects, \(effect) effect$estimate[2]))
  se[, 1]        <- unlist(lapply(effects, \(effect) effect$std.error[2]))
  p_values[, 1]  <- unlist(lapply(effects, \(effect) effect$p.value[2]))
  low[, 1]       <- unlist(lapply(effects, \(effect) effect$conf.low[2]))
  high[, 1]      <- unlist(lapply(effects, \(effect) effect$conf.high[2]))

  result <- list(
    estimator = "matchit",
    estimates = estimates,
    se = se,
    alpha = 0.95,
    low = low,
    high = high,
    p_values = p_values
  )
  class(result) <- "direct"
  result
}

#' @importFrom stats pnorm qnorm sd
theta_direct_tmle <- function(task, fluctuations, g) {
  theta <- colMeans(fluctuations$Qtilde)

  theta_eif <- eif_direct(task$data[[task$trt]], task$data[[task$outcome]], task$trt_levels, g, fluctuations$Qtilde, theta)

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("direct")
  colnames(p_values) <- "direct"
  estimates[, 1] <- theta

  se[, 1] <- apply(theta_eif, 2, sd) / sqrt(task$n)

  alpha <- 0.95
  low  <- estimates + qnorm((1 - alpha) / 2) * se
  high <- estimates + qnorm(1 - (1 - alpha) / 2) * se

  p_values <- 2 * pnorm(abs((estimates[, 1] - 1) / se[, 1]), lower.tail = FALSE)

  result <- list(
    estimator = "TMLE",
    estimates = estimates,
    se = se,
    alpha = 0.95,
    low = low,
    high = high,
    p_values = p_values
  )
  class(result) <- "direct"
  result
}

theta_direct_sub <- function(task, Qtilde) {
  theta <- colMeans(Qtilde)

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("direct")

  estimates[, 1] <- theta

  result <- list(
    estimator = "sub",
    estimates = estimates,
    se = se,
    alpha = 0.95,
    low = low,
    high = high,
    p_values = p_values
  )
  class(result) <- "direct"
  result
}

theta_direct_pw <- function(task, g) {
  theta <- colMeans(task$data[[task$outcome]] * task$trt_indicator / g$treatment_probs)

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("direct")

  estimates[, 1] <- theta

  result <- list(
    estimator = "pw",
    estimates = estimates,
    se = se,
    alpha = 0.95,
    low = low,
    high = high,
    p_values = p_values
  )
  class(result) <- "direct"
  result
}
