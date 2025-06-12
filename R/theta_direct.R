eif_direct <- function(a, y, trt_levels, g, riesz, Qtilde, theta) {
  eif <- matrix(nrow = length(a), ncol = length(trt_levels))
  colnames(eif) <- trt_levels

  for(trt_level in trt_levels) {
    if(is.null(g)) {
      eif[, trt_level] <- (a == trt_level) * riesz$rr[, trt_level] * (y - Qtilde[, trt_level]) + Qtilde[, trt_level] - theta[trt_level]
    }
    else {
      eif[, trt_level] <- (a == trt_level) / g$treatment_probs[, trt_level] * (y - Qtilde[, trt_level]) + Qtilde[, trt_level] - theta[trt_level]
    }
  }
  eif
}

#' @importFrom stats pnorm qnorm sd
theta_direct_onestep <- function(task, Qtilde, g, riesz) {
  theta <- colMeans(Qtilde)
  theta_eif <- eif_direct(task$data[[task$trt]], task$data[[task$outcome]], task$trt_levels, g, riesz, Qtilde, theta)

  # Apply one-step correction
  theta <- theta + colMeans(theta_eif)
  theta_eif <- eif_direct(task$data[[task$trt]], task$data[[task$outcome]], task$trt_levels, g, riesz, Qtilde, theta)

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

theta_direct_weightit <- function(task, weights) {
  theta <- numeric(length(task$trt_levels))
  names(theta) <- task$trt_levels
  for(trt_level in task$trt_levels) {
    trt_indicator <- task$trt_indicator[, trt_level]
    theta[trt_level] <- mean((weights[[trt_level]]$weights * task$data[[task$outcome]])[trt_indicator])
  }

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("direct")

  estimates[, 1] <- theta

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

#' @importFrom stats pnorm qnorm sd
theta_direct_tmle <- function(task, fluctuations, g, riesz) {
  theta <- colMeans(fluctuations$Qtilde)

  theta_eif <- eif_direct(task$data[[task$trt]], task$data[[task$outcome]], task$trt_levels, g, riesz, fluctuations$Qtilde, theta)

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

theta_direct_pw <- function(task, g, riesz) {
  if(is.null(g)) {
    theta <- colMeans(task$data[[task$outcome]] * task$trt_indicator * riesz$rr)
  }
  else {
    theta <- colMeans(task$data[[task$outcome]] * task$trt_indicator / g$treatment_probs)
  }

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
