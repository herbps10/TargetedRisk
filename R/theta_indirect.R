eif_theta1 <- function(a, y, trt_prop, Qtilde, trt_indicator, theta) {
  eif <- matrix(nrow = nrow(trt_indicator), ncol = ncol(trt_indicator))
  colnames(eif) <- colnames(trt_indicator)

  for(trt_level in colnames(trt_indicator)) {
    eif[, trt_level] <- (a == trt_level) / trt_prop[, trt_level] * (y - Qtilde[, trt_level])
  }
  eif
}

eif_theta2 <- function(a, y, trt_prop, g, Qtilde, trt_indicator, theta) {
  eif <- matrix(nrow = nrow(trt_indicator), ncol = ncol(trt_indicator))
  colnames(eif) <- colnames(trt_indicator)

  for(trt_level in colnames(trt_indicator)) {
    eif[, trt_level] <- 1 / trt_prop[, trt_level] * (
      g$treatment_probs[, trt_level] * (y - Qtilde[, trt_level]) + (a == trt_level) * (Qtilde[, trt_level] - theta[trt_level])
    )
  }
  eif
}

eif_er <- function(theta1, theta2, eif1, eif2) {
  eif1 - eif2
}

eif_smr <- function(theta1, theta2, eif1, eif2) {
  matrix(1 / theta2, ncol = ncol(eif1), nrow = nrow(eif1), byrow = TRUE) * eif1 -
    matrix(theta1 / theta2^2, ncol = ncol(eif2), nrow = nrow(eif2), byrow = TRUE) * eif2
}


#' @importFrom stats pnorm qnorm sd
theta_indirect_tmle <- function(task, trt_prop, fluctuations, g) {
  theta1 <- colMeans(matrix(fluctuations$ybar, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator / trt_prop)
  theta2 <- colMeans(matrix(fluctuations$Qtilde, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator / trt_prop)
  thetaER  <- theta1 - theta2
  thetaSMR <- theta1 / theta2

  theta1_eif <- eif_theta1(task$data[[task$trt]], task$data[[task$outcome]], trt_prop, fluctuations$ybar, task$trt_indicator, theta1)
  theta2_eif <- eif_theta2(task$data[[task$trt]], task$data[[task$outcome]], trt_prop, g, fluctuations$Qtilde, task$trt_indicator, theta2)
  thetaER_eif <- eif_er(theta1, theta2, theta1_eif, theta2_eif)
  thetaSMR_eif <- eif_smr(theta1, theta2, theta1_eif, theta2_eif)

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 4)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 2)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("theta1", "theta2", "ER", "SMR")
  colnames(p_values) <- c("ER", "SMR")
  estimates[, 1] <- theta1
  estimates[, 2] <- theta2
  estimates[, 3] <- thetaER
  estimates[, 4] <- thetaSMR

  se[, 1] <- apply(theta1_eif, 2, sd) / sqrt(task$n)
  se[, 2] <- apply(theta2_eif, 2, sd) / sqrt(task$n)
  se[, 3] <- apply(thetaER_eif, 2, sd) / sqrt(task$n)
  se[, 4] <- apply(thetaSMR_eif, 2, sd) / sqrt(task$n)

  alpha <- 0.95
  low  <- estimates + qnorm((1 - alpha) / 2) * se
  high <- estimates + qnorm(1 - (1 - alpha) / 2) * se

  p_values[, 1] <- 2 * pnorm(abs((estimates[, 3]) / se[, 3]), lower.tail = FALSE)
  p_values[, 2] <- 2 * pnorm(abs((estimates[, 4] - 1) / se[, 4]), lower.tail = FALSE)

  result <- list(
    estimator = "TMLE",
    estimates = estimates,
    se = se,
    alpha = 0.95,
    low = low,
    high = high,
    p_values = p_values
  )
  class(result) <- "smr"
  result
}


#' @importFrom stats pnorm qnorm sd
theta_indirect_onestep <- function(task, trt_prop, ybar, Qtilde, g) {
  theta1 <- colMeans(matrix(ybar, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator / trt_prop)
  theta2 <- colMeans(matrix(Qtilde, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator / trt_prop)
  thetaER  <- theta1 - theta2
  thetaSMR <- theta1 / theta2

  Qtilde_mat <- matrix(Qtilde, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE)
  colnames(Qtilde_mat) <- task$trt_levels

  ybar_mat <- matrix(ybar, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE)
  colnames(ybar_mat) <- task$trt_levels

  theta1_eif <- eif_theta1(task$data[[task$trt]], task$data[[task$outcome]], trt_prop, ybar_mat, task$trt_indicator, theta1)
  theta2_eif <- eif_theta2(task$data[[task$trt]], task$data[[task$outcome]], trt_prop, g, Qtilde_mat, task$trt_indicator, theta2)
  thetaER_eif <- eif_er(theta1, theta2, theta1_eif, theta2_eif)
  thetaSMR_eif <- eif_smr(theta1, theta2, theta1_eif, theta2_eif)

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 4)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 2)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("theta1", "theta2", "ER", "SMR")
  colnames(p_values) <- c("ER", "SMR")
  estimates[, 1] <- theta1
  estimates[, 2] <- theta2
  estimates[, 3] <- thetaER
  estimates[, 4] <- thetaSMR

  se[, 1] <- apply(theta1_eif, 2, sd) / sqrt(task$n)
  se[, 2] <- apply(theta2_eif, 2, sd) / sqrt(task$n)
  se[, 3] <- apply(thetaER_eif, 2, sd) / sqrt(task$n)
  se[, 4] <- apply(thetaSMR_eif, 2, sd) / sqrt(task$n)

  alpha <- 0.95
  low  <- estimates + qnorm((1 - alpha) / 2) * se
  high <- estimates + qnorm(1 - (1 - alpha) / 2) * se

  p_values[, 1] <- 2 * pnorm(abs((estimates[, 3]) / se[, 3]), lower.tail = FALSE)
  p_values[, 2] <- 2 * pnorm(abs((estimates[, 4] - 1) / se[, 4]), lower.tail = FALSE)

  result <- list(
    estimator = "Onestep",
    estimates = estimates,
    se = se,
    alpha = 0.95,
    low = low,
    high = high,
    p_values = p_values
  )
  class(result) <- "smr"
  result
}

theta_indirect_sub <- function(task, ybar, trt_prop, Qtilde) {
  theta1 <- colSums(matrix(ybar, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator) / colSums(task$trt_indicator)
  theta2 <- colSums(matrix(Qtilde, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator) / colSums(task$trt_indicator)
  thetaER <- theta1 - theta2
  thetaSMR <- theta1 / theta2

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 4)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 2)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("theta1", "theta2", "ER", "SMR")

  estimates[, 1] <- theta1
  estimates[, 2] <- theta2
  estimates[, 3] <- thetaER
  estimates[, 4] <- thetaSMR

  result <- list(
    estimator = "sub",
    estimates = estimates,
    se = se,
    alpha = 0.95,
    low = low,
    high = high,
    p_values = p_values
  )
  class(result) <- "smr"
  result
}

theta_indirect_pw <- function(task, ybar, trt_prop, g) {
  theta1 <- colSums(matrix(ybar, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator) / colSums(task$trt_indicator)
  theta2 <- colSums(matrix(task$data[[task$outcome]], nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * g$treatment_probs) / colSums(task$trt_indicator)
  thetaER <- theta1 - theta2
  thetaSMR <- theta1 / theta2

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 4)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 2)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("theta1", "theta2", "ER", "SMR")

  estimates[, 1] <- theta1
  estimates[, 2] <- theta2
  estimates[, 3] <- thetaER
  estimates[, 4] <- thetaSMR

  p_values[, 1] <- 2 * pnorm(abs((estimates[, 3]) / se[, 3]), lower.tail = FALSE)
  p_values[, 2] <- 2 * pnorm(abs((estimates[, 4] - 1) / se[, 4]), lower.tail = FALSE)

  result <- list(
    estimator = "pw",
    estimates = estimates,
    se = se,
    alpha = 0.95,
    low = low,
    high = high,
    p_values = p_values
  )
  class(result) <- "smr"
  result
}

#' @importFrom marginaleffects avg_predictions
theta_indirect_matchit <- function(task, ybar, trt_prop, matches) {
  ybar_mat <- matrix(ybar, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE)
  colnames(ybar_mat) <- task$trt_levels
  theta1 <- colSums(matrix(ybar, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator) / colSums(task$trt_indicator)
  theta1_eif <- eif_theta1(task$data[[task$trt]], task$data[[task$outcome]], trt_prop, ybar_mat, task$trt_indicator, theta1)

  effects <- lapply(matches$outcome_fits, marginaleffects::avg_predictions, variables = "A", newdata = subset(A == 1))

  theta2 <- unlist(lapply(effects, \(effect) effect$estimate[1]))

  thetaER <- theta1 - theta2
  thetaSMR <- theta1 / theta2

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 4)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 4)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("theta1", "theta2", "ER", "SMR")

  estimates[, 1] <- theta1
  estimates[, 2] <- theta2
  estimates[, 3] <- thetaER
  estimates[, 4] <- thetaSMR


  se[, 1] <- apply(theta1_eif, 2, sd) / sqrt(task$n)
  se[, 2] <- unlist(lapply(effects, \(effect) effect$std.error[1]))

  se[, 3] <- sqrt(se[, 1]^2 + se[, 2]^2)

  alpha <- 0.95
  low  <- estimates + qnorm((1 - alpha) / 2) * se
  high <- estimates + qnorm(1 - (1 - alpha) / 2) * se

  result <- list(
    estimator = "matchit",
    estimates = estimates,
    se = se,
    alpha = 0.95,
    low = low,
    high = high,
    p_values = p_values
  )
  class(result) <- "smr"
  result
}


#' @importFrom marginaleffects avg_predictions
theta_indirect_weightit <- function(task, ybar, trt_prop, weights) {
  ybar_mat <- matrix(ybar, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE)
  colnames(ybar_mat) <- task$trt_levels
  theta1 <- colSums(matrix(ybar, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator) / colSums(task$trt_indicator)
  theta1_eif <- eif_theta1(task$data[[task$trt]], task$data[[task$outcome]], trt_prop, ybar_mat, task$trt_indicator, theta1)

  effects <- lapply(weights$outcome_fits, marginaleffects::avg_predictions, variables = "A", newdata = subset(A == 1))

  theta2 <- unlist(lapply(effects, \(effect) effect$estimate[1]))

  thetaER <- theta1 - theta2
  thetaSMR <- theta1 / theta2

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 4)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 4)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("theta1", "theta2", "ER", "SMR")

  estimates[, 1] <- theta1
  estimates[, 2] <- theta2
  estimates[, 3] <- thetaER
  estimates[, 4] <- thetaSMR


  se[, 1] <- apply(theta1_eif, 2, sd) / sqrt(task$n)
  se[, 2] <- unlist(lapply(effects, \(effect) effect$std.error[1]))

  se[, 3] <- sqrt(se[, 1]^2 + se[, 2]^2)

  alpha <- 0.95
  low  <- estimates + qnorm((1 - alpha) / 2) * se
  high <- estimates + qnorm(1 - (1 - alpha) / 2) * se

  result <- list(
    estimator = "matchit",
    estimates = estimates,
    se = se,
    alpha = 0.95,
    low = low,
    high = high,
    p_values = p_values
  )
  class(result) <- "smr"
  result
}

