eif_theta1 <- function(a, y, trt_prop, Qtilde, trt_indicator, theta) {
  eif <- matrix(nrow = nrow(trt_indicator), ncol = ncol(trt_indicator))
  colnames(eif) <- colnames(trt_indicator)

  for(trt_level in colnames(trt_indicator)) {
    eif[, trt_level] <- (a == trt_level) / trt_prop[, trt_level] * (y - Qtilde[, trt_level])
  }
  eif
}

eif_theta2 <- function(a, y, trt_prop, g, riesz, Qtilde, trt_indicator, theta) {
  eif <- matrix(nrow = nrow(trt_indicator), ncol = ncol(trt_indicator))
  colnames(eif) <- colnames(trt_indicator)

  for(trt_level in colnames(trt_indicator)) {
    if(is.null(g)) {
      eif[, trt_level] <- 1 / trt_prop[, trt_level] * (
        riesz$rr[, trt_level] * (y - Qtilde[, trt_level]) + (a == trt_level) * (Qtilde[, trt_level] - theta[trt_level])
      )
    }
    else {
      eif[, trt_level] <- 1 / trt_prop[, trt_level] * (
        g$treatment_probs[, trt_level] * (y - Qtilde[, trt_level]) + (a == trt_level) * (Qtilde[, trt_level] - theta[trt_level])
      )
    }
  }
  eif
}

eif_smr <- function(theta1, theta2, eif1, eif2) {
  matrix(1 / theta2, ncol = ncol(eif1), nrow = nrow(eif1), byrow = TRUE) * eif1 -
    matrix(theta1 / theta2^2, ncol = ncol(eif2), nrow = nrow(eif2), byrow = TRUE) * eif2
}

#' @importFrom stats pnorm qnorm sd
theta_indirect_tmle <- function(task, trt_prop, fluctuations, g, riesz) {
  theta1 <- colSums(matrix(fluctuations$ybar, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator) / colSums(task$trt_indicator)
  theta2 <- colSums(matrix(fluctuations$Qtilde, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator) / colSums(task$trt_indicator)
  thetaSMR <- theta1 / theta2

  theta1_eif <- eif_theta1(task$data[[task$trt]], task$data[[task$outcome]], trt_prop, fluctuations$ybar, task$trt_indicator, theta1)
  theta2_eif <- eif_theta2(task$data[[task$trt]], task$data[[task$outcome]], trt_prop, g, riesz, fluctuations$Qtilde, task$trt_indicator, theta2)
  thetaSMR_eif <- eif_smr(theta1, theta2, theta1_eif, theta2_eif)

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 3)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("theta1", "theta2", "SMR")
  colnames(p_values) <- "SMR"
  estimates[, 1] <- theta1
  estimates[, 2] <- theta2
  estimates[, 3] <- thetaSMR

  se[, 1] <- apply(theta1_eif, 2, sd)   / sqrt(task$n)
  se[, 2] <- apply(theta2_eif, 2, sd)   / sqrt(task$n)
  se[, 3] <- apply(thetaSMR_eif, 2, sd) / sqrt(task$n)

  alpha <- 0.95
  low  <- estimates + qnorm((1 - alpha) / 2) * se
  high <- estimates + qnorm(1 - (1 - alpha) / 2) * se

  p_values <- 2 * pnorm(abs((estimates[, 3] - 1) / se[, 3]), lower.tail = FALSE)

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

theta_indirect_sub <- function(task, ybar, trt_prop, Qtilde) {
  theta1 <- colSums(matrix(ybar, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator) / colSums(task$trt_indicator)
  theta2 <- colSums(matrix(Qtilde, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator) / colSums(task$trt_indicator)
  thetaSMR <- theta1 / theta2

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 3)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("theta1", "theta2", "SMR")

  estimates[, 1] <- theta1
  estimates[, 2] <- theta2
  estimates[, 3] <- thetaSMR

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

theta_indirect_pw <- function(task, ybar, trt_prop, g, riesz) {
  theta1 <- colSums(matrix(ybar, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator) / colSums(task$trt_indicator)
  if(is.null(g)) {
    theta2 <- colSums(matrix(task$data[[task$outcome]], nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * riesz$rr) / colSums(task$trt_indicator)
  }
  else {
    theta2 <- colSums(matrix(task$data[[task$outcome]], nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * g$treatment_probs) / colSums(task$trt_indicator)
  }
  thetaSMR <- theta1 / theta2

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 3)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("theta1", "theta2", "SMR")

  estimates[, 1] <- theta1
  estimates[, 2] <- theta2
  estimates[, 3] <- thetaSMR

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
