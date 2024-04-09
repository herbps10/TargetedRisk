eif_psi1 <- function(w, a, y, trt_prop, g, Qtilde, trt_indicator, theta) {
  eif <- matrix(nrow = nrow(trt_indicator), ncol = ncol(trt_indicator))
  colnames(eif) <- colnames(trt_indicator)

  for(trt_level in colnames(trt_indicator)) {
    eif[, trt_level] <- (a == trt_level) / trt_prop[, trt_level] * (y - Qtilde[, trt_level])
  }
  eif
}

eif_psi2 <- function(w, a, y, trt_prop, g, Qtilde, trt_indicator, theta) {
  N <- length(w)
  eif <- matrix(nrow = nrow(trt_indicator), ncol = ncol(trt_indicator))
  colnames(eif) <- colnames(trt_indicator)

  for(trt_level in colnames(trt_indicator)) {
    eif[, trt_level] <- 1 / trt_prop[, trt_level] * (
      g[, trt_level] * (y - Qtilde[, trt_level]) + (a == trt_level) * (Qtilde[, trt_level] - theta[trt_level])
    )
  }
  eif
}

eif_smr <- function(psi1, psi2, eif1, eif2) {
  matrix(1 / psi2, ncol = ncol(eif1), nrow = nrow(eif1), byrow = TRUE) * eif1 -
    matrix(psi1 / psi2^2, ncol = ncol(eif2), nrow = nrow(eif2), byrow = TRUE) * eif2
}

#' @importFrom stats pnorm qnorm sd
theta_tmle <- function(task, trt_prop, g, fluctuations) {
  psi1 <- colSums(matrix(fluctuations$ybar, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator) / colSums(task$trt_indicator)
  psi2 <- colSums(matrix(fluctuations$Qtilde, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator) / colSums(task$trt_indicator)
  psiSMR <- psi1 / psi2

  psi1_eif <- eif_psi1(task$data[, task$baseline], task$data[[task$trt]], task$data[[task$outcome]], trt_prop, g, fluctuations$ybar, task$trt_indicator, psi1)
  psi2_eif <- eif_psi2(task$data[, task$baseline], task$data[[task$trt]], task$data[[task$outcome]], trt_prop, g, fluctuations$Qtilde, task$trt_indicator, psi2)
  psiSMR_eif <- eif_smr(psi1, psi2, psi1_eif, psi2_eif)

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 3)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("psi1", "psi2", "SMR")
  colnames(p_values) <- "SMR"
  estimates[, 1] <- psi1
  estimates[, 2] <- psi2
  estimates[, 3] <- psiSMR

  se[, 1] <- apply(psi1_eif, 2, sd)   / sqrt(task$n)
  se[, 2] <- apply(psi2_eif, 2, sd)   / sqrt(task$n)
  se[, 3] <- apply(psiSMR_eif, 2, sd) / sqrt(task$n)

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

theta_sub <- function(task, ybar, trt_prop, Qtilde) {
  psi1 <- colSums(matrix(ybar, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator) / colSums(task$trt_indicator)
  psi2 <- colSums(matrix(Qtilde, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator) / colSums(task$trt_indicator)
  psiSMR <- psi1 / psi2

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 3)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("psi1", "psi2", "SMR")

  estimates[, 1] <- psi1
  estimates[, 2] <- psi2
  estimates[, 3] <- psiSMR

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

theta_pw <- function(task, ybar, trt_prop, g) {
  psi1 <- colSums(matrix(ybar, nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * task$trt_indicator) / colSums(task$trt_indicator)
  psi2 <- colSums(matrix(task$data[[task$outcome]], nrow = nrow(task$trt_indicator), ncol = ncol(task$trt_indicator), byrow = FALSE) * g) / colSums(task$trt_indicator)
  psiSMR <- psi1 / psi2

  estimates <- se <- low <- high <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 3)
  p_values <- matrix(NA_real_, nrow = length(task$trt_levels), ncol = 1)
  rownames(estimates) <- rownames(se) <- rownames(low) <- rownames(high) <- rownames(p_values) <- task$trt_levels
  colnames(estimates) <- colnames(se) <- colnames(low) <- colnames(high) <- c("psi1", "psi2", "SMR")

  estimates[, 1] <- psi1
  estimates[, 2] <- psi2
  estimates[, 3] <- psiSMR

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
