tmle_direct <- function(task, Qtilde, g) {
  results <- list()
  for(fold_index in seq_along(task$cv)) {
    results[[fold_index]] <- estimate_tmle_direct(
      task$get_fold(fold_index),
      task$outcome,
      task$trt,
      task$trt_levels,
      task$outcome_type,
      extract_fold(Qtilde$predicted_outcomes, task$cv, fold_index),
      extract_fold(g$treatment_probs, task$cv, fold_index)
    )
  }
  combine_tmle_direct(results, task$data, task$trt_levels, task$cv)
}

combine_tmle_direct <- function(results, data, trt_levels, cv) {
  Qtilde_fluctuation <- matrix(0, nrow = nrow(data), ncol = length(trt_levels))

  colnames(Qtilde_fluctuation) <- trt_levels

  for(fold_index in seq_along(cv)) {
    Qtilde_fluctuation[cv[[fold_index]]$validation_set, ] <- results[[fold_index]]$Qtilde_fluctuation
  }
  list(
    Qtilde = Qtilde_fluctuation
  )
}

#' @importFrom stats coef glm plogis qlogis
estimate_tmle_direct <- function(data, outcome, trt, trt_levels, outcome_type, Qtilde, g = NULL) {
  Qtilde_fluctuation  <- matrix(nrow = nrow(data$validation), ncol = length(trt_levels))

  colnames(Qtilde_fluctuation) <- trt_levels

  for(trt_level in trt_levels) {
    clever_covariate_train <- (data$training[[trt]] == trt_level) / g$training[, trt_level]
    clever_covariate_valid <- 1 / g$validation[, trt_level]

    if(outcome_type == "binomial") {
      fit <- glm(data$training[[outcome]] ~ -1 + clever_covariate_train + offset(qlogis(Qtilde$training[, trt_level])), family = "binomial")
      Qtilde_fluctuation[, trt_level] <- plogis(qlogis(Qtilde$validation[, trt_level]) + coef(fit) * clever_covariate_valid)
    }
    else {
      fit <- glm(data$training[[outcome]] ~ -1 + clever_covariate_train + offset(Qtilde$training[, trt_level]))
      eps <- coef(fit)
      Qtilde_fluctuation[, trt_level] <- Qtilde$validation[, trt_level] + eps * clever_covariate_valid
    }
  }

  list(Qtilde_fluctuation = Qtilde_fluctuation)
}
