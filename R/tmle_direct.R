tmle_direct <- function(task, Qtilde, g = NULL, riesz = NULL) {
  results <- list()
  for(fold_index in seq_along(task$cv)) {
    results[[fold_index]] <- estimate_tmle_direct(
      task$get_fold(fold_index),
      task$outcome,
      task$trt,
      task$trt_levels,
      extract_fold(Qtilde$predicted_outcomes, task$cv, fold_index),
      if(is.null(g)) { NULL } else { extract_fold(g$treatment_probs, task$cv, fold_index) },
      if(is.null(riesz)) { NULL } else { extract_fold(riesz$rr, task$cv, fold_index) }
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
estimate_tmle_direct <- function(data, outcome, trt, trt_levels, Qtilde, g = NULL, riesz = NULL) {
  Qtilde_fluctuation  <- matrix(nrow = nrow(data$validation), ncol = length(trt_levels))

  colnames(Qtilde_fluctuation) <- trt_levels

  for(trt_level in trt_levels) {
    if(is.null(g)) {
      clever_covariate_train <- (data$training[[trt]] == trt_level) * riesz$training[, trt_level]
      clever_covariate_valid <- (data$validation[[trt]] == trt_level) * riesz$validation[, trt_level]
    }
    else {
      clever_covariate_train <- (data$training[[trt]] == trt_level) / g$training[, trt_level]
      clever_covariate_valid <- (data$validation[[trt]] == trt_level) / g$validation[, trt_level]
    }

    fit <- glm(data$training[[outcome]] ~ -1 + clever_covariate_train + offset(qlogis(Qtilde$training[, trt_level])), family = "binomial")
    Qtilde_fluctuation[, trt_level] <- plogis(qlogis(Qtilde$validation[, trt_level]) + coef(fit) * clever_covariate_valid)
  }

  list(Qtilde_fluctuation = Qtilde_fluctuation)
}
