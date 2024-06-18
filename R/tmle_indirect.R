tmle_indirect <- function(task, ybar, trt_prop, Qtilde, g = NULL, riesz = NULL) {
  results <- list()
  for(fold_index in seq_along(task$cv)) {
    results[[fold_index]] <- estimate_tmle(
      task$get_fold(fold_index),
      task$outcome,
      task$trt,
      task$trt_levels,
      extract_fold(ybar, task$cv, fold_index),
      extract_fold(trt_prop, task$cv, fold_index),
      extract_fold(Qtilde$predicted_outcomes, task$cv, fold_index),
      if(is.null(g)) { NULL } else { extract_fold(g$treatment_probs, task$cv, fold_index) },
      if(is.null(riesz)) { NULL } else { extract_fold(riesz$rr, task$cv, fold_index) }
    )
  }
  combine_tmle_indirect(results, task$data, task$trt_levels, task$cv)
}

combine_tmle_indirect <- function(results, data, trt_levels, cv) {
  ybar_fluctuation   <- matrix(0, nrow = nrow(data), ncol = length(trt_levels))
  Qtilde_fluctuation <- matrix(0, nrow = nrow(data), ncol = length(trt_levels))

  colnames(ybar_fluctuation) <- trt_levels
  colnames(Qtilde_fluctuation) <- trt_levels

  for(fold_index in seq_along(cv)) {
    ybar_fluctuation[cv[[fold_index]]$validation_set, ] <- results[[fold_index]]$ybar_fluctuation
    Qtilde_fluctuation[cv[[fold_index]]$validation_set, ] <- results[[fold_index]]$Qtilde_fluctuation
  }
  list(
    ybar = ybar_fluctuation,
    Qtilde = Qtilde_fluctuation
  )
}

#' @importFrom stats coef glm plogis qlogis
estimate_tmle <- function(data, outcome, trt, trt_levels, ybar, trt_prop, Qtilde, g = NULL, riesz = NULL) {
  ybar_fluctuation <- matrix(nrow = nrow(data$validation), ncol = length(trt_levels))
  Qtilde_fluctuation  <- matrix(nrow = nrow(data$validation), ncol = length(trt_levels))

  colnames(ybar_fluctuation) <- trt_levels
  colnames(Qtilde_fluctuation) <- trt_levels

  for(trt_level in trt_levels) {
    # Psi 1
    clever_covariate_train1 <- as.numeric(data$training[[trt]] == trt_level) / trt_prop$training[, trt_level]
    clever_covariate_valid1 <- as.numeric(data$validation[[trt]] == trt_level) / trt_prop$validation[, trt_level]
    fit1 <- glm(data$training[[outcome]] ~ -1 + clever_covariate_train1 + offset(qlogis(ybar$training[, 1])), family = "binomial")
    ybar_fluctuation[, trt_level] <- plogis(qlogis(ybar$validation) + coef(fit1) * clever_covariate_valid1)

    # Psi 2
    if(is.null(g)) {
      clever_covariate_train2 <- riesz$training[, trt_level]
      clever_covariate_valid2 <- riesz$validation[, trt_level]
    }
    else {
      clever_covariate_train2 <- g$training[, trt_level]
      clever_covariate_valid2 <- g$validation[, trt_level]
    }

    fit2 <- glm(data$training[[outcome]] ~ -1 + clever_covariate_train2 + offset(qlogis(Qtilde$training)), family = "binomial")
    Qtilde_fluctuation[, trt_level] <- plogis(qlogis(Qtilde$validation) + coef(fit2) * clever_covariate_valid2)
  }

  list(ybar_fluctuation = ybar_fluctuation, Qtilde_fluctuation = Qtilde_fluctuation)
}
