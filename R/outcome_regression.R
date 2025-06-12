outcome_regression <- function(task, learners, include_treatment, full_fits, learner_folds) {
  results <- list()
  for(fold_index in seq_along(task$cv)) {
    results[[fold_index]] <- estimate_outcome_regression(
      task$get_fold(fold_index),
      task$trt,
      task$baseline,
      task$outcome,
      task$outcome_type,
      task$trt_levels,
      learners,
      include_treatment,
      full_fits,
      learner_folds
    )
  }
  combine_outcome_regressions(results, task$data, task$trt_levels, include_treatment, task$cv)
}

combine_outcome_regressions <- function(results, data, trt_levels, include_treatment, cv) {
  predicted_outcomes <- matrix(0, nrow = nrow(data), ncol = ifelse(include_treatment, length(trt_levels), 1))
  if(include_treatment == TRUE) {
    colnames(predicted_outcomes) <- trt_levels
  }
  for(fold_index in seq_along(cv)) {
    predicted_outcomes[cv[[fold_index]]$validation_set, ] <- results[[fold_index]]$predicted_outcomes
  }
  list(
    predicted_outcomes = predicted_outcomes
  )
}

#' @importFrom stats predict
estimate_outcome_regression <- function(data, trt, baseline, outcome, outcome_type, trt_levels, learners, include_treatment, full_fits, learner_folds) {
  fits <- list()
  predicted_outcomes <- matrix(0, nrow = nrow(data$validation), ncol = 1)

  if(include_treatment == TRUE) {
    onehot <- model.matrix(~-1 + trt, data = data.frame(trt = factor(data$training[[trt]], levels = trt_levels)))
    train <- cbind(data$training, onehot)
    set <- c(c(baseline, outcome), colnames(onehot))
  }
  else {
    train <- data$training
    set <- c(baseline, outcome)
  }

  fit <- superlearner(
    data = train[, set, drop = FALSE],
    outcome = outcome,
    outcome_type = outcome_type,
    learners = learners,
    learner_folds = learner_folds
  )

  if(full_fits) {
    fits[[1]] <- fit
  }

  if(include_treatment == TRUE) {
    predicted_outcomes <- matrix(0, ncol = length(trt_levels), nrow = nrow(data$validation))
    colnames(predicted_outcomes) <- trt_levels
    for(trt_level in trt_levels) {
      data$validation[[trt]] <- trt_level
      onehot <- model.matrix(~-1 + trt, data = data.frame(trt = factor(data$validation[[trt]], levels = trt_levels)))
      valid <- cbind(data$validation, onehot)
      predicted_outcomes[, trt_level] <- predict(fit, valid)
      predicted_outcomes[, trt_level] <- ifelse(predicted_outcomes[, trt_level] == 0, 0.0001, ifelse(predicted_outcomes[, trt_level] == 0, 0.9999, predicted_outcomes[, trt_level]))
    }
  }
  else {
    predicted_outcomes <- predict(fit, data$validation)
    predicted_outcomes <- ifelse(predicted_outcomes == 0, 0.0001, ifelse(predicted_outcomes == 1, 0.9999, predicted_outcomes))
  }

  list(predicted_outcomes = predicted_outcomes, fits = fits)
}
