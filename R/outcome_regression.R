outcome_regression <- function(task, learners, full_fits, learner_folds) {
  results <- list()
  for(fold_index in seq_along(task$cv)) {
    results[[fold_index]] <- estimate_outcome_regression(
      task$get_fold(fold_index),
      task$baseline,
      task$outcome,
      task$outcome_type,
      learners,
      full_fits,
      learner_folds
    )
  }
  combine_outcome_regressions(results, task$data, task$trt_levels, task$cv)
}

combine_outcome_regressions <- function(results, data, trt_levels, cv) {
  predicted_outcomes <- matrix(0, nrow = nrow(data), ncol = 1)
  for(fold_index in seq_along(cv)) {
    predicted_outcomes[cv[[fold_index]]$validation_set, ] <- results[[fold_index]]$predicted_outcomes
  }
  list(
    predicted_outcomes = predicted_outcomes
  )
}

estimate_outcome_regression <- function(data, baseline, outcome, outcome_type, learners, full_fits, learner_folds) {
  fits <- list()
  predicted_outcomes <- matrix(0, nrow = nrow(data$validation), ncol = 1)

  fit <- superlearner(
    data = data$training[, c(baseline, outcome)],
    outcome = outcome,
    outcome_type = outcome_type,
    learners = learners,
    learner_folds = learner_folds
  )

  if(full_fits) {
    fits[[trt_level]] <- fit
  }

  predicted_outcomes <- predict(fit, data$validation)
  predicted_outcomes <- ifelse(predicted_outcomes == 0, 0.0001, ifelse(predicted_outcomes == 1, 0.9999, predicted_outcomes))

  list(predicted_outcomes = predicted_outcomes, fits = fits)
}
