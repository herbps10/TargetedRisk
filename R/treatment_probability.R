#' @importFrom future future
treatment_probability <- function(task, learners, full_fits, learner_folds, verbose = FALSE) {
  results <- list()
  for(fold_index in seq_along(task$cv)) {
    if(verbose == TRUE) cat(paste0("Fold: ", fold_index, "\n"))
    results[[fold_index]] <- future::future({
      estimate_treatment_probability(
        task$get_fold(fold_index),
        task$baseline,
        task$trt,
        task$trt_levels,
        learners,
        full_fits,
        learner_folds,
        verbose
      )
    }, packages = c("TargetedRisk", "mlr3extralearners"))
  }

  results <- future::value(results)

  combine_treatment_probabilities(results, task$data, task$trt_levels, task$cv)
}

combine_treatment_probabilities <- function(results, data, trt_levels, cv) {
  treatment_probs <- matrix(0, nrow = nrow(data), ncol = length(trt_levels))
  colnames(treatment_probs) <- trt_levels

  for(fold_index in seq_along(cv)) {
    treatment_probs[cv[[fold_index]]$validation_set, ] <- results[[fold_index]]$treatment_probs
  }
  list(
    treatment_probs = treatment_probs
  )
}

estimate_treatment_probability <- function(data, baseline, trt, trt_levels, learners, full_fits, learner_folds, verbose = FALSE) {
  fits <- list()
  treatment_probs <- matrix(0, nrow = nrow(data$validation), ncol = length(trt_levels))
  colnames(treatment_probs) <- trt_levels
  for(trt_level in trt_levels) {
    if(verbose == TRUE) cat(paste0("Treatment: ", trt_level, "\n"))
    data$training$trt_indicator <- as.numeric(data$training[[trt]] == trt_level)
    fit <- superlearner(
      data = data$training[, c(baseline, "trt_indicator")],
      outcome = "trt_indicator",
      outcome_type = "binomial",
      learners = learners,
      learner_folds = learner_folds
    )
    print(fit)

    if(full_fits) {
      fits[[trt_level]] <- fit
    }

    treatment_probs[, trt_level] <- predict(fit, data$validation)
    treatment_probs[treatment_probs[, trt_level] == 0, trt_level] <- min(treatment_probs[treatment_probs[, trt_level] > 0, trt_level])
  }

  list(treatment_probs = treatment_probs, fits = fits)
}
