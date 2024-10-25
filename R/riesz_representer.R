riesz_representer <- function(task, learners, full_fits, learner_folds, parameter = "smr", method = "superriesz", torch_params = list(), verbose = FALSE) {
  results <- list()
  for(fold_index in seq_along(task$cv)) {
    if(verbose) cat("Fold", fold_index, "\n")
    results[[fold_index]] <- estimate_riesz_representer(
      task$get_fold(fold_index),
      task$baseline,
      task$trt,
      task$trt_levels,
      method,
      parameter,
      learners,
      full_fits,
      learner_folds,
      torch_params,
      verbose
    )
  }
  combine_riesz_representer(results, task$data, task$trt_levels, task$cv)
}

combine_riesz_representer <- function(results, data, trt_levels, cv) {
  representer <- matrix(0, nrow = nrow(data), ncol = length(trt_levels))
  colnames(representer) <- trt_levels

  for(fold_index in seq_along(cv)) {
    representer[cv[[fold_index]]$validation_set, ] <- results[[fold_index]]$rr
  }
  list(
    rr = representer
  )
}

#' @importFrom SuperRiesz super_riesz
estimate_riesz_representer <- function(data, baseline, trt, trt_levels, method, parameter, learners, full_fits, learner_folds, torch_params, verbose) {
  fits <- list()
  rr <- matrix(0, nrow = nrow(data$validation), ncol = length(trt_levels))
  colnames(rr) <- trt_levels

  if(method == "torch") {
    rr <- estimate_riesz_representer_torch(data, baseline, trt, trt_levels, parameter, torch_params)
  }
  else if(method == "superriesz") {
    for(trt_level in trt_levels) {
      if(verbose) cat("Treatment:", trt_level, "\n")
      trt_indicator <- as.numeric(data$training[[trt]] == trt_level)

      if(parameter == "smr") {
        m   <- \(alpha, data) alpha(data()) * data("weight")[, 1]
        set <- baseline
        data_natural       <- data$training[, set]
        data_shifted       <- data$training[, set]
        data_natural_valid <- data$validation[, set]
        data_shifted_valid <- data$validation[, set]
      }
      else {
        m   <- \(alpha, data) alpha(data("shifted"))
        set <- c(trt, baseline)

        data_natural              <- data$training[, set]
        data_natural[[trt]]       <- as.numeric(data$training[[trt]] == trt_level)

        data_shifted              <- data$training[, set]
        data_shifted[[trt]]       <- 1

        data_natural_valid        <- data$validation[, set]
        data_natural_valid[[trt]] <- as.numeric(data$validation[[trt]] == trt_level)

        data_shifted_valid        <- data$validation[, set]
        data_shifted_valid[[trt]] <- 1
      }

      fit <- super_riesz(
        data = data_natural,
        library = learners,
        alternatives = list(shifted = data_shifted, weight = data.frame(weight = trt_indicator)),
        m = m,
        folds = learner_folds
      )

      if(parameter == "smr") {
        rr[, trt_level] <- predict(fit, data_natural_valid)
      }
      else {
        rr[, trt_level] <- predict(fit, data_shifted_valid)
      }
      #if(any(rr[, trt_level] > 100)) browser()
    }
  }

  list(rr = rr)
}
