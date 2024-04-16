#' @importFrom R6 R6Class
tsmr_Task <- R6::R6Class(
  "tsmr_Task",
  public = list(
    data = NULL,
    n = NULL,
    trt = NULL,
    trt_levels = NULL,
    trt_prop = NULL,
    outcome = NULL,
    baseline = NULL,
    outcome_type = NULL,
    folds = NULL,
    cv = NULL,
    trt_indicator = NULL,
    initialize = function(data, trt, outcome, baseline, outcome_type, folds, use_strata = TRUE) {
      data$tsmr_id <- 1:nrow(data)

      self$data <- data
      self$n <- nrow(data)
      self$trt <- trt
      self$outcome <- outcome
      self$baseline <- baseline
      self$outcome_type <- outcome_type
      self$folds <- folds
      self$trt_levels <- sort(unique(data[[trt]]))
      self$trt_prop <- table(data[[trt]])[self$trt_levels] / self$n

      self$trt_indicator <- matrix(nrow = nrow(data), ncol = length(self$trt_levels))
      colnames(self$trt_indicator) <- self$trt_levels
      for(trt_level in self$trt_levels) {
        self$trt_indicator[, trt_level] <- self$data[[trt]] == trt_level
      }

      self$cv <- cv_setup(self$data, data$tsmr_id, data[[self$trt]], folds, use_strata)
    },
    get_fold = function(fold_index) {
      list(
        training   = self$data[self$cv[[fold_index]]$training_set,   , drop = FALSE],
        validation = self$data[self$cv[[fold_index]]$validation_set, , drop = FALSE]
      )
    }
  )
)
