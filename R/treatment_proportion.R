treatment_proportion <- function(task) {
  result <- matrix(0, nrow = nrow(task$data), ncol = length(task$trt_levels))
  colnames(result) <- task$trt_levels

  for(fold_index in seq_along(task$cv)) {
    data <- task$get_fold(fold_index)
    for(trt_level in task$trt_levels) {
      result[task$cv[[fold_index]]$validation_set, trt_level] <- mean(data$training[, task$trt] == trt_level)
    }
  }
  result
}
