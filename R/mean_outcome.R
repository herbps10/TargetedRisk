mean_outcome <- function(task) {
  result <- matrix(0, nrow = nrow(task$data), ncol = 1)

  for(fold_index in seq_along(task$cv)) {
    data <- task$get_fold(fold_index)
    for(trt_level in task$trt_levels) {
      valid <- task$cv[[fold_index]]$validation_set
      ybar <- mean(data$training[[task$outcome]][data$training[[task$trt]] == trt_level])
      if(ybar == 1) ybar <- 0.9999
      if(ybar == 0) ybar <- 1 - 0.9999
      result[valid, 1][task$data[[task$trt]][valid] == trt_level] <- ybar
    }
  }
  result
}
