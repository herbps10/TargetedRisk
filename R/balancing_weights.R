#' @importFrom WeightIt weightit
#' @importFrom stats as.formula
balancing_weights <- function(task, method = "glm") {
  fits <- list()

  for(trt_level in task$trt_levels) {
    f <- as.formula(paste0("task$trt_indicator[, trt_level] ~", paste0(task$baseline, collapse = " + ")))
    fits[[trt_level]] <- weightit(formula = f, data = task$data, estimand = "ATE", method = method)
  }
  fits
}
