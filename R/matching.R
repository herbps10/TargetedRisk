#' @importFrom MatchIt matchit
#' @importFrom stats as.formula
matching <- function(task, method = "nearerst", distance = "glm") {
  fits <- list()

  for(trt_level in task$trt_levels) {
    f <- as.formula(paste0("task$trt_indicator[, trt_level] ~", paste0(task$baseline, collapse = " + ")))
    fits[[trt_level]] <- matchit(formula = f, data = task$data, method = method, distance = distance)
  }
  fits
}
