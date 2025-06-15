#' @importFrom MatchIt matchit
#' @importFrom stats as.formula lm
matching <- function(task, method = "nearerst", distance = "glm", estimand = "ATT") {
  fits <- list()
  outcome_fits <- list()

  for(trt_level in task$trt_levels) {
    f <- as.formula(paste0("task$trt_indicator[, trt_level] ~", paste0(task$baseline, collapse = " + ")))
    fits[[trt_level]] <- matchit(formula = f, data = task$data, method = method, distance = distance, estimand = estimand)

    df <- data.frame(Y = task$data[[task$outcome]], A = task$trt_indicator[, trt_level])
    outcome_fits[[trt_level]] <- lm(Y ~ A, data = df, weights = fits[[trt_level]]$weights)
  }
  list(fits = fits, outcome_fits = outcome_fits)
}
