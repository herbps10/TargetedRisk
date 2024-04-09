#' @importFrom generics tidy
#' @export
generics::tidy

#' @export
tidy.smr <- function(x, ...) {
  out <- data.frame(estimator = x$estimator,
                    parameter = rep(c("psi1", "psi2", "SMR"), each = nrow(x$estimates)),
                    trt       = rep(rownames(x$estimates), 3),
                    estimate  = c(x$estimates[, 1], x$estimates[, 2], x$estimates[, 3]),
                    std.error = c(x$se[, 1], x$se[, 2], x$se[, 3]),
                    conf.low  = c(x$low[, 1], x$low[, 2], x$low[, 3]),
                    conf.high = c(x$high[, 1], x$high[, 2], x$high[, 3]))
  class(out) <- c("tbl_df", "tbl", "data.frame")
  out
}
