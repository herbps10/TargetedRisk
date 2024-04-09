#' @importFrom generics tidy
#' @export
generics::tidy

#' @export
tidy.smr <- function(x, ...) {
  out <- data.frame(estimator = x$estimator,
                    trt = rownames(x$estimates),
                    estimate = x$estimates[, 3],
                    std.error = x$se[, 3],
                    conf.low = x$low[, 3],
                    conf.high = x$high[, 3])
  class(out) <- c("tbl_df", "tbl", "data.frame")
  out
}
