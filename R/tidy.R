#' @importFrom generics tidy
#' @export
generics::tidy

#' @export
tidy.smr <- function(x, ...) {
  out <- data.frame(estimator = x$estimator,
                    parameter = rep(c("psi1", "psi2", "ER", "SMR"), each = nrow(x$estimates)),
                    trt       = rep(rownames(x$estimates), 4),
                    estimate  = c(x$estimates[, 1], x$estimates[, 2], x$estimates[, 3], x$estimates[, 4]),
                    std.error = c(x$se[, 1], x$se[, 2], x$se[, 3], x$se[, 4]),
                    conf.low  = c(x$low[, 1], x$low[, 2], x$low[, 3], x$low[, 4]),
                    conf.high = c(x$high[, 1], x$high[, 2], x$high[, 3], x$high[, 4]))
  class(out) <- c("tbl_df", "tbl", "data.frame")
  out
}

#' @export
tidy.direct <- function(x, ...) {
  out <- data.frame(estimator = x$estimator,
                    parameter = "direct",
                    trt       = rownames(x$estimates),
                    estimate  = x$estimates[, 1],
                    std.error = x$se[, 1],
                    conf.low  = x$low[, 1],
                    conf.high = x$high[, 1])
  class(out) <- c("tbl_df", "tbl", "data.frame")
  out
}
