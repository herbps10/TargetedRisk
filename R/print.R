#' @importFrom glue glue
print_smr_table <- function(x, param) {
  f <- function(x, digits = 3) round(x, digits)
  if(param == "SMR") {
    cat("Trt \t Est \t SE \t 95% CI \t p\n")
  }
  else {
    cat("Trt \t Est \t SE \t 95% CI\n")
  }
  for(trt in seq_along(rownames(x$estimates))) {
    text <- glue::glue("{trt} \t {f(x$estimates[trt, param])} \t {f(x$se[trt, param])} \t ({f(x$low[trt, param])}, {f(x$high[trt, param])})")
    if(param == "SMR") {
      text <- paste0(text, glue::glue("\t {f(x$p_values[trt], 3)}"))
    }
    text <- paste0(text, "\n")
    cat(text)
  }
}

#' @importFrom cli cli_text
#' @export
print.smr <- function(x, ...) {
  estimator <- ifelse(x$estimator == "TMLE", "TMLE", ifelse(x$estimator == "sub", "Substitution", "Probability weighted"))
  cli::cli_text("{.strong Indirect Standardization Estimator}: {estimator}")
  cat("\n")
  cli::cli_text("{.strong Psi1}")
  print_smr_table(x, "psi1")
  cat("\n")
  cli::cli_text("{.strong Psi2}")
  print_smr_table(x, "psi2")
  cat("\n")
  cli::cli_text("{.strong Standardized Mortality Ratio (SMR)}")
  print_smr_table(x, "SMR")
  cat("\n")
}

#' @importFrom glue glue
print_direct_table <- function(x, param) {
  f <- function(x, digits = 3) round(x, digits)
  if(param == "SMR") {
    cat("Trt \t Est \t SE \t 95% CI \t p\n")
  }
  else {
    cat("Trt \t Est \t SE \t 95% CI\n")
  }
  for(trt in seq_along(rownames(x$estimates))) {
    text <- glue::glue("{trt} \t {f(x$estimates[trt, param])} \t {f(x$se[trt, param])} \t ({f(x$low[trt, param])}, {f(x$high[trt, param])})")
    if(param == "SMR") {
      text <- paste0(text, glue::glue("\t {f(x$p_values[trt], 3)}"))
    }
    text <- paste0(text, "\n")
    cat(text)
  }
}

#' @importFrom cli cli_text
#' @export
print.direct <- function(x, ...) {
  estimator <- ifelse(x$estimator == "TMLE", "TMLE", ifelse(x$estimator == "sub", "Substitution", "Probability weighted"))
  cli::cli_text("{.strong Direct Standardization Estimator}: {estimator}")
  cat("\n")
  cli::cli_text("{.strong Psi}")
  print_direct_table(x, "direct")
  cat("\n")
}
