#' Targeted Minimum Loss-Based Estimator for Direct Standardization
#'
#' @param data \[\code{data.frame}\]\cr
#' A \code{data.frame} containing all baseline, treatment, and outcome variables.
#' @param trt \[\code{character}\]\cr
#' Column name of treatment variable.
#' @param outcome \[\code{character}\]\cr
#' Column name of outcome variable.
#' @param baseline \[\code{character}\]\cr
#' Vector of column names of baseline variables.
#' @param outcome_type \[\code{character}\]\cr
#' Outcome variable type: binomial (binary) or continuous.
#' @param folds \[\code{integer(1)}\]\cr
#' Number of folds for top level cross-fitting.
#' @param learners_trt \[\code{character}\]\cr
#' Vector of learners to include in SuperLearner library for estimating treatment assignment mechanism.
#' @param learners_outcome \[\code{character}\]\cr
#' Vector of learners to include in SuperLearner library for estimating outcome regression.
#' @param verbose \[\code{logical}]\cr
#' Whether to print information messages during fitting
#' @param control \[\code{standardization_control}\]\cr
#' Additional tuning parameters for controlling fitting. Specify using \link{standardization_control}.
#'
#' @return A list of class \code{smr}
#'
#' @export
direct_tmle <- function(data, trt, outcome, baseline, outcome_type = c("binomial", "continuous"), trt_method = "default", folds = 5, learners_trt = c("mean", "glm"), learners_outcome = c("mean", "glm"), verbose = FALSE, control = standardization_control(), torch_params = list()) {
  if(length(outcome_type) > 1) outcome_type <- outcome_type[1]

  task <- tsmr_Task$new(
    data = data,
    trt = trt,
    outcome = outcome,
    baseline = baseline,
    outcome_type = outcome_type,
    folds = folds
  )

  g <- NULL
  riesz <- NULL
  if(verbose == TRUE) cat("Starting treatment fitting\n")
  if(trt_method == "default") {
    g <- treatment_probability(task, learners_trt, control$.return_full_fits, control$.learners_trt_folds, verbose)
  }
  else if(tolower(trt_method) == "superriesz") {
    riesz <- riesz_representer(task, learners_trt, control$.return_full_fits, control$.learners_trt_folds, parameter = "direct", method = "superriesz")
  }
  else if(tolower(trt_method) == "torch") {
    riesz <- riesz_representer(task, learners_trt, control$.return_full_fits, control$.learners_trt_folds, parameter = "direct", method = "torch", torch_params = control$.torch_params)
  }

  if(verbose == TRUE) cat("Starting outcome fitting\n")
  Qtilde <- outcome_regression(task, learners_outcome, include_treatment = TRUE, control$.return_full_fits, control$.learners_outcome_folds)

  fluctuations <- tmle_direct(task, Qtilde, g, riesz)

  theta <- theta_direct_tmle(task, fluctuations, g, riesz)
  theta$g <- g
  theta$riesz <- riesz
  theta$Qtilde <- Qtilde

  theta
}

#' Substitution Estimator for Direct Standardization
#'
#' @param data \[\code{data.frame}\]\cr
#' A \code{data.frame} containing all baseline, treatment, and outcome variables.
#' @param trt \[\code{character}\]\cr
#' Column name of treatment variable.
#' @param outcome \[\code{character}\]\cr
#' Column name of outcome variable.
#' @param baseline \[\code{character}\]\cr
#' Vector of column names of baseline variables.
#' @param outcome_type \[\code{character}\]\cr
#' Outcome variable type: binomial (binary) or continuous.
#' @param folds \[\code{integer(1)}\]\cr
#' Number of folds for top level cross-fitting.
#' @param learners \[\code{character}\]\cr
#' Vector of learners to include in SuperLearner library for estimating outcome regression.
#' @param control \[\code{standardization_control}\]\cr
#' Additional tuning parameters for controlling fitting. Specify using \link{standardization_control}.
#'
#' @return A list of class \code{smr}
#'
#' @export
direct_sub <- function(data, trt, outcome, baseline, outcome_type = "binomial", folds = 5, learners = c("mean", "glm"), control = standardization_control()) {
  task <- tsmr_Task$new(
    data = data,
    trt = trt,
    outcome = outcome,
    baseline = baseline,
    outcome_type = outcome_type,
    folds = folds
  )

  Qtilde <- outcome_regression(task, learners, control$.return_full_fits, control$.learners_outcome_folds, include_treatment = TRUE)

  theta_direct_sub(task, Qtilde$predicted_outcomes)
}

#' Probability Weighted Estimator for Direct Standardization
#'
#' @param data \[\code{data.frame}\]\cr
#' A \code{data.frame} containing all baseline, treatment, and outcome variables.
#' @param trt \[\code{character}\]\cr
#' Column name of treatment variable.
#' @param outcome \[\code{character}\]\cr
#' Column name of outcome variable.
#' @param baseline \[\code{character}\]\cr
#' Vector of column names of baseline variables.
#' @param outcome_type \[\code{character}\]\cr
#' Outcome variable type: binomial (binary) or continuous.
#' @param folds \[\code{integer(1)}\]\cr
#' Number of folds for top level cross-fitting.
#' @param learners \[\code{character}\]\cr
#' Vector of learners to include in SuperLearner library for estimating treatment assignment mechanism.
#' @param control \[\code{standardization_control}\]\cr
#' Additional tuning parameters for controlling fitting. Specify using \link{standardization_control}.
#'
#' @return A list of class \code{smr}
#'
#' @export
direct_pw <- function(data, trt, outcome, baseline, outcome_type = "binomial", trt_method = "default", folds = 5, learners = c("mean", "glm"), control = standardization_control()) {
  task <- tsmr_Task$new(
    data = data,
    trt = trt,
    outcome = outcome,
    baseline = baseline,
    outcome_type = outcome_type,
    folds = folds
  )

  g <- NULL
  riesz <- NULL
  if(trt_method == "default") {
    g <- treatment_probability(task, learners, control$.return_full_fits, control$.learners_trt_folds)
  }
  else if(tolower(trt_method) == "superriesz") {
    riesz <- riesz_representer(task, learners, control$.return_full_fits, control$.learners_trt_folds, parameter = "direct", method = "superriesz")
  }
  else if(tolower(trt_method) == "torch") {
    riesz <- riesz_representer(task, learners, control$.return_full_fits, control$.learners_trt_folds, parameter = "direct", method = "torch", torch_params = control$.torch_params)
  }

  theta <- theta_direct_pw(task, g, riesz)
  theta$g <- g
  theta$riesz <- riesz
  theta
}
