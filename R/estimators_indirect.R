#' Targeted Minimum Loss-Based Estimator for Indirect Standardization
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
#' @param trt_method \[\code{character}\]\cr
#' Method for estimating treatment assignment mechanism. One of "default", "SuperRiesz", or "torch".
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
indirect_tmle <- function(data, trt, outcome, baseline, outcome_type = c("binomial", "continuous"), folds = 5, trt_method = "default", learners_trt = c("mean", "glm"), learners_outcome = c("mean", "glm"), Qtilde = NULL, g = NULL, verbose = FALSE, control = standardization_control()) {
  if(length(outcome_type) > 1) outcome_type <- outcome_type[1]

  task <- tsmr_Task$new(
    data = data,
    trt = trt,
    outcome = outcome,
    baseline = baseline,
    outcome_type = outcome_type,
    folds = folds
  )

  ybar <- mean_outcome(task)
  trt_prop <- treatment_proportion(task)

  riesz <- NULL
  if(verbose == TRUE) cat("Starting treatment fitting\n")
  if(is.null(g)) {
    if(trt_method == "default") {
      g <- treatment_probability(task, learners_trt, control$.return_full_fits, control$.learners_trt_folds)
    }
    else if(tolower(trt_method) == "superriesz") {
      riesz <- riesz_representer(task, learners_trt, control$.return_full_fits, control$.learners_trt_folds, parameter = "smr", method = "superriesz")
    }
    else if(tolower(trt_method) == "torch") {
      riesz <- riesz_representer(task, learners_trt, control$.return_full_fits, control$.learners_trt_folds, parameter = "smr", method = "torch", torch_params = control$.torch_params)
    }
  }
  else {
    g <- list(treatment_probs = g)
  }

  if(verbose == TRUE) cat("Starting outcome fitting\n")
  if(is.null(Qtilde)) {
    Qtilde <- outcome_regression(task, learners_outcome, include_treatment = FALSE, control$.return_full_fits, control$.learners_outcome_folds)
  }
  else {
    Qtilde <- list(predicted_outcomes = Qtilde)
  }

  fluctuations <- tmle_indirect(task, ybar, trt_prop, Qtilde, g, riesz)

  theta <- theta_indirect_tmle(task, trt_prop, fluctuations, g, riesz)
  theta$g <- g
  theta$riesz <- riesz
  theta$Qtilde <- Qtilde

  theta
}

#' Substitution Estimator for Indirect Standardizatio
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
indirect_sub <- function(data, trt, outcome, baseline, outcome_type = "binomial", folds = 5, learners = c("mean", "glm"), control = standardization_control()) {
  task <- tsmr_Task$new(
    data = data,
    trt = trt,
    outcome = outcome,
    baseline = baseline,
    outcome_type = outcome_type,
    folds = folds
  )

  ybar <- mean_outcome(task)
  trt_prop <- treatment_proportion(task)
  Qtilde <- outcome_regression(task, learners, include_treatment = FALSE, control$.return_full_fits, control$.learners_outcome_folds)


  theta_indirect_sub(task, ybar, trt_prop, Qtilde$predicted_outcomes)
}

#' Probability Weighted Estimator for Indirect Standardizatio
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
#' @param trt_method \[\code{character}\]\cr
#' Method for estimating treatment assignment mechanism. One of "default", "SuperRiesz", or "torch".
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
indirect_pw <- function(data, trt, outcome, baseline, outcome_type = "binomial", trt_method = "default", folds = 5, learners = c("mean", "glm"), control = standardization_control()) {
  task <- tsmr_Task$new(
    data = data,
    trt = trt,
    outcome = outcome,
    baseline = baseline,
    outcome_type = outcome_type,
    folds = folds
  )

  ybar <- mean_outcome(task)
  trt_prop <- treatment_proportion(task)

  g <- NULL
  riesz <- NULL
  if(trt_method == "default") {
    g <- treatment_probability(task, learners, control$.return_full_fits, control$.learners_trt_folds)
  }
  else if(tolower(trt_method) == "superriesz") {
    riesz <- riesz_representer(task, learners, control$.return_full_fits, control$.learners_trt_folds, parameter = "smr", method = "superriesz")
  }
  else if(tolower(trt_method) == "torch") {
    riesz <- riesz_representer(task, learners, control$.return_full_fits, control$.learners_trt_folds, parameter = "smr", method = "torch", torch_params = control$.torch_params)
  }

  theta <- theta_indirect_pw(task, ybar, trt_prop, g, riesz)
  theta$g <- g
  theta$riesz <- riesz
  theta
}
