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
indirect_tmle <- function(data, trt, outcome, baseline, outcome_type = c("binomial"), folds = 5, learners_trt = c("mean", "glm"), learners_outcome = c("mean", "glm"), Qtilde = NULL, g = NULL, verbose = FALSE, control = standardization_control()) {
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

  if(verbose == TRUE) cat("Starting treatment fitting\n")
  if(is.null(g)) {
    g <- treatment_probability(task, learners_trt, control$.return_full_fits, control$.learners_trt_folds, verbose)
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

  fluctuations <- tmle_indirect(task, ybar, trt_prop, Qtilde, g)

  theta <- theta_indirect_tmle(task, trt_prop, fluctuations, g)
  theta$g <- g
  theta$Qtilde <- Qtilde

  theta
}

#' One-step Estimator for Indirect Standardization
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
indirect_onestep <- function(data, trt, outcome, baseline, outcome_type = c("binomial"), folds = 5, learners_trt = c("mean", "glm"), learners_outcome = c("mean", "glm"), Qtilde = NULL, g = NULL, verbose = FALSE, control = standardization_control()) {
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

  if(verbose == TRUE) cat("Starting treatment fitting\n")
  if(is.null(g)) {
    g <- treatment_probability(task, learners_trt, control$.return_full_fits, control$.learners_trt_folds, verbose)
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

  theta <- theta_indirect_onestep(task, trt_prop, ybar, Qtilde$predicted_outcomes, g)
  theta$g <- g
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
indirect_pw <- function(data, trt, outcome, baseline, outcome_type = "binomial", folds = 5, learners = c("mean", "glm"), control = standardization_control()) {
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

  g <- treatment_probability(task, learners, control$.return_full_fits, control$.learners_trt_folds)

  theta <- theta_indirect_pw(task, ybar, trt_prop, g)
  theta$g <- g
  theta
}
