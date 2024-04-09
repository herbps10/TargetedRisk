#' Targeted Minimum Loss-Based Estimator of the Standardized Mortality Ratio
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
#' @param control \[\code{smr_control}\]\cr
#' Additional tuning parameters for controlling fitting. Specify using \link{smr_control}.
#'
#' @return A list of class \code{smr}
#'
#' @export
smr_tmle <- function(data, trt, outcome, baseline, outcome_type = c("binomial", "continuous"), folds = 5, learners_trt = c("mean", "glm"), learners_outcome = c("mean", "glm"), control = smr_control()) {
  task <- tsmr_Task$new(
    data = data,
    trt = trt,
    outcome = outcome,
    baseline = baseline,
    outcome_type = outcome_type,
    folds = folds
  )

  if(length(outcome_type) > 1) outcome_type <- outcome_type[1]

  ybar <- mean_outcome(task)
  trt_prop <- treatment_proportion(task)
  g <- treatment_probability(task, learners_trt, control$.return_full_fits, control$.learners_trt_folds)
  Qtilde <- outcome_regression(task, learners_outcome, control$.return_full_fits, control$.learners_outcome_folds)
  fluctuations <- tmle(task, ybar, trt_prop, g, Qtilde)

  theta <- theta_tmle(task, trt_prop, g$treatment_probs, fluctuations)
  theta$g <- g

  theta
}

#' Substitution Estimator of the Standardized Mortality Ratio
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
#' @param control \[\code{smr_control}\]\cr
#' Additional tuning parameters for controlling fitting. Specify using \link{smr_control}.
#'
#' @return A list of class \code{smr}
#'
#' @export
smr_sub <- function(data, trt, outcome, baseline, outcome_type = "binomial", folds = 5, learners = c("mean", "glm"), control = smr_control()) {
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
  Qtilde <- outcome_regression(task, learners, control$.return_full_fits, control$.learners_outcome_folds)

  theta_sub(task, ybar, trt_prop, Qtilde$predicted_outcomes)
}

#' Probability Weighted Estimator of the Standardized Mortality Ratio
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
#' @param control \[\code{smr_control}\]\cr
#' Additional tuning parameters for controlling fitting. Specify using \link{smr_control}.
#'
#' @return A list of class \code{smr}
#'
#' @export
smr_pw <- function(data, trt, outcome, baseline, outcome_type = "binomial", folds = 5, learners = c("mean", "glm"), control = smr_control()) {
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

  theta_pw(task, ybar, trt_prop, g$treatment_probs)
}
