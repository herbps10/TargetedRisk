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
#' @param Qtilde \[\code{matrix}\]\cr
#' Optional precomputed Qtilde matrix
#' @param g \[\code{matrix}\]\cr
#' Optional precomputed probability of treatment matrix
#' @param verbose \[\code{logical}]\cr
#' Whether to print information messages during fitting
#' @param control \[\code{standardization_control}\]\cr
#' Additional tuning parameters for controlling fitting. Specify using \link{standardization_control}.
#'
#' @return A list of class \code{smr}
#'
#' @export
direct_tmle <- function(data, trt, outcome, baseline, outcome_type = "binomial", folds = 5, learners_trt = c("mean", "glm"), learners_outcome = c("mean", "glm"), Qtilde = NULL, g = NULL, verbose = FALSE, control = standardization_control()) {
  if(length(outcome_type) > 1) outcome_type <- outcome_type[1]

  task <- tsmr_Task$new(
    data = data,
    trt = trt,
    outcome = outcome,
    baseline = baseline,
    outcome_type = outcome_type,
    folds = folds
  )

  if(verbose == TRUE) cat("Starting treatment fitting\n")
  if(is.null(g)) {
    g <- treatment_probability(task, learners_trt, control$.return_full_fits, control$.learners_trt_folds, verbose)
  }
  else {
    g <- list(treatment_probs = g)
  }

  if(verbose == TRUE) cat("Starting outcome fitting\n")
  if(is.null(Qtilde)) {
    Qtilde <- outcome_regression(task, learners_outcome, include_treatment = TRUE, control$.return_full_fits, control$.learners_outcome_folds)
  }
  else {
    Qtilde <- list(predicted_outcomes = Qtilde)
  }

  fluctuations <- tmle_direct(task, Qtilde, g)

  theta <- theta_direct_tmle(task, fluctuations, g)
  theta$g <- g
  theta$Qtilde <- Qtilde

  theta
}

#' One-step Estimator for Direct Standardization
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
#' @param Qtilde \[\code{matrix}\]\cr
#' Optional precomputed Qtilde matrix
#' @param g \[\code{matrix}\]\cr
#' Optional precomputed probability of treatment matrix
#' @param verbose \[\code{logical}]\cr
#' Whether to print information messages during fitting
#' @param control \[\code{standardization_control}\]\cr
#' Additional tuning parameters for controlling fitting. Specify using \link{standardization_control}.
#'
#' @return A list of class \code{smr}
#'
#' @export
direct_onestep <- function(data, trt, outcome, baseline, outcome_type = "binomial", folds = 5, learners_trt = c("mean", "glm"), learners_outcome = c("mean", "glm"), Qtilde = NULL, g = NULL, verbose = FALSE, control = standardization_control()) {
  if(length(outcome_type) > 1) outcome_type <- outcome_type[1]

  task <- tsmr_Task$new(
    data = data,
    trt = trt,
    outcome = outcome,
    baseline = baseline,
    outcome_type = outcome_type,
    folds = folds
  )

  if(verbose == TRUE) cat("Starting treatment fitting\n")
  if(is.null(g)) {
    g <- treatment_probability(task, learners_trt, control$.return_full_fits, control$.learners_trt_folds, verbose)
  }
  else {
    g <- list(treatment_probs = g)
  }

  if(verbose == TRUE) cat("Starting outcome fitting\n")
  if(is.null(Qtilde)) {
    Qtilde <- outcome_regression(task, learners_outcome, include_treatment = TRUE, control$.return_full_fits, control$.learners_outcome_folds)
  }
  else {
    Qtilde <- list(predicted_outcomes = Qtilde)
  }

  theta <- theta_direct_onestep(task, Qtilde$predicted_outcomes, g)
  theta$g <- g
  theta$Qtilde <- Qtilde

  theta
}

#' WeightIt Estimator for Direct Standardization
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
#' @param method \[\code{character}\]\cr
#' Method used within WeightIt for weights estimation. See WeightIt documentation for available methods.
#' @param verbose \[\code{logical}]\cr
#' Whether to print information messages during fitting
#' @param control \[\code{standardization_control}\]\cr
#' Additional tuning parameters for controlling fitting. Specify using \link{standardization_control}.
#'
#' @return A list of class \code{smr}
#'
#' @export
direct_weightit <- function(data, trt, outcome, baseline, outcome_type = "binomial", method = "glm", verbose = FALSE, control = standardization_control()) {
  if(length(outcome_type) > 1) outcome_type <- outcome_type[1]

  task <- tsmr_Task$new(
    data = data,
    trt = trt,
    outcome = outcome,
    baseline = baseline,
    outcome_type = outcome_type,
    folds = 1
  )

  weights <- balancing_weights(task, method, estimand = "ATE")
  theta <- theta_direct_weightit(task, weights)

  theta$weights <- weights

  theta
}

#' MatchIt Estimator for Direct Standardization
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
#' @param method \[\code{character}\]\cr
#' Method used within MatchIt to select matches. See MatchIt documentation for available methods.
#' @param distance \[\code{character}\]\cr
#' Distance method used by MatchIt. See MatchIt documentation for available distance methods.
#' @param verbose \[\code{logical}]\cr
#' Whether to print information messages during fitting
#' @param control \[\code{standardization_control}\]\cr
#' Additional tuning parameters for controlling fitting. Specify using \link{standardization_control}.
#'
#' @return A list of class \code{smr}
#'
#' @export
direct_matchit <- function(data, trt, outcome, baseline, outcome_type = "binomial", method = "nearest", distance = "glm", verbose = FALSE, control = standardization_control()) {
  if(length(outcome_type) > 1) outcome_type <- outcome_type[1]

  task <- tsmr_Task$new(
    data = data,
    trt = trt,
    outcome = outcome,
    baseline = baseline,
    outcome_type = outcome_type,
    folds = 1
  )

  matches <- matching(task, method, distance, estimand = "ATE")
  theta <- theta_direct_matchit(task, matches)

  theta$matches <- matches

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
direct_pw <- function(data, trt, outcome, baseline, outcome_type = "binomial", folds = 5, learners = c("mean", "glm"), control = standardization_control()) {
  task <- tsmr_Task$new(
    data = data,
    trt = trt,
    outcome = outcome,
    baseline = baseline,
    outcome_type = outcome_type,
    folds = folds
  )

  g <- treatment_probability(task, learners, control$.return_full_fits, control$.learners_trt_folds)

  theta <- theta_direct_pw(task, g)
  theta$g <- g
  theta
}
