smr_tmle <- function(data, trt, outcome, baseline, outcome_type = "binomial", folds = 5, learners_trt = c("mean", "glm"), learners_outcome = c("mean", "glm"), control = smr_control()) {
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
  g <- treatment_probability(task, learners_trt, control$.return_full_fits, control$.learners_trt_folds)
  Qtilde <- outcome_regression(task, learners_outcome, control$.return_full_fits, control$.learners_outcome_folds)
  fluctuations <- tmle(task, ybar, trt_prop, g, Qtilde)

  theta <- theta_tmle(task, trt_prop, g$treatment_probs, fluctuations)
  theta$g <- g

  theta
}

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
