superlearner <- function(data, outcome, learners, learner_folds = 5, outcome_type = "binomial") {
  mlr3superlearner::mlr3superlearner(
    data, outcome, learners, outcome_type, discrete = FALSE, folds = learner_folds
  )
}
