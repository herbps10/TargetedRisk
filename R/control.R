#' Algorithm tuning parameters for estimating standardized mortality ratios
#'
#' @param .return_full_fits \[\code{logical(1)}\]\cr
#' Whether to return full SuperLearner fits
#' @param .learners_trt_folds \[\code{integer(1)}\]\cr
#' Number of cross-fitting folds to use in SuperLearner for estimating treatment assignment mechanism
#' @param .learners_outcome_folds \[\code{integer(1)}\]\cr
#' Number of cross-fitting folds to use in SuperLearner for estimating outcome regression.
#'
#' @export
smr_control = function(.return_full_fits = FALSE, .learners_trt_folds = 5, .learners_outcome_folds = 5) {
  list(
    .return_full_fits = .return_full_fits,
    .learners_trt_folds = .learners_trt_folds,
    .learners_outcome_folds = .learners_outcome_folds
  )
}
