#' @importFrom origami make_folds
cv_setup <- function(data, id, folds) {
  cv <- make_folds(data, cluster_ids = id, V = folds)
  if(folds == 1) {
    cv[[1]]$training_set <- cv[[1]]$validation_set
  }
  cv
}

extract_fold <- function(x, cv, cv_index) {
  list(
    training   = x[cv[[cv_index]]$training_set, , drop = FALSE],
    validation = x[cv[[cv_index]]$validation_set, , drop = FALSE]
  )
}

