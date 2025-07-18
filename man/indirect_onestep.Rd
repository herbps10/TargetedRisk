% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimators_indirect.R
\name{indirect_onestep}
\alias{indirect_onestep}
\title{One-step Estimator for Indirect Standardization}
\usage{
indirect_onestep(
  data,
  trt,
  outcome,
  baseline,
  outcome_type = c("binomial"),
  folds = 5,
  learners_trt = c("mean", "glm"),
  learners_outcome = c("mean", "glm"),
  Qtilde = NULL,
  g = NULL,
  verbose = FALSE,
  control = standardization_control()
)
}
\arguments{
\item{data}{\[\code{data.frame}\]\cr
A \code{data.frame} containing all baseline, treatment, and outcome variables.}

\item{trt}{\[\code{character}\]\cr
Column name of treatment variable.}

\item{outcome}{\[\code{character}\]\cr
Column name of outcome variable.}

\item{baseline}{\[\code{character}\]\cr
Vector of column names of baseline variables.}

\item{outcome_type}{\[\code{character}\]\cr
Outcome variable type: binomial (binary) or continuous.}

\item{folds}{\[\code{integer(1)}\]\cr
Number of folds for top level cross-fitting.}

\item{learners_trt}{\[\code{character}\]\cr
Vector of learners to include in SuperLearner library for estimating treatment assignment mechanism.}

\item{learners_outcome}{\[\code{character}\]\cr
Vector of learners to include in SuperLearner library for estimating outcome regression.}

\item{Qtilde}{\[\code{matrix}\]\cr
Optional precomputed Qtilde matrix}

\item{g}{\[\code{matrix}\]\cr
Optional precomputed probability of treatment matrix}

\item{verbose}{\[\code{logical}]\cr
Whether to print information messages during fitting}

\item{control}{\[\code{standardization_control}\]\cr
Additional tuning parameters for controlling fitting. Specify using \link{standardization_control}.}
}
\value{
A list of class \code{smr}
}
\description{
One-step Estimator for Indirect Standardization
}
