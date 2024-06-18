# `TargetedSMR`
`TargetedSMR` is an `R` package that provides methods for Targeted estimation of indirectly standardized mortality ratios (SMRs).

The interface of the package is similar to that of the `lmtp` package for longitudinal modified treatment policies.

# Installation

To install the development version from Github:
```
remotes::install_github("herbps10/TargetedSMR")
```

# Usage

Simulate test data:
```r
data <- simulate_providers(N = 1e3, providers = 5, covariates = 5)
```

Estimate direct standardization parameters with TMLE:
```r
direct_results <- direct_tmle(
  data,
  outcome = "Y",
  trt = "A",
  baseline = c("W1", "W2", "W3", "W4", "W5"),
  learners_trt = c("mean", "glm"),
  learners_outcome = c("mean", "glm")
)
```

Estimate indirect standardization parameters with TMLE:
```r
indirect_results <- indirect_tmle(
  data,
  outcome = "Y",
  trt = "A",
  baseline = c("W1", "W2", "W3", "W4", "W5"),
  learners_trt = c("mean", "glm"),
  learners_outcome = c("mean", "glm")
)
```
