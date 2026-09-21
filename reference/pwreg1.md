# A test version of pwreg

This function is a test version of the `pwreg` function. It is intended
to simplify the code and make computations more efficient.

## Usage

``` r
pwreg1(
  id,
  time,
  status,
  Z,
  wfun = NULL,
  Y = NULL,
  strata = NULL,
  fixedL = TRUE,
  eps = 1e-06,
  maxiter = 50
)
```

## Arguments

- id:

  A vector of subject identifiers.

- time:

  A vector of event times.

- status:

  A vector of event indicators.

- Z:

  A matrix or data frame (tibble) of covariates.

- wfun:

  A function for the win function (Default is standard Pocock's win
  function).

- Y:

  An optional matrix or data frame (tibble) of continuous response
  variables.

- strata:

  A vector of strata indicators.

- fixedL:

  A logical value indicating whether the frailty variance is fixed.

- eps:

  A numeric value for the convergence criterion.

- maxiter:

  A numeric value for the maximum number of iterations.
