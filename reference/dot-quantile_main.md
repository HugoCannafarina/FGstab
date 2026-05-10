# Simulate event times for the event of interest via quantile inversion

Draws event times from the sub-distribution of the event of interest
using the inverse-CDF method under the Fine-Gray model.

## Usage

``` r
.quantile_main(X_main, beta, p_main, main_evt_rate)
```

## Arguments

- X_main:

  Numeric matrix of covariates for individuals with the event of
  interest.

- beta:

  Numeric vector of true coefficients (length p).

- p_main:

  Numeric vector of marginal probabilities of the main event for the
  same individuals.

- main_evt_rate:

  Numeric. Overall rate of the main event in the population.

## Value

Numeric vector of simulated event times.
