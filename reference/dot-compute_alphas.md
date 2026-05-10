# Compute the regularisation path

Finds lambda_max via bisection and builds a log-spaced grid. Truncates
at the first lambda where the model fails to converge.

## Usage

``` r
.compute_alphas(X, y, n_alphas, penalty.factor)
```

## Arguments

- X:

  Numeric matrix (n x p).

- y:

  Two-column matrix (time, status).

- n_alphas:

  Integer. Desired path length.

- penalty.factor:

  Numeric vector (length p).

## Value

Numeric vector of lambda values.
