# Binary search for lambda_max

Finds the smallest lambda such that no feature is selected.

## Usage

``` r
.get_alpha_max(X, y, penalty.factor, eps = 0.01, iter_max = 40L)
```

## Arguments

- X:

  Numeric matrix (n x p).

- y:

  Two-column matrix (time, status).

- penalty.factor:

  Numeric vector (length p).

- eps:

  Numeric. Bisection tolerance. Default `0.01`.

- iter_max:

  Integer. Maximum iterations. Default `40`.

## Value

Numeric scalar: lambda_max.
