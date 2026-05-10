# Fit LASSO Fine-Gray and return binary selection matrix

Fit LASSO Fine-Gray and return binary selection matrix

## Usage

``` r
.fit_lasso_fg(X, y, alphas)
```

## Arguments

- X:

  Numeric matrix (n_sub x p).

- y:

  Two-column matrix (time, status).

- alphas:

  Numeric vector of lambda values.

## Value

Integer matrix (n_alphas x p): 1 = selected, 0 = not, NA =
non-convergence.
