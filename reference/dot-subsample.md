# One subsampling iteration

Splits data in half and fits LASSO Fine-Gray on each half.

## Usage

``` r
.subsample(X, y, alphas)
```

## Arguments

- X:

  Numeric matrix (n x p).

- y:

  Two-column matrix (time, status).

- alphas:

  Numeric vector of regularisation values.

## Value

Array (n_alphas x 2 x p) of binary selection indicators. `NA` for
non-converged fits.
