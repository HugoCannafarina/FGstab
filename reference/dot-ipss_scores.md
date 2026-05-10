# Compute IPSS scores along the regularisation path

Compute IPSS scores along the regularisation path

## Usage

``` r
.ipss_scores(
  stability_paths,
  B,
  alphas,
  average_selected,
  ipss_function,
  delta,
  cutoff
)
```

## Arguments

- stability_paths:

  Matrix (n_alphas x p) of selection frequencies.

- B:

  Integer. Number of subsampling iterations.

- alphas:

  Numeric vector of regularisation values.

- average_selected:

  Numeric vector of average selected features per alpha.

- ipss_function:

  Character. One of `"h1"`, `"h2"`, `"h3"`.

- delta:

  Numeric. Integration weight exponent.

- cutoff:

  Numeric. Early-stop threshold for the integrated bound.

## Value

A list with `scores` (numeric vector, length p), `integral` (numeric
scalar), and `stop_index` (integer).
