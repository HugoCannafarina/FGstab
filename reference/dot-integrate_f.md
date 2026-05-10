# Numerical integration along the regularisation path

Weighted integral of `values` w.r.t. `alphas` on a log scale. Stops
early if `cutoff` is exceeded.

## Usage

``` r
.integrate_f(values, alphas, delta = 1, cutoff = NULL)
```

## Arguments

- values:

  Numeric vector (length n_alphas).

- alphas:

  Numeric vector of lambda values (decreasing).

- delta:

  Numeric. Weight exponent (`1` = uniform on log scale).

- cutoff:

  Numeric or `NULL`. Early-stop threshold.

## Value

List with `output` (numeric scalar) and `stop_index` (integer).
