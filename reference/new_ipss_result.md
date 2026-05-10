# Constructor for ipss_result objects

Constructor for ipss_result objects

## Usage

``` r
new_ipss_result(selected_features, stability_paths, efp_scores, q_values)
```

## Arguments

- selected_features:

  Integer vector of selected feature indices.

- stability_paths:

  Matrix of stability paths (n_alphas x p).

- efp_scores:

  Named list of EFP scores (one per feature).

- q_values:

  Named numeric vector of q-values derived from EFP scores.

## Value

An object of class `ipss_result`.
