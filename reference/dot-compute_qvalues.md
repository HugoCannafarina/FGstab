# Compute q-values from EFP scores

For each feature, the q-value is the minimum FDR achievable when
including that feature in the selected set.

## Usage

``` r
.compute_qvalues(efp_scores)
```

## Arguments

- efp_scores:

  Named list of EFP scores (one per feature).

## Value

Named numeric vector of q-values.
