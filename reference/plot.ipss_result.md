# Plot stability paths from an ipss_result object

Displays the selection frequency of each feature along the
regularisation path. Features that are ultimately selected are
highlighted in colour; all others are drawn in light grey. The top-k
most stable features (by maximum selection frequency) are labelled on
the right margin.

## Usage

``` r
# S3 method for class 'ipss_result'
plot(
  x,
  n_highlight = 5L,
  selected_color = "#E63946",
  alpha_unselected = 0.4,
  ...
)
```

## Arguments

- x:

  An `ipss_result` object returned by
  [`FGstab`](https://hugocannafarina.github.io/FGstab/reference/FGstab.md).

- n_highlight:

  Integer. Number of top features to label on the right margin, ranked
  by maximum selection frequency. Default `5`.

- selected_color:

  Character. Colour used for selected feature paths. Default `"#E63946"`
  (red).

- alpha_unselected:

  Numeric in `[0, 1]`. Transparency of unselected feature paths. Default
  `0.4`.

- ...:

  Additional arguments (ignored).

## Value

A `ggplot` object (invisible). The plot is also printed.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- FGstab(X, y, n_alphas = 15, target_fdr = 0.2)
plot(res)
plot(res, n_highlight = 10)
} # }
```
