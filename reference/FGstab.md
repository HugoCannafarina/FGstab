# Variable Selection via IPSS for the Penalised Fine-Gray Model

Implements Integrated Path Stability Selection (IPSS) (Melikechi &
Miller, 2026) for the Fine-Gray competing-risks regression model. The
method repeatedly subsamples the data, fits a LASSO-penalised Fine-Gray
model on each subsample across a regularisation path, and aggregates
selection frequencies into stability paths. EFP scores and q-values are
derived from the integrated paths to control either the expected number
of false positives (`target_fp`) or the false discovery rate
(`target_fdr`). When `y[, 2]` contains only 0 and 1 (no competing
event), the Fine-Gray model reduces to the Cox model; `FGstab` then
performs IPSS for penalised Cox regression. See
[`sim_survival`](https://hugocannafarina.github.io/FGstab/reference/sim_survival.md)
to simulate data in that setting.

## Usage

``` r
FGstab(
  X,
  y,
  n_alphas = 15L,
  B = 50L,
  target_fp = NULL,
  target_fdr = 0.1,
  ipss_function = "h2",
  n_jobs = detectCores() - 1L
)
```

## Arguments

- X:

  Numeric matrix of covariates, dimensions n x p. Used as-is (no
  internal standardisation).

- y:

  Two-column numeric matrix. The first column must contain the event
  time (observed follow-up time). The second column must contain the
  event type indicator: 0 = censored, 1 = event of interest, 2 =
  competing event (omit 2 for a standard survival / Cox setting).

- n_alphas:

  Integer. Number of regularisation values (lambda) on the path. Default
  `15`.

- B:

  Integer. Number of subsampling iterations. Default `50`.

- target_fp:

  Numeric or `NULL`. Maximum expected number of false positives.
  Features with EFP score \\\leq\\ `target_fp` are selected. Ignored if
  `target_fdr` is also provided. Default `NULL`.

- target_fdr:

  Numeric in \\(0, 1)\\ or `NULL`. Maximum false discovery rate.
  Features with q-value \\\leq\\ `target_fdr` are selected. Takes
  priority over `target_fp`. Default 0.1.

- ipss_function:

  Character string. IPSS bounding function controlling conservatism. One
  of `"h1"`, `"h2"` (default), or `"h3"`. `"h3"` is the least
  conservative.

- n_jobs:

  Integer. Number of cores for parallel subsampling. Defaults to
  `detectCores() - 1`.

## Value

An object of class `ipss_result` with fields:

- `selected_features`:

  Integer vector of selected feature indices. Empty if neither
  `target_fp` nor `target_fdr` is specified.

- `stability_paths`:

  Matrix (n_alphas x p) of selection frequencies along the
  regularisation path.

- `efp_scores`:

  Named list of EFP scores. Lower values indicate more reliably selected
  features.

- `q_values`:

  Named numeric vector of q-values derived from EFP scores.
  Interpretable as an upper bound on the FDR.

## References

- Melikechi, O., & Miller, J. W. (2026). Integrated path stability
  selection. *Journal of the American Statistical Association*,
  121(553), 454–464.

- Fine, J. P. & Gray, R. J. (1999). A proportional hazards model for the
  subdistribution of a competing risk. *JASA*, 94(446), 496–509.

## See also

[`print.ipss_result`](https://hugocannafarina.github.io/FGstab/reference/print.ipss_result.md),
[`summary.ipss_result`](https://hugocannafarina.github.io/FGstab/reference/summary.ipss_result.md),
[`plot.ipss_result`](https://hugocannafarina.github.io/FGstab/reference/plot.ipss_result.md),
[`sim_competing`](https://hugocannafarina.github.io/FGstab/reference/sim_competing.md),
[`sim_survival`](https://hugocannafarina.github.io/FGstab/reference/sim_survival.md)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
dat <- sim_competing(n_ind = 150, n_pred = 10, n_var = 500,
                     design = "independent")
y <- as.matrix(dat$surv[, c("Tobs", "evtype")])

# Return scores only (inspect before thresholding)
res <- FGstab(dat$X, y, n_alphas = 15)
print(res)
plot(res)

# Control expected false positives
res_fp <- FGstab(dat$X, y, n_alphas = 15, target_fp = 1)
print(res_fp)

# Control FDR
res_fdr <- FGstab(dat$X, y, n_alphas = 15, target_fdr = 0.2)
print(res_fdr)
} # }

```
