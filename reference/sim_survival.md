# Simulate survival data under the Cox proportional hazards model

Generates a synthetic dataset for a standard survival model with no
competing risks, following the Cox proportional hazards model with an
exponential baseline hazard. Covariates \\X \in \mathbb{R}^{n \times
p}\\ are drawn from a multivariate normal distribution \\X \sim
\mathcal{N}(0, \Sigma)\\ whose covariance structure \\\Sigma\\ is
controlled by the `design` argument (same options as
[`sim_competing`](https://hugocannafarina.github.io/FGstab/reference/sim_competing.md)).

## Usage

``` r
sim_survival(
  n_ind = 250L,
  n_pred = 10L,
  n_var = 5000L,
  censor_rate = 0.3,
  design = c("independent", "toeplitz", "block"),
  rho = 0.9
)
```

## Arguments

- n_ind:

  Integer. Number of individuals. Default `250`.

- n_pred:

  Integer. Number of truly predictive covariates (\\\|\mathcal{S}\|\\).
  Default `10`.

- n_var:

  Integer. Total number of covariates (\\p\\). Default `5000`.

- censor_rate:

  Numeric in \\\[0, 1)\\. Target proportion of censored observations.
  Default `0.3`.

- design:

  Character string. Covariance structure of \\\Sigma\\. One of
  `"independent"`, `"toeplitz"`, or `"block"`. See the *Covariance
  structures* section of
  [`sim_competing`](https://hugocannafarina.github.io/FGstab/reference/sim_competing.md)
  for details.

- rho:

  Numeric in \\(0, 1)\\. Decay parameter for the `"toeplitz"` design.
  Ignored for other designs. Default `0.9`.

## Value

A named list with four elements:

- `X`:

  Numeric matrix of covariates (`n_ind` x `n_var`), centred and scaled.

- `surv`:

  Data frame with columns `Tobs` (observed time), `evtype` (0 =
  censored, 1 = event), and `id`.

- `true_pred_index`:

  Integer vector of length `n_pred`: indices of the truly predictive
  covariates (\\\mathcal{S}\\).

- `beta`:

  Numeric vector of true coefficients (length `n_var`; non-zero only at
  `true_pred_index`).

## Data-generating mechanism

True coefficients \\\beta \in \mathbb{R}^p\\ are sparse: `n_pred`
indices are drawn uniformly at random and assigned coefficients from
\\\mathcal{U}(-1, 1)\\; all other coefficients are zero. Event times are
drawn from an exponential distribution with rate \\\exp(X_i^\top
\beta)\\: \$\$T_i \sim \mathrm{Exp}\\\left(\exp(X_i^\top
\beta)\right).\$\$ This corresponds to a Cox model with unit exponential
baseline hazard. Censoring is applied by thresholding at the
`censor_rate` quantile of the observed times.

The output is formatted identically to
[`sim_competing`](https://hugocannafarina.github.io/FGstab/reference/sim_competing.md)
with `evtype` \\\in \\0, 1\\\\ (no competing event), so it can be passed
directly to
[`FGstab`](https://hugocannafarina.github.io/FGstab/reference/FGstab.md),
which then operates as IPSS for the Cox model.

## See also

[`sim_competing`](https://hugocannafarina.github.io/FGstab/reference/sim_competing.md),
[`FGstab`](https://hugocannafarina.github.io/FGstab/reference/FGstab.md)

## Examples

``` r
dat <- sim_survival(n_ind = 100, n_pred = 5, n_var = 200,
                    design = "independent")
str(dat)
#> List of 4
#>  $ X              : num [1:100, 1:200] 0.4636 0.352 -0.0817 -0.1719 1.0143 ...
#>   ..- attr(*, "scaled:center")= num [1:200] 0.1684 0.1789 0.063 -0.0212 0.1375 ...
#>   ..- attr(*, "scaled:scale")= num [1:200] 1.033 0.895 0.991 0.932 1.044 ...
#>  $ surv           :'data.frame': 100 obs. of  3 variables:
#>   ..$ Tobs  : num [1:100] 0.334 13.755 11.546 0.616 0.273 ...
#>   ..$ evtype: int [1:100] 1 0 0 1 1 0 1 0 1 1 ...
#>   ..$ id    : int [1:100] 1 2 3 4 5 6 7 8 9 10 ...
#>  $ true_pred_index: int [1:5] 127 159 182 37 1
#>  $ beta           : num [1:200] 0.316 0 0 0 0 ...

if (FALSE) { # \dontrun{
y <- as.matrix(dat$surv[, c("Tobs", "evtype")])
res <- FGstab(dat$X, y, n_alphas = 10, target_fdr = 0.2)
print(res)
} # }
```
