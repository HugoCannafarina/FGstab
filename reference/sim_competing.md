# Simulate competing risks data under the Fine-Gray model

Generates a synthetic dataset for a competing risks model with two event
types (event of interest and competing event), following the
sub-distribution hazard model of Fine and Gray (1999). Covariates \\X
\in \mathbb{R}^{n \times p}\\ are drawn from a multivariate normal
distribution \\X \sim \mathcal{N}(0, \Sigma)\\ whose covariance
structure \\\Sigma\\ is controlled by the `design` argument.

## Usage

``` r
sim_competing(
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
  structures* section for details.

- rho:

  Numeric in \\(0, 1)\\. Decay parameter for the `"toeplitz"` design.
  Ignored for other designs. Default `0.9`.

## Value

A named list with four elements:

- `X`:

  Numeric matrix of covariates (`n_ind` x `n_var`), centred and scaled.

- `surv`:

  Data frame with columns `Tobs` (observed time), `evtype` (0 =
  censored, 1 = event of interest, 2 = competing event), and `id`.

- `true_pred_index`:

  Integer vector of length `n_pred`: indices of the truly predictive
  covariates (\\\mathcal{S}\\).

- `beta`:

  Numeric vector of true coefficients (length `n_var`; non-zero only at
  `true_pred_index`).

## Covariance structures

- **Independent** (`"independent"`):

  \\\Sigma = I_p\\. Each covariate is drawn i.i.d. from \\\mathcal{N}(0,
  1)\\.

- **Toeplitz** (`"toeplitz"`):

  \$\$\Sigma\_{jk} = \rho^{\|j - k\|}, \quad j, k = 1, \ldots, p.\$\$
  Correlations decay exponentially with the distance between indices.
  The parameter `rho` \\\in (0, 1)\\ controls the decay rate: values
  close to 1 yield strongly correlated neighbouring covariates, values
  close to 0 approach the independent design.

- **Block** (`"block"`):

  \$\$ \Sigma\_{jk} = \begin{cases} 1 & \text{if } j = k, \\ 0.5 &
  \text{if } j \neq k \text{ and } \lfloor j / B \rfloor = \lfloor k / B
  \rfloor, \\ 0 & \text{otherwise,} \end{cases} \$\$ where \\B = 10\\ is
  the block size. Covariates within the same block share a fixed
  pairwise correlation of 0.5; covariates in different blocks are
  independent.

## Data-generating mechanism

True coefficients \\\beta \in \mathbb{R}^p\\ are sparse: `n_pred`
indices are drawn uniformly at random and assigned coefficients from
\\\mathcal{U}(-1, 1)\\; all other coefficients are zero. The marginal
probability of the event of interest for individual \\i\\ is \$\$ \pi_i
= 1 - (1 - \pi_0)^{\exp(X_i^\top \beta)}, \$\$ where \\\pi_0 = 0.75\\ is
the baseline event rate. Event times for the event of interest are drawn
via quantile inversion of the sub-distribution function; competing event
times are drawn from \\\mathrm{Exp}(\text{rate} = 3)\\. Censoring is
applied by thresholding at the `censor_rate` quantile of the observed
times.

## References

Fine, J. P. and Gray, R. J. (1999). A proportional hazards model for the
subdistribution of a competing risk. *Journal of the American
Statistical Association*, **94**(446), 496–509.

## Examples

``` r
# Independent covariates
dat <- sim_competing(n_ind = 100, n_pred = 5, n_var = 200, design = "independent")
str(dat)
#> List of 4
#>  $ X              : num [1:100, 1:200] -1.3996 0.1756 -2.3866 -0.0727 0.5241 ...
#>   ..- attr(*, "scaled:center")= num [1:200] 0.0708 0.1149 -0.03 0.1886 -0.1766 ...
#>   ..- attr(*, "scaled:scale")= num [1:200] 1.051 1.014 1.051 0.992 1.031 ...
#>  $ surv           :'data.frame': 100 obs. of  3 variables:
#>   ..$ Tobs  : num [1:100] 0.02171 0.2301 0.7918 0.00742 0.07869 ...
#>   ..$ evtype: int [1:100] 2 1 0 2 1 2 2 0 1 1 ...
#>   ..$ id    : int [1:100] 1 2 3 4 5 6 7 8 9 10 ...
#>  $ true_pred_index: int [1:5] 109 184 129 3 37
#>  $ beta           : num [1:200] 0 0 0.606 0 0 ...

# Toeplitz covariance (moderate correlation)
dat_toep <- sim_competing(n_ind = 150, n_pred = 8, n_var = 300,
                     design = "toeplitz", rho = 0.7)

# Block covariance
dat_block <- sim_competing(n_ind = 150, n_pred = 8, n_var = 300,
                      design = "block")
```
