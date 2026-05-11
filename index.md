# FGstab

**Stability-based variable selection for the Fine-Gray competing risks
model**

FGstab selects relevant predictors from high-dimensional competing risks
or survival data, with statistical guarantees on the number of false
discoveries. It combines the penalised Fine-Gray model (Fu, Parikh &
Zhou, 2017) with Integrated Path Stability Selection (Melikechi &
Miller, 2026), extending the stability selection framework of
Meinshausen & Bühlmann (2010) to competing risks regression. Also
supports standard survival data with no competing event.

------------------------------------------------------------------------

## Installation

### 1. **IMPORTANT – Install `fastcmprsk` first**

This package **requires** `fastcmprsk`, which is not available on
CRAN.  
You must install it from GitHub **before** using FGstab:

``` r

# Install devtools if needed
install.packages("devtools")

# Install fastcmprsk from GitHub (MANDATORY)
devtools::install_github("erickawaguchi/fastcmprsk")
```

### 2. **Install `FGstab`**

``` r

devtools::install_github("hugocannafarina/FGstab")
```

------------------------------------------------------------------------

## Quick start

``` r

library(FGstab)
set.seed(123)

# Simulate data: 250 individuals, 10 true predictors, 500 total variables
dat <- sim_competing(n_ind = 250, n_pred = 10, n_var = 500,
                     design = "independent")
y <- cbind(time = dat$surv$Tobs, status = dat$surv$evtype)

# Select features controlling the FDR
res <- FGstab(dat$X, y, target_fdr = 0.1)
print(res)
plot(res)

# Select features controlling the expected number of false positives
res_fp <- FGstab(dat$X, y, target_fp = 1)
print(res_fp)

# Standard survival data (no competing event) — reduces to Cox
dat_cox <- sim_survival(n_ind = 250, n_pred = 10, n_var = 500,
                        design = "independent")
y_cox <- cbind(time = dat_cox$surv$Tobs, status = dat_cox$surv$evtype)
res_cox <- FGstab(dat_cox$X, y_cox, target_fdr = 0.1)
print(res_cox)
```

------------------------------------------------------------------------

## Package structure

    FGstab/
    ├── R/
    │   ├── stability_selection.R   # FGstab(), S3 class ipss_result, helpers
    │   ├── sim_data.R              # sim_competing(), sim_survival()
    │   └── plot_stability_paths.R  # plot.ipss_result()
    ├── tests/testthat/
    │   └── test-FGstab.R
    ├── DESCRIPTION
    ├── NAMESPACE
    └── README.md

------------------------------------------------------------------------

## Author

Hugo Cannafarina — <hugo.cannafarina@free.fr>

------------------------------------------------------------------------

## References

- Melikechi, O., & Miller, J. W. (2026). Integrated path stability
  selection. Journal of the American Statistical Association, 121(553),
  454–464.
- Meinshausen, N. & Bühlmann, P. (2010). Stability selection. *Journal
  of the Royal Statistical Society: Series B*, 72(4), 417–473.
- Fu, Z., Parikh, C. R. & Zhou, B. (2017). Penalized variable selection
  in competing risks regression. *Lifetime Data Analysis*, 23(3),
  353–376.
- Fine, J. P. & Gray, R. J. (1999). A proportional hazards model for the
  subdistribution of a competing risk. *JASA*, 94(446), 496–509.
