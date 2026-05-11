# FGstab

**Stability-based variable selection for the Fine-Gray competing risks
model**

FGstab provides reliable sparse variable selection for high-dimensional
competing risks and survival data. It implements Integrated Path
Stability Selection (IPSS) (Melikechi & Miller, 2026) for the penalised
Fine-Gray model, extending the stability selection framework of
Meinshausen & Bühlmann (2010) to competing risks regression. The method
aggregates selection frequencies along the LASSO regularisation path and
bounds the expected number of false positives and false discovery rate
via path integration. It also supports standard survival data with no
competing event, in which case it reduces to IPSS for the penalised Cox
model.

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
set.seed(1)

# Simulate data: 150 individuals, 8 true predictors, 300 total variables
dat <- sim_competing(n_ind = 150, n_pred = 8, n_var = 300,
                     design = "independent")
y <- cbind(time = dat$surv$Tobs, status = dat$surv$evtype)

# Run IPSS — returns scores without automatic selection
res <- FGstab(dat$X, y, n_alphas = 15)
print(res)
plot(res)

# Select features controlling the expected number of false positives
res_fp <- FGstab(dat$X, y, n_alphas = 15, target_fp = 1)
print(res_fp)

# Select features controlling the FDR
res_fdr <- FGstab(dat$X, y, n_alphas = 15, target_fdr = 0.2)
summary(res_fdr)
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
- Fine, J. P. & Gray, R. J. (1999). A proportional hazards model for the
  subdistribution of a competing risk. *JASA*, 94(446), 496–509.
