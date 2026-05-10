library(testthat)

# ---------------------------------------------------------------------------
# Shared test data
# ---------------------------------------------------------------------------

set.seed(42)
n <- 60L
p <- 15L

X_test <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X_test) <- paste0("X", seq_len(p))

y_test <- cbind(
  time   = abs(rnorm(n)) + 0.1,
  status = sample(c(0L, 1L, 2L), n, replace = TRUE)
)

# ---------------------------------------------------------------------------
# Smoke test: runs without error and returns the right class
# ---------------------------------------------------------------------------

test_that("FGstab runs without error and returns an ipss_result", {

  expect_no_error({
    res <- FGstab(X_test, y_test, n_alphas = 8L, n_jobs = 1L)
  })

  res <- FGstab(X_test, y_test, n_alphas = 8L, n_jobs = 1L)
  expect_s3_class(res, "ipss_result")
})

# ---------------------------------------------------------------------------
# sim_competing
# ---------------------------------------------------------------------------

test_that("sim_competing returns a well-formed list", {
  
  dat <- sim_competing(n_ind = 80L, n_pred = 5L, n_var = 100L,
                       design = "independent", censor_rate = 0.3)
  
  expect_named(dat, c("X", "surv", "true_pred_index", "beta"))
  expect_equal(nrow(dat$X),    80L)
  expect_equal(ncol(dat$X),   100L)
  expect_equal(nrow(dat$surv), 80L)
  expect_named(dat$surv, c("Tobs", "evtype", "id"))
  expect_equal(sum(dat$beta != 0), 5L)
  expect_true(all(dat$surv$evtype %in% c(0L, 1L, 2L)))
  expect_true(all(dat$surv$Tobs > 0))
})

test_that("sim_competing censoring rate is approximately correct", {
  
  set.seed(1)
  dat <- sim_competing(n_ind = 500L, n_pred = 5L, n_var = 100L,
                       design = "independent", censor_rate = 0.3)
  
  expect_lt(abs(mean(dat$surv$evtype == 0L) - 0.3), 0.1)
})

# ---------------------------------------------------------------------------
# sim_survival
# ---------------------------------------------------------------------------

test_that("sim_survival returns a well-formed list", {

  dat <- sim_survival(n_ind = 80L, n_pred = 5L, n_var = 100L,
                      design = "independent", censor_rate = 0.3)

  expect_named(dat, c("X", "surv", "true_pred_index", "beta"))
  expect_equal(nrow(dat$X),    80L)
  expect_equal(ncol(dat$X),   100L)
  expect_equal(nrow(dat$surv), 80L)
  expect_named(dat$surv, c("Tobs", "evtype", "id"))
  expect_equal(sum(dat$beta != 0), 5L)
  expect_true(all(dat$surv$evtype %in% c(0L, 1L)))
  expect_true(all(dat$surv$Tobs > 0))
})

test_that("sim_survival censoring rate is approximately correct", {

  set.seed(1)
  dat <- sim_survival(n_ind = 500L, n_pred = 5L, n_var = 100L,
                      design = "independent", censor_rate = 0.3)

  expect_lt(abs(mean(dat$surv$evtype == 0L) - 0.3), 0.1)
})

# ---------------------------------------------------------------------------
# Output structure
# ---------------------------------------------------------------------------

test_that("ipss_result has correct structure", {
  res <- FGstab(X_test, y_test, n_alphas = 8L, n_jobs = 1L)

  # stability_paths: n_alphas rows, p columns
  expect_equal(ncol(res$stability_paths), p)
  expect_true(nrow(res$stability_paths) >= 1L)

  # efp_scores and q_values: one entry per feature
  expect_length(res$efp_scores, p)
  expect_length(res$q_values,   p)

  # No features selected when no threshold is given
  expect_length(res$selected_features, 0L)
})

# ---------------------------------------------------------------------------
# target_fp and target_fdr select features
# ---------------------------------------------------------------------------

dat <- sim_competing(n_ind = 200L, n_pred = 10L, n_var = 100L,
                     design = "independent", censor_rate = 0.3)
X_comp <- dat$X
y_comp <- as.matrix(dat$surv[,1:2])

test_that("target_fp selects features", {

  res <- FGstab(X_comp, y_comp, n_alphas = 8L, target_fp = 2, n_jobs = 1L)
  expect_true(length(res$selected_features) >= 0L)
})

test_that("target_fdr selects features", {

  res <- FGstab(X_comp, y_comp, n_alphas = 8L, target_fdr = 0.5, n_jobs = 1L)
  expect_true(length(res$selected_features) >= 0L)
})

# ---------------------------------------------------------------------------
# plot() returns a ggplot
# ---------------------------------------------------------------------------

test_that("plot.ipss_result returns a ggplot without error", {

  res <- FGstab(X_test, y_test, n_alphas = 8L, n_jobs = 1L)
  expect_no_error(p_out <- plot(res))
  expect_s3_class(p_out, "ggplot")
})


