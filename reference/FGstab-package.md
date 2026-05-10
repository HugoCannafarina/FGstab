# FGstab: Stability-Based Variable Selection for the Fine-Gray Competing Risks Model

Implements Integrated Path Stability Selection (IPSS) for penalised
Fine-Gray competing risks regression in high-dimensional settings. The
method repeatedly subsamples the data, fits a LASSO-penalised Fine-Gray
model across a regularisation path, and aggregates selection frequencies
into stability paths. EFP scores and q-values derived from path
integration allow control of the expected number of false positives or
the false discovery rate. The package also provides tools for simulating
competing risks data under several covariance structures and for
visualising stability paths.

## Author

**Maintainer**: Hugo Cannafarina <hugo.cannafarina@free.fr>
