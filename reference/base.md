# Compute the baseline parameters needed for sample size calculation for standard win ratio test

Compute the baseline parameters \\\zeta_0^2\\ and
\\\boldsymbol\delta_0\\ needed for sample size calculation for standard
win ratio test (see
[`WRSS`](https://lmaowisc.github.io/WR/reference/WRSS.md)). The
calculation is based on a Gumbel–Hougaard copula model for survival time
\\D^{(a)}\\ and nonfatal event time \\T^{(a)}\\ for group \\a\\ (1:
treatment; 0: control): \$\${P}(D^{(a)}\>s, T^{(a)}\>t)
=\exp\left(-\left\[\left\\\exp(a\xi_1)\lambda_Ds\right\\^\kappa+
\left\\\exp(a\xi_2)\lambda_Ht\right\\^\kappa\right\]^{1/\kappa}\right),\$\$
where \\\xi_1\\ and \\\xi_2\\ are the component-wise log-hazard ratios
to be used as effect size in
[`WRSS`](https://lmaowisc.github.io/WR/reference/WRSS.md). We also
assume that patients are recruited uniformly over the period \\\[0,
\tau_b\]\\ and followed until time \\\tau\\ (\\\tau\geq\tau_b\\), with
an exponential loss-to-follow-up hazard \\\lambda_L\\.

## Usage

``` r
base(lambda_D, lambda_H, kappa, tau_b, tau, lambda_L, N = 1000, seed = 12345)
```

## Arguments

- lambda_D:

  Baseline hazard \\\lambda_D\\ for death.

- lambda_H:

  Baseline hazard \\\lambda_H\\ for nonfatal event.

- kappa:

  Gumbel–Hougaard copula correlation parameter \\\kappa\\.

- tau_b:

  Length of the initial (uniform) accrual period \\\tau_b\\.

- tau:

  Total length of follow-up \\\tau\\.

- lambda_L:

  Exponential hazard rate \\\lambda_L\\ for random loss to follow-up.

- N:

  Simulated sample size for monte-carlo integration.

- seed:

  Seed for monte-carlo simulation.

## Value

A list containing real number `zeta2` for \\\zeta_0^2\\ and bivariate
vector `delta` for \\\boldsymbol\delta_0\\.

## References

Mao, L., Kim, K. and Miao, X. (2021). Sample size formula for general
win ratio analysis. Biometrics, https://doi.org/10.1111/biom.13501.

## See also

[`gumbel.est`](https://lmaowisc.github.io/WR/reference/gumbel.est.md),
[`WRSS`](https://lmaowisc.github.io/WR/reference/WRSS.md)

## Examples

``` r
# see the example for WRSS
```
