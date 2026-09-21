# Compute the sample size for standard win ratio test

Compute the sample size for standard win ratio test.

## Usage

``` r
WRSS(xi, bparam, q = 0.5, alpha = 0.05, side = 2, power = 0.8)
```

## Arguments

- xi:

  A bivariate vector of hypothesized component-wise
  (treatment-to-control) log-hazard ratios under the Gumbel–Hougaard
  copula model described in
  [`base`](https://lmaowisc.github.io/WR/reference/base.md).

- bparam:

  A list containing baseline parameters `zeta2` for \\\zeta_0^2\\ and
  `delta` for \\\boldsymbol\delta_0\\; Can directly use the output of
  [base](https://lmaowisc.github.io/WR/reference/base.md).

- q:

  Proportion of patients assigned to treatment.

- alpha:

  Type I error rate.

- side:

  2-sided or 1-sided test.

- power:

  Target power.

## Value

A list containing `n`, the computed sample size.

## References

Mao, L., Kim, K. and Miao, X. (2021). Sample size formula for general
win ratio analysis. Biometrics, https://doi.org/10.1111/biom.13501.

## See also

[`gumbel.est`](https://lmaowisc.github.io/WR/reference/gumbel.est.md),
[`base`](https://lmaowisc.github.io/WR/reference/base.md)

## Examples

``` r
# The following is not run in package checking to save time.
if (FALSE) { # \dontrun{
## load the package and pilot dataset
library(WR)
head(hfaction_cpx9)
dat<-hfaction_cpx9
## subset to control group
pilot<-dat[dat$trt_ab==0,]

## get the data ready for gumbel.est()
id<-pilot$patid
## convert time from month to year
time<-pilot$time/12
status<-pilot$status
## compute the baseline parameters for the Gumbel--Hougaard
## copula for death and hospitalization
gum<-gumbel.est(id, time, status)

## get the baseline parameters
lambda_D<-gum$lambda_D
lambda_H<-gum$lambda_H
kappa<-gum$kappa
## set up design parameters and use base()
## to calculate bparam for WRSS()
# max follow-up 4 years
tau<-4
# 3 years of initial accrual
tau_b<-3
# loss to follow-up rate
lambda_L=0.05
# compute the baseline parameters
bparam<-base(lambda_D,lambda_H,kappa,tau_b,tau,lambda_L)
bparam

## sample size with power=0.8 under hazard ratios
## 0.9 and 0.8 for death and hospitalization, respectively.
WRSS(xi=log(c(0.9,0.8)),bparam=bparam,q=0.5,alpha=0.05,
    power=0.8)$n
## sample size under the same set-up but with power 0.9
WRSS(xi=log(c(0.9,0.8)),bparam=bparam,q=0.5,alpha=0.05,
    power=0.9)$n
} # }
```
