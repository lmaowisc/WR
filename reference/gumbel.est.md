# Estimate baseline parameters in the Gumbel–Hougaard model for sample size calculation using pilot data

Estimate baseline parameters in the Gumbel–Hougaard model described in
[`base`](https://lmaowisc.github.io/WR/reference/base.md) for sample
size calculation using pilot study data.

## Usage

``` r
gumbel.est(id, time, status)
```

## Arguments

- id:

  A vector of unique patient identifiers.

- time:

  A numeric vector of event times.

- status:

  A vector of event type variable; 2 = nonfatal event, 1 = death, and 0
  = censoring.

## Value

A list containing `lambda_D` for \\\lambda_D\\, `lambda_H` for
\\\lambda_H\\, and `kappa` for \\\kappa\\ in the Gumbel–Hougaard model.

## References

Mao, L., Kim, K. and Miao, X. (2021). Sample size formula for general
win ratio analysis. Biometrics, https://doi.org/10.1111/biom.13501.

## See also

[`base`](https://lmaowisc.github.io/WR/reference/base.md),
[`WRSS`](https://lmaowisc.github.io/WR/reference/WRSS.md)

## Examples

``` r
# see the example for WRSS
```
