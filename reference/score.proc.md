# Computes the standardized score processes

Computes the standarized score processes for the covariates.

## Usage

``` r
score.proc(obj, t = NULL)
```

## Arguments

- obj:

  an object of class pwreg.

- t:

  a vector containing times. If not specified, the function will use all
  unique event times from the data.

## Value

An object of class `pwreg.score` consisting of `t:` a vector of times;
and `score:` a matrix whose rows are the standardized score processes as
a function of `t`.

## References

Mao, L. and Wang, T. (2020). A class of proportional win-fractions
regression models for composite outcomes. Biometrics, 10.1111/biom.13382

## See also

[`pwreg`](https://lmaowisc.github.io/WR/reference/pwreg.md),
[`print.pwreg`](https://lmaowisc.github.io/WR/reference/print.pwreg.md)

## Examples

``` r
library(WR)
head(non_ischemic)
#>   ID time status trt_ab age sex Black.vs.White Other.vs.White   bmi bipllvef
#> 1  1  221      2      0  62   1              0              0 25.18    32.24
#> 2  1  383      0      0  62   1              0              0 25.18    32.24
#> 3  2   23      2      0  75   1              1              0 22.96    21.71
#> 4  2 1400      0      0  75   1              1              0 22.96    21.71
#> 5  5    7      2      0  48   1              1              0 34.37    22.97
#> 6  5   10      1      0  48   1              1              0 34.37    22.97
#>   hyperten COPD diabetes acei betab smokecurr
#> 1        0    0        0    0     1         1
#> 2        0    0        0    0     1         1
#> 3        1    0        0    0     1         0
#> 4        1    0        0    0     1         0
#> 5        1    0        0    0     1         0
#> 6        1    0        0    0     1         0

# Randomly sample 200 subjects from non_ischemic data
id_unique <-unique(non_ischemic$ID)
set.seed(2019)
id_sample <- sample(id_unique, 200)
non_ischemic_reduce <- non_ischemic[non_ischemic$ID %in% id_sample, ]

# Use the reduced non_ischemic data for analysis
nr <- nrow(non_ischemic_reduce)
p <- ncol(non_ischemic_reduce)-3
ID <- non_ischemic_reduce[,"ID"]
time <- non_ischemic_reduce[,"time"]
status <- non_ischemic_reduce[,"status"]
Z <- as.matrix(non_ischemic_reduce[,4:(3+p)],nr,p)
pwreg.obj <- pwreg(time=time,status=status,Z=Z,ID=ID)
score.obj <- score.proc(pwreg.obj)
#plot the standardized score process for the first covariate
plot(score.obj, k = 1)
```
