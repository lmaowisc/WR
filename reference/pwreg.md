# Fit a standard proportional win-fractions (PW) regression model

Fit a standard proportional win-fractions (PW) regression model.

## Usage

``` r
pwreg(
  ID,
  time,
  status,
  Z,
  rho = 0,
  strata = NULL,
  fixedL = TRUE,
  eps = 1e-04,
  maxiter = 50
)
```

## Arguments

- ID:

  a vector of unique subject-level identifiers.

- time:

  a vector of event times.

- status:

  a vector of event type labels. 0: censoring, 1:death and 2: non-fatal
  event.

- Z:

  a matrix or a vector of covariates.

- rho:

  a non-negative number as the power of the survival function used in
  the weight. Default (`rho=0`) is recommended. If there is a \`strata\`
  argument, then \`rho\` is ignored.

- strata:

  a vector of stratifying variable if a stratified model is desired.

- fixedL:

  logical variable indicating which variance estimator to be used. If
  \`TRUE\`, the type I variance estimator (for a small number strata) is
  used; otherwise the type II variance estimator (for a large number
  strata) is used.

- eps:

  precision for the convergence of Newton-Raphson algorithm.

- maxiter:

  maximum number of iterations allow for the Newton-Raphson algorithm.

## Value

An object of class `pwreg` with the following components. `beta`:a
vector of estimated regression coefficients. `Var`:estimated covariance
matrix for `beta`. `conv:` boolean variable indicating whether the
algorithm converged within the maximum number of iterations.

## References

Mao, L. and Wang, T. (2020). A class of proportional win-fractions
regression models for composite outcomes. Biometrics, 10.1111/biom.13382

Wang, T. and Mao, L. (2021+). Stratified Proportional Win-fractions
Regression Analysis.

## See also

[`score.proc`](https://lmaowisc.github.io/WR/reference/score.proc.md),
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
id_unique <-unique(non_ischemic$ID)

# Randomly sample 200 subjects from non_ischemic data
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
## unstratified analysis
pwreg.obj <- pwreg(time=time,status=status,Z=Z,ID=ID)
print(pwreg.obj)
#> Call:
#> pwreg(ID = ID, time = time, status = status, Z = Z)
#> 
#> Proportional win-fractions regression models for priority-adjusted composite endpoint
#> 
#>     (Mao and Wang, 2020):
#> 
#> Total number of pairs: 19900 
#> Wins-losses on death:  1306 (6.6%) 
#> Wins-losses on non-fatal event:  15345 (77.1%) 
#> Indeterminate pairs 3249 (16.3%) 
#> 
#> Newton-Raphson algorithm converged in 5 iterations.
#> 
#> Overall test: chisq test with 13 degrees of freedom; 
#>  Wald statistic 25.1 with p-value 0.02230495 
#> 
#> Estimates for Regression parameters:
#>                  Estimate         se z.value p.value  
#> trt_ab          0.1897070  0.2017570  0.9403 0.34708  
#> age            -0.0107734  0.0087005 -1.2383 0.21562  
#> sex            -0.4367173  0.2053482 -2.1267 0.03344 *
#> Black.vs.White -0.1674782  0.2262604 -0.7402 0.45918  
#> Other.vs.White -0.3121187  0.7777347 -0.4013 0.68819  
#> bmi            -0.0352119  0.0158922 -2.2157 0.02671 *
#> bipllvef        0.0273218  0.0123319  2.2155 0.02672 *
#> hyperten       -0.3656759  0.2391985 -1.5288 0.12633  
#> COPD           -0.3244369  0.3315767 -0.9785 0.32784  
#> diabetes        0.0348824  0.2128076  0.1639 0.86980  
#> acei           -0.2307814  0.2835244 -0.8140 0.41566  
#> betab          -0.2781857  0.4866345 -0.5717 0.56756  
#> smokecurr       0.0002984  0.2631628  0.0011 0.99910  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> 
#> Point and interval estimates for the win ratios:
#>                Win Ratio 95% lower CL 95% higher CL
#> trt_ab         1.2088953    0.8140532     1.7952486
#> age            0.9892844    0.9725575     1.0062990
#> sex            0.6461541    0.4320593     0.9663374
#> Black.vs.White 0.8457951    0.5428401     1.3178269
#> Other.vs.White 0.7318946    0.1593821     3.3609149
#> bmi            0.9654009    0.9357939     0.9959445
#> bipllvef       1.0276985    1.0031568     1.0528406
#> hyperten       0.6937276    0.4340930     1.1086516
#> COPD           0.7229343    0.3774507     1.3846417
#> diabetes       1.0354980    0.6823498     1.5714169
#> acei           0.7939130    0.4554456     1.3839147
#> betab          0.7571562    0.2917168     1.9652127
#> smokecurr      1.0002984    0.5972072     1.6754604
#> 
if (FALSE) { # \dontrun{
## stratified PW by sex
sex<-Z[,3]
## take out sex from the covariate matrix
Z1<-Z[,-3]
pwreg.obj1 <- pwreg(time=time,status=status,Z=Z1,ID=ID,strata=sex)
print(pwreg.obj1)
} # }
```
