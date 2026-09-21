# Proportional win-fractions (PW) regression of composite endpoints of death and nonfatal event

This vignette demonstrates the use of the `WR` package in fitting the
proportional win-fractions (PW) regression model for prioritized
composite endpoints consisting of death and a nonfatal event (Mao and
Wang, 2020, *Biometrics*). The PW model can be viewed as the regression
version of the two-sample win ratio proposed by Pocock et al. (2012).

## MODEL SPECIFICATION

Let $`D`$ denote the survival time, $`T`$ time to the first nonfatal
event like hospitalization, and $`\boldsymbol Z`$ a $`p`$-vector of
covariates. The composite outcome is $`\boldsymbol Y=(D, T)`$ with $`D`$
prioritized over $`T`$. If the $`i`$th and $`j`$th patients are both
followed up to time $`t`$, define
``` math
\mathcal W(\boldsymbol Y_i, \boldsymbol Y_j)(t) = I(\mbox{subject $i$ wins against subject $j$ by $t$}) = I(D_j < D_i \wedge t) + I(D_i \wedge D_j > t, T_j < T_i \wedge t),
```
where $`a\wedge b=\min(a,b)`$. This ‘’win indicator’’ uses the
sequential comparison rule of Pocock et al. (2012).

Then, the covariate-specific win ratio at $`t`$ is
``` math
WR(t;\boldsymbol Z_i, \boldsymbol Z_j):=\frac{E\{\mathcal W(\boldsymbol Y_i, \boldsymbol Y_j)(t)\mid \boldsymbol Z_i, \boldsymbol Z_j\}}{E\{\mathcal W(\boldsymbol Y_j, \boldsymbol Y_i)(t)\mid \boldsymbol Z_i, \boldsymbol Z_j\}}.
```
For example, if subject $`i`$ is from the treatment arm with
$`\boldsymbol Z_i=1`$ and subject $`j`$ is from the treatment arm with
$`\boldsymbol Z_j=0`$, then $`WR(t;\boldsymbol Z_i, \boldsymbol Z_j)`$
is precisely the estimand of Pocock’s win ratio comparing the treatment
to the control when all subjects in both arms are followed to time
$`t`$.

The PW model specifies that \$\$\mbox{PW Model: }\hspace{3mm}
WR(t;\boldsymbol Z_i, \boldsymbol Z_j)=\exp\\\boldsymbol\beta^{\rm
T}(\boldsymbol Z_i -\boldsymbol Z_j)\\.\$\$ Clearly, the PW model
assumes that win ratio is invariant to the follow-up time
(proportionality assumption). Under the model, the components of
$`\boldsymbol\beta`$ can be interpreted as the log win ratios associated
with unit increases in the corresponding covariates. Under the PW model,
we can obtain consistent estimates for the parameter
$`\boldsymbol\beta`$ based on censored data regardless of distribution
of the censoring time $`C`$ as long as
``` math
(C\perp \boldsymbol Y)\mid \boldsymbol Z.
```

## BASIC SYNTAX

The input data must be in the “long format”, with an `ID` vector
containing unique patient-level identifiers. In addition, we need a
`time` vector containing the event times and a `status` vector
indicating the corresponding cause of the event. The vector `status`
should be coded as `1`=death; `2`=non-fatal event; `0`=censoring. In the
case of recurrent non-fatal event, multiple rows with `status`=2 are
allowed. However, by nature of the method, only time to the first
episode will be used. Finally, we need a covariate matrix `Z` with the
same row as `ID`. Each column of `Z` represents a covariate. All
covariates need to be time-constant.

The main function to fit the PW model is `pwreg(time,status, Z, ID)`.
The function returns an object of class `pwreg` with a `beta` vector for
$`\widehat{\boldsymbol\beta}`$ and a `Var` matrix for
$`\text{var}(\widehat{\boldsymbol\beta})`$. For details, refer to
documentation of the `WR` package.

## AN EXAMPLE WITH THE HF-ACTION TRIAL

We consider a dataset from the HF-ACTION study consisting of 451
non-ischemic heart failure patients. The study was conducted between
April 2003 through Feb 2007 at 82 sites in the USA, Canada, and France
(O’Connor et al., 2009). The study objective was to assess the effect of
adding aerobic exercise training to usual care on the patient’s CV
outcomes. The primary endpoint was a composite of all-cause death and
all-cause hospitalization.

We first load the `WR` package and the analysis dataset `non_ischemic`.

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
```

We re-label the covariates with informative names.

``` r

colnames(non_ischemic)[4:16]=c(
  "Training vs Usual","Age (year)","Male vs Female","Black vs White", 
  "Other vs White", "BMI","LVEF","Hypertension","COPD","Diabetes",
  "ACE Inhibitor","Beta Blocker", "Smoker"
)
```

Compute the sample size the median length of follow-up.

``` r

# sample size
length(unique(non_ischemic$ID))
#> [1] 451
# median length of follow-up time
median(non_ischemic$time[non_ischemic$status<2])/30.5
#> [1] 31.63934
```

So we indeed have $`n=451`$ unique patients with a median follow-up of
31.6 months.

Next, we use the
[`pwreg()`](https://lmaowisc.github.io/WR/reference/pwreg.md) function
to fit the PW model:

``` r

# get the number of rows and number of covariates.
nr <- nrow(non_ischemic)
p <- ncol(non_ischemic)-3

# extract ID, time, status and covariates matrix Z from the data.
# note that: ID, time and status should be column vector.
# covariatesZ should be (nr, p) matrix.
ID <- non_ischemic[,"ID"]
time <- non_ischemic[,"time"]
status <- non_ischemic[,"status"]
Z <- as.matrix(non_ischemic[,4:(3+p)],nr,p)


# pass the parameters into the function
pwreg.obj <- pwreg(time=time,status=status,Z=Z,ID=ID)
print(pwreg.obj)
#> Call:
#> pwreg(ID = ID, time = time, status = status, Z = Z)
#> 
#> Proportional win-fractions regression models for priority-adjusted composite endpoint
#> 
#>     (Mao and Wang, 2020):
#> 
#> Total number of pairs: 101475 
#> Wins-losses on death:  7644 (7.5%) 
#> Wins-losses on non-fatal event:  78387 (77.2%) 
#> Indeterminate pairs 15444 (15.2%) 
#> 
#> Newton-Raphson algorithm converged in 5 iterations.
#> 
#> Overall test: chisq test with 13 degrees of freedom; 
#>  Wald statistic 24.9 with p-value 0.02392931 
#> 
#> Estimates for Regression parameters:
#>                     Estimate         se z.value p.value  
#> Training vs Usual  0.1906687  0.1264658  1.5077 0.13164  
#> Age (year)        -0.0128306  0.0057285 -2.2398 0.02510 *
#> Male vs Female    -0.1552923  0.1294198 -1.1999 0.23017  
#> Black vs White    -0.3026335  0.1461330 -2.0709 0.03836 *
#> Other vs White    -0.3565390  0.3424360 -1.0412 0.29779  
#> BMI               -0.0181310  0.0097582 -1.8580 0.06316 .
#> LVEF               0.0214905  0.0086449  2.4859 0.01292 *
#> Hypertension      -0.0318291  0.1456217 -0.2186 0.82698  
#> COPD              -0.4023069  0.2066821 -1.9465 0.05159 .
#> Diabetes           0.0703990  0.1419998  0.4958 0.62006  
#> ACE Inhibitor     -0.1068201  0.1571317 -0.6798 0.49662  
#> Beta Blocker      -0.5344979  0.3289319 -1.6250 0.10417  
#> Smoker            -0.0602350  0.1682826 -0.3579 0.72039  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> 
#> Point and interval estimates for the win ratios:
#>                   Win Ratio 95% lower CL 95% higher CL
#> Training vs Usual 1.2100585    0.9444056     1.5504374
#> Age (year)        0.9872513    0.9762288     0.9983983
#> Male vs Female    0.8561648    0.6643471     1.1033663
#> Black vs White    0.7388699    0.5548548     0.9839127
#> Other vs White    0.7000951    0.3578286     1.3697431
#> BMI               0.9820323    0.9634287     1.0009952
#> LVEF              1.0217231    1.0045572     1.0391823
#> Hypertension      0.9686721    0.7281543     1.2886357
#> COPD              0.6687755    0.4460178     1.0027865
#> Diabetes          1.0729362    0.8122757     1.4172433
#> ACE Inhibitor     0.8986873    0.6604773     1.2228110
#> Beta Blocker      0.5859634    0.3075270     1.1164977
#> Smoker            0.9415433    0.6770144     1.3094312
```

The output consists of three parts. The first part presents some
descriptive statistics on the proportions of win-loss status among all
$`{n\choose 2}=101,475`$ pairs. According to the output, $`7.5\%`$ of
them are determined by death; $`77.2\%`$ by hospitalization, and the
remaining $`7.2\%`$ are indeterminate. It also reports an overall (Wald)
test with $`p`$-value 0.024, suggesting that, at the conventional 0.05
level, the 13 covariates are significantly associated with the composite
outcome.

The second part presents a table for the estimates and standard errors
of the regression coefficient, along with their corresponding
$`p`$-value for testing the coefficient being zero. The third part is
perhaps the most informative, tabulating the estimated win ratios
(exponential of the regression coefficients) and their associated
$`95\%`$ confidence intervals. We can see that a patient in exercise
training is $`21\%`$ more likely to have a better priority-adjusted
composite outcome than one in usual care. However, this difference is
statistically not significant. In addition, younger age, white race,
higher LVEF are significantly associated with more favorable outcomes
than otherwise, while the beneficial effects of low BMI and absence of
COPD history are border-line significant.

To assess the effect of race on the composite outcome, we test the null
hypothesis
``` math
H_0:\beta_4=\beta_5=0.
```
We conduct a 2-df Chi-square Wald test based on
$`(\widehat\beta_4,\widehat\beta_5)^{T}`$:

``` r

#extract estimates of (\beta_4,\beta_5)
beta <- matrix(pwreg.obj$beta[4:5])
#extract estimated covariance matrix for (\beta_4,\beta_5)
Sigma <- pwreg.obj$Var[4:5,4:5]
#compute chisq statistic in quadratic form
chistats <- t(beta) %*% solve(Sigma) %*% beta  

#compare the Wald statistic with the reference
# distribution of chisq(2) to obtain the p-value
1 - pchisq(chistats, df=2)
#>           [,1]
#> [1,] 0.1016988
```

The $`p`$-value is 0.102. So the overall effect of race on the composite
outcome is non-significant.

Finally, we use the
[`score.proc()`](https://lmaowisc.github.io/WR/reference/score.proc.md)
function to plot the standardized score process for each covariate:

``` r

score.obj <- score.proc(pwreg.obj)
print(score.obj)
#> This object contains two components:
#>  't': an l-vector of times
#>  'score': a p-by-l matrix whose k'th row is the standardized score process for the k'th covariate
#>           as a function of t
#> 
#> Use 'plot(object,k=k)' to plot the k'th score process.

oldpar <- par(mfrow = par("mfrow"))
par(mfrow = c(4,4))
for(i in c(1:13)){
  plot(score.obj, k = i)
}
par(oldpar)
```

![](PW_reg_files/figure-html/unnamed-chunk-6-1.png)

Most curves are fairly patternless with suprema well bounded by 2. So we
conclude that the proportionality assumption approximately holds.

## References

- Mao, L. and Wang, T. (2020). A class of proportional win-fractions
  regression models for composite outcomes. *Biometrics*,
  <https://doi.org/10.1111/biom.13382>.

- O’Connor, C. M., Whellan, D. J., Lee, K. L., Keteyian, S. J.,
  Cooper, L. S., Ellis, S. J., Leifer, E. S., Kraus, W. E., Kitzman, D.
  W., Blumenthal, J. A. et al. (2009). Efficacy and safety of exercise
  training in patients with chronic heart failure: HF-ACTION randomized
  controlled trial. *Journal of the American Medical Association*, 301,
  1439–1450.

- Pocock, S., Ariti, C., Collier, T., and Wang, D. (2012). The win
  ratio: a new approach to the analysis of composite endpoints in
  clinical trials based on clinical priorities. *European Heart
  Journal*, 33, 176–182.
