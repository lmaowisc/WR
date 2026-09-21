# Generalized win ratio tests

Perform stratified two-sample test of possibly recurrent nonfatal event
and death using the recommended last-event assisted win ratio (LWR),
and/or naive win ratio (NWR) and first-event assisted win ratio (FWR)
(Mao et al., 2022). The LWR and FWR reduce to the standard win ratio of
Pocock et al. (2012).

## Usage

``` r
WRrec(ID, time, status, trt, strata = NULL, naive = FALSE)
```

## Arguments

- ID:

  A vector of unique patient identifiers.

- time:

  A numeric vector of event times.

- status:

  A vector of event type variable; 2 = recurrent event, 1 = death, and 0
  = censoring.

- trt:

  A vector of binary treatment indicators.

- strata:

  A vector of categorical variable for strata; Default is NULL, which
  leads to unstratified analysis.

- naive:

  If TRUE, results for NWR and FWR will be provided in addition to LWR;
  Default is FALSE, which gives LWR only.

## Value

An object of class `WRrec`, which contains the following elements.

- theta:

  A bivariate vector of win/loss fractions by LWR.

- log.WR, se:

  Log-win ratio estimate and its standard error by LWR.

- pval:

  \\p\\-value by the LWR test.

- theta.naive:

  A bivariate vector of win/loss fractions by NWR.

- log.WR.naive, se.naive:

  Log-win ratio estimate and its standard error by NWR.

- theta.FI:

  A bivariate vector of win/loss fractions by FWR.

- log.WR.FI, se.FI:

  Log-win ratio estimate and its standard error by FWR.

- ...:

## References

Mao, L., Kim, K. and Li, Y. (2022). On recurrent-event win ratio.
Statistical Methods in Medical Research, under review.

Pocock, S., Ariti, C., Collier, T., and Wang, D. (2012). The win ratio:
a new approach to the analysis of composite endpoints in clinical trials
based on clinical priorities. European Heart Journal, 33, 176–182.

## See also

[`print.WRrec`](https://lmaowisc.github.io/WR/reference/print.WRrec.md).

## Examples

``` r
## load the HF-ACTION trial data
library(WR)
head(hfaction_cpx9)
#>        patid       time status trt_ab age60
#> 1 HFACT00001  7.2459016      2      0     1
#> 2 HFACT00001 12.5573770      0      0     1
#> 3 HFACT00002  0.7540984      2      0     1
#> 4 HFACT00002  4.2950820      2      0     1
#> 5 HFACT00002  4.7540984      2      0     1
#> 6 HFACT00002 45.9016393      0      0     1
dat<-hfaction_cpx9
## Comparing exercise training to usual care by LWR, FWR, and NWR
obj<-WRrec(ID=dat$patid,time=dat$time,status=dat$status,
          trt=dat$trt_ab,strata=dat$age60,naive=TRUE)
## print the results
obj
#> Call:
#> WRrec(ID = dat$patid, time = dat$time, status = dat$status, trt = dat$trt_ab, 
#>     strata = dat$age60, naive = TRUE)
#> 
#>             N Rec. Event Death Med. Follow-up
#> Control   221        571    57       28.62295
#> Treatment 205        451    36       27.57377
#> 
#> Analysis of last-event-assisted WR (LWR; recommended), first-event-assisted WR (FWR), and naive WR (NWR):
#>     Win prob Loss prob WR (95% CI)*      p-value
#> LWR 50.4%    38.2%     1.32 (1.05, 1.66) 0.0189 
#> FWR 50.4%    38.3%     1.32 (1.04, 1.66) 0.0202 
#> NWR 47%      35%       1.34 (1.05, 1.72) 0.0193 
#> -----
#> *Note: The scale of WR should be interpreted with caution as it depends on 
#> censoring distribution without modeling assumptions.
```
