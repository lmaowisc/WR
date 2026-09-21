# A subset of the HF-ACTION study data on non-ischemic heart failure patients with full covariate measurement.

These are a subset of the data on 451 non-ischemic patients in the
HF-ACTION study will complete baseline covariates.

## Usage

``` r
non_ischemic
```

## Format

A data frame with 751 rows and 16 variables:

- ID:

  subject IDs

- time:

  event times (days)

- status:

  event status; 0:censoring, 1:death, 2:hospitalization

- trt_ab:

  treatment indicator: 1=exercise training; 0=usual care

- age:

  patient age in years

- sex:

  1=female; 2=male

- Black.vs.White:

  1=black; 0=otherwise

- Other.vs.White:

  1=race other than black or white; 0=otherwise

- bmi:

  body mass index

- bipllvef:

  (biplane) left-ventricular ejection fraction

- hyperten:

  indicator for history of hypertension

- COPD:

  indicator for history of COPD

- diabetes:

  indicator for history of diabetes

- acei:

  indicator for current use of ACE inhibitors

- betab:

  indicator for current use of beta blockers

- smokecurr:

  indicator for current smoker

## References

O'Connor, C. M., Whellan, D. J., Lee, K. L., Keteyian, S. J., Cooper, L.
S., Ellis, S. J., Leifer, E. S., Kraus, W. E., Kitzman, D. W.,
Blumenthal, J. A. et al. (2009). Efficacy and safety of exercise
training in patients with chronic heart failure: HF-ACTION randomized
controlled trial. Journal of the American Medical Association, 301,
1439–1450.
