# A subset of the German Breast Cancer study data

These are a subset of the German Breast Cancer study data.

## Usage

``` r
gbc
```

## Format

A data frame with 985 rows and 12 variables:

- id:

  subject IDs

- time:

  event times (months)

- status:

  event status; 0:censoring, 1:death, 2:cancer recurrence

- hormone:

  treatment indicator: 1=Hormone therapy; 2=standard therapy

- age:

  age at diagnosis (years)

- menopause:

  menopausal Status; 1=No; 2=Yes

- size:

  tumor size

- grade:

  tumor grade, 1-3

- nodes:

  number of nodes involved

- prog_recp:

  number of progesterone receptors

- estrg_recp:

  number of estrogen receptors

## References

Sauerbrei, W., Royston, P., Bojar, H., Schmoor, C. and Schumacher, M.
(1999). Modelling the effects of standard prognostic factors in
node-positive breast cancer. German Breast Cancer Study Group (GBSG).
British Journal of Cancer, 79, 1752–1760.

Hosmer, D.W. and Lemeshow, S. and May, S. (2008) Applied Survival
Analysis: Regression Modeling of Time to Event Data: Second Edition,
John Wiley and Sons Inc., New York, NY
