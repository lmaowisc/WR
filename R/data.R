#' A subset of the German Breast Cancer study data
#'
#' @description These are a subset of the German Breast Cancer study data.
#' @format A data frame with 985 rows and 12 variables:
#' \describe{
#'   \item{id}{subject IDs}
#'   \item{time}{event times (months)}
#'   \item{status}{event status; 0:censoring, 1:death, 2:cancer recurrence}
#'   \item{hormone}{treatment indicator: 1=Hormone therapy; 2=standard therapy}
#'   \item{age}{age at diagnosis (years)}
#'   \item{menopause}{menopausal Status; 1=No; 2=Yes}
#'   \item{size}{tumor size}
#'   \item{grade}{tumor grade, 1-3}
#'   \item{nodes}{number of nodes involved}
#'   \item{prog_recp}{number of progesterone receptors}
#'   \item{estrg_recp}{number of estrogen receptors}
#'   }
#' @references Sauerbrei, W., Royston, P., Bojar, H., Schmoor, C. and Schumacher, M.
#' (1999). Modelling the effects of standard prognostic factors in node-positive
#' breast cancer. German Breast Cancer Study Group (GBSG). British Journal of Cancer,
#' 79, 1752–1760.
#' @references Hosmer, D.W. and Lemeshow, S. and May, S. (2008) Applied Survival
#' Analysis: Regression Modeling of Time to Event Data: Second Edition, John Wiley
#' and Sons Inc., New York, NY
"gbc"


#' A subset of the HF-ACTION study data on non-ischemic heart failure patients with full
#' covariate measurement.
#'
#' @description These are a subset of the data on 451 non-ischemic patients in the HF-ACTION study
#' will complete baseline covariates.
#' @format A data frame with 751 rows and 16 variables:
#' \describe{
#'   \item{ID}{subject IDs}
#'   \item{time}{event times (days)}
#'   \item{status}{event status; 0:censoring, 1:death, 2:hospitalization}
#'   \item{trt_ab}{treatment indicator: 1=exercise training; 0=usual care}
#'   \item{age}{patient age in years}
#'   \item{sex}{1=female; 2=male}
#'   \item{Black.vs.White}{1=black; 0=otherwise}
#'   \item{Other.vs.White}{1=race other than black or white; 0=otherwise}
#'   \item{bmi}{body mass index}
#'   \item{bipllvef}{(biplane) left-ventricular ejection fraction}
#'   \item{hyperten}{indicator for history of hypertension}
#'   \item{COPD}{indicator for history of COPD}
#'   \item{diabetes}{indicator for history of diabetes}
#'   \item{acei}{indicator for current use of ACE inhibitors}
#'   \item{betab}{indicator for current use of beta blockers}
#'   \item{smokecurr}{indicator for current smoker}
#'   }
#' @references O'Connor, C. M., Whellan, D. J., Lee, K. L., Keteyian, S. J., Cooper, L. S., Ellis, S. J.,
#' Leifer, E. S., Kraus, W. E., Kitzman, D. W., Blumenthal, J. A. et al. (2009). Efficacy and
#' safety of exercise training in patients with chronic heart failure: HF-ACTION randomized
#' controlled trial. Journal of the American Medical Association, 301, 1439--1450.
"non_ischemic"


#' A subset of the HF-ACTION study data on high-risk non-ischemic heart failure patients
#'
#' @description These are data on a subgroup of 426 high-risk non-ischemic patients in the HF-ACTION study.
#' @format A data frame with 1,448 rows and 5 variables:
#' \describe{
#'   \item{patid}{patient ID}
#'   \item{time}{event times (months)}
#'   \item{status}{event status; 0:censoring, 1:death, 2:hospitalization}
#'   \item{trt_ab}{treatment indicator: 1: exercise training, 0: usual care}
#'   \item{age60}{1: 60 years or older, 0: less than 60 years old}
#'   }
#' @references O'Connor, C. M., Whellan, D. J., Lee, K. L., Keteyian, S. J., Cooper, L. S., Ellis, S. J.,
#' Leifer, E. S., Kraus, W. E., Kitzman, D. W., Blumenthal, J. A. et al. (2009). Efficacy and
#' safety of exercise training in patients with chronic heart failure: HF-ACTION randomized
#' controlled trial. Journal of the American Medical Association, 301, 1439--1450.
"hfaction_cpx9"
