


########################################################################
#             Functions for composite outcomes (under Gumbel model)   ##
########################################################################

# The delta function using numerical integration
# the output is a vector of (delta_1,delta_2)

delta.fun = function(baseline){
  lambda_D=baseline$lambda_D
  lambda_H=baseline$lambda_H
  kappa=baseline$kappa
  tau_b=baseline$tau_b
  tau=baseline$tau
  lambda_L=baseline$lambda_L

  # the input "x" is a vector containing two entries (s,t)
  SD = function(x){
    exp(-lambda_D*x[2])
  }

  g_tilde = function(x){
    if(x[1] > tau - tau_b){
      2/tau_b^2*exp(-2*lambda_L*x[1])*(tau-x[1])*(lambda_L*(tau-x[1])+1)
    }else{
      2*lambda_L*exp(-2*lambda_L*x[1])
    }
  }
  delta1 = adaptIntegrate(function(x){
    tem = ((lambda_D*x[1])^kappa+
             (lambda_H*x[2])^kappa)^(1/kappa)
    H0 = exp(-tem)
    lambda1= tem^(1-2*kappa)*(1-kappa)*(lambda_D*lambda_H*x[1])^kappa*x[2]^(kappa-1)

    return(ifelse(x[2]<=x[1],1,0)*(lambda_D*SD(x)^2 + H0^2*lambda1)*g_tilde(x))
  },
  c(0,0), c(tau,tau))

  delta2 = adaptIntegrate(function(x){
    tem = ((lambda_D*x[1])^kappa+
             (lambda_H*x[2])^kappa)^(1/kappa)
    H0 = exp(-tem)
    lambda2 =  (1-kappa)*tem^(1-2*kappa) *
      (lambda_H)^(2*kappa)*x[2]^(2*kappa-1)+
      kappa * tem^(1-kappa) *
      x[2]^(kappa-1)*lambda_H^kappa

    return(ifelse(x[2]<=x[1],1,0) * (H0^2*lambda2)*g_tilde(x))
  },
  c(0,0), c(tau,tau))

  return(c(delta1 = delta1$integral,delta2 = delta2$integral))
}





###############################################################
# Win function for composite outcomes
# see argument "winfun" in WR.anal()
#############################################################
Wp = function(y1,y2){
  ifelse(y2[1]<min(c(y1[1],y1[3],y2[3])),1,0) +
    ifelse((min(c(y1[1],y2[1])) > min(c(y1[3], y2[3])) ) &
             (y2[2] < min(c(y1[2],y1[3],y2[3]))),1,0)
}

## This is the win-loss difference function
win.comp = function(y1,y2){Wp(y1,y2) - Wp(y2,y1)}



# monte-carlo integration of zeta_0^2

zeta2.fun=function(baseline,N=2000,seed){


  lambda_D=baseline$lambda_D
  lambda_H=baseline$lambda_H
  kappa=baseline$kappa
  tau_b=baseline$tau_b
  tau=baseline$tau
  lambda_L=baseline$lambda_L


  #simulate death, hosp, and censoring
  # Death
  set.seed(seed)
  outcome=rgumbel(N,alpha=kappa,dim=2)
  D=-log(outcome[,1])/lambda_D
  H=-log(outcome[,2])/lambda_H
  set.seed(seed)
  Ca=tau_b*runif(N)+tau-tau_b
  set.seed(seed)
  L=rexp(N)/lambda_L
  C=pmin(Ca,L)



  dat=cbind(D,H,C)

  tmp=combn(1:nrow(dat),2,function(x) win.comp(dat[x[1],],dat[x[2],]))



  WL.matrix=matrix(0,N,N)

  WL.matrix[lower.tri(WL.matrix)]=tmp

  WL.matrix=t(WL.matrix)-WL.matrix

  w0=sum(abs(WL.matrix))/(2*N*(N-1))
  ranks=rowSums(WL.matrix)/(N-1)
  zeta2=sum(ranks^2)/(N-3)

  return(list(w0=w0,zeta2=zeta2))


}

#######################################################
# Function to compute the baseline parameters
# in the sample size formula of Pocock's win ratio
######################################################




#' Compute the baseline parameters needed for sample size calculation for
#' standard win ratio test
#'
#' @description Compute the baseline parameters \eqn{\zeta_0^2} and \eqn{\boldsymbol\delta_0}
#' needed for sample size calculation for standard win ratio test (see \code{\link{WRSS}}).
#' The calculation is based
#' on a Gumbel--Hougaard copula model for survival time \eqn{D^{(a)}} and nonfatal event
#' time \eqn{T^{(a)}} for group \eqn{a} (1: treatment; 0: control):
#'\deqn{{P}(D^{(a)}>s, T^{(a)}>t) =\exp\left(-\left[\left\{\exp(a\xi_1)\lambda_Ds\right\}^\kappa+
#' \left\{\exp(a\xi_2)\lambda_Ht\right\}^\kappa\right]^{1/\kappa}\right),}
#' where \eqn{\xi_1} and \eqn{\xi_2} are the component-wise log-hazard ratios to be used
#' as effect size in \code{\link{WRSS}}.
#' We also assume that patients are recruited uniformly over the period \eqn{[0, \tau_b]}
#' and followed until time \eqn{\tau}  (\eqn{\tau\geq\tau_b}), with an exponential
#' loss-to-follow-up hazard \eqn{\lambda_L}.
#' @param lambda_D Baseline hazard \eqn{\lambda_D} for death.
#' @param lambda_H Baseline hazard \eqn{\lambda_H} for nonfatal event.
#' @param kappa Gumbel--Hougaard copula correlation parameter \eqn{\kappa}.
#' @param tau_b Length of the initial (uniform) accrual period \eqn{\tau_b}.
#' @param tau Total length of follow-up \eqn{\tau}.
#' @param lambda_L Exponential hazard rate \eqn{\lambda_L} for random loss to follow-up.
#' @param N Simulated sample size for monte-carlo integration.
#' @param seed Seed for monte-carlo simulation.
#' @return A list containing real number \code{zeta2} for \eqn{\zeta_0^2}
#' and bivariate vector \code{delta} for \eqn{\boldsymbol\delta_0}.
#' @seealso \code{\link{gumbel.est}}, \code{\link{WRSS}}
#' @import gumbel
#' @import cubature
#' @importFrom utils combn tail
#' @importFrom stats rexp runif
#' @keywords WRSS
#' @examples
#' # see the example for WRSS
#' @export
#' @references Mao, L., Kim, K. and Miao, X. (2021). Sample size formula for general win ratio analysis.
#' Biometrics, https://doi.org/10.1111/biom.13501.
base=function(lambda_D,lambda_H,kappa,tau_b,tau,lambda_L,N=1000,seed=12345){
  baseline=list(lambda_D=lambda_D,lambda_H=lambda_H,kappa=kappa,tau_b=tau_b,
                tau=tau,lambda_L=lambda_L)

  obj=zeta2.fun(baseline,N,seed)
  zeta2=obj$zeta2
  w0=obj$w0

  delta=delta.fun(baseline)

  return(list(zeta2=zeta2,w0=w0,delta=delta))

}



#' Compute the sample size for standard win ratio test
#'
#' @description Compute the sample size for standard win ratio test.
#'
#' @param xi A bivariate vector of hypothesized component-wise (treatment-to-control) log-hazard
#' ratios under the Gumbel--Hougaard copula model described in \code{\link{base}}.
#' @param bparam A list containing baseline parameters \code{zeta2} for \eqn{\zeta_0^2}
#' and \code{delta} for \eqn{\boldsymbol\delta_0}; Can directly use the output of \link{base}.
#' @param q Proportion of patients assigned to treatment.
#' @param alpha Type I error rate.
#' @param side 2-sided or 1-sided test.
#' @param power Target power.
#' @return A list containing \code{n}, the computed sample size.
#' @seealso \code{\link{gumbel.est}}, \code{\link{base}}
#' @importFrom stats qnorm
#' @keywords WRSS
#' @export
#' @references Mao, L., Kim, K. and Miao, X. (2021). Sample size formula for general win ratio analysis.
#' Biometrics, https://doi.org/10.1111/biom.13501.
#' @export
#' @examples
#' # The following is not run in package checking to save time.
#' \dontrun{
#'## load the package and pilot dataset
#'library(WR)
#'head(hfaction_cpx9)
#'dat<-hfaction_cpx9
#'## subset to control group
#'pilot<-dat[dat$trt_ab==0,]
#'
#'## get the data ready for gumbel.est()
#'id<-pilot$patid
#'## convert time from month to year
#'time<-pilot$time/12
#'status<-pilot$status
#'## compute the baseline parameters for the Gumbel--Hougaard
#'## copula for death and hospitalization
#'gum<-gumbel.est(id, time, status)
#'
#'## get the baseline parameters
#'lambda_D<-gum$lambda_D
#'lambda_H<-gum$lambda_H
#'kappa<-gum$kappa
#'## set up design parameters and use base()
#'## to calculate bparam for WRSS()
#'# max follow-up 4 years
#'tau<-4
#'# 3 years of initial accrual
#'tau_b<-3
#'# loss to follow-up rate
#'lambda_L=0.05
#'# compute the baseline parameters
#'bparam<-base(lambda_D,lambda_H,kappa,tau_b,tau,lambda_L)
#'bparam
#'
#'## sample size with power=0.8 under hazard ratios
#'## 0.9 and 0.8 for death and hospitalization, respectively.
#'WRSS(xi=log(c(0.9,0.8)),bparam=bparam,q=0.5,alpha=0.05,
#'     power=0.8)$n
#'## sample size under the same set-up but with power 0.9
#'WRSS(xi=log(c(0.9,0.8)),bparam=bparam,q=0.5,alpha=0.05,
#'     power=0.9)$n
#' }
WRSS=function(xi,bparam,q=0.5,alpha=0.05,side=2,power=0.8){

  xi=as.matrix(xi)
  #sample size formula
  za=qnorm(1-alpha/side)
  zb=qnorm(power)

  zeta2=bparam$zeta2
  delta=bparam$delta

  n=zeta2*(za+zb)^2/(q*(1-q)*(t(xi)%*%delta)^2)

  return(list(zeta2=zeta2,delta=delta,n=n,xi=xi))

}



#===========================================================================================
# Estimating the baseline parameters lambda_D, lambda_H
# and kappa from pilot study data
#===========================================================================================

#' Estimate baseline parameters in the Gumbel--Hougaard model for sample size
#' calculation using pilot data
#'
#' @description Estimate baseline parameters in the Gumbel--Hougaard model
#' described in \code{\link{base}} for sample size calculation using pilot study data.
#'
#' @param id A vector of unique patient identifiers.
#' @param time A numeric vector of event times.
#' @param status A vector of event type variable; 2 = nonfatal event, 1 = death,
#' and 0 = censoring.
#' @return A list containing \code{lambda_D} for \eqn{\lambda_D},
#' \code{lambda_H} for \eqn{\lambda_H}, and \code{kappa} for \eqn{\kappa}
#' in the Gumbel--Hougaard model.
#' @seealso \code{\link{base}}, \code{\link{WRSS}}
#' @examples
#' # see the example for WRSS
#' @export
#' @keywords WRSS
#' @references Mao, L., Kim, K. and Miao, X. (2021). Sample size formula for general win ratio analysis.
#' Biometrics, https://doi.org/10.1111/biom.13501.
gumbel.est=function(id, time, status){

  data=data.frame(id,time,status)
  #Sort the data by time within each id
  o=order(data$id,data$time)
  data=data[o,]

  ######################################
  #Composite (time-to-first) event rate#
  ######################################
  #subset to first event data
  #Sort the data by time within each id

  # get the first row for each id
  data.CE=data[!duplicated(data$id),]
  # event rate
  lambda_CE=sum(data.CE$status>0)/sum(data.CE$time)
  #########################################
  # Cause-specific hazard rate for hosp
  ##########################################
  lambda_sH=sum(data.CE$status==2)/sum(data.CE$time)


  ######################
  #     Death rate    ##
  ######################
  data.D=data[data$status<2,]
  # event rate
  lambda_D=sum(data.D$status)/sum(data.D$time)




  kappa=log(1-lambda_sH/lambda_CE)/log(lambda_D/lambda_CE)

  lambda_H=(lambda_CE^kappa-lambda_D^kappa)^(1/kappa)

  return(list(lambda_D=lambda_D,lambda_H=lambda_H,kappa=kappa))

}


