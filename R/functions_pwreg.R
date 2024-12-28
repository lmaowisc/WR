#' @export
extract.times <- function(df){
  # Given a wide-format data, transform it to a long-format longitudinal data
  p <- ncol(df)-3
  time <- df[,"time"]
  tau <- max(time)
  status <- df[,"status"]
  Z <- matrix(df[,4:(3+p)],nrow(df),p)
  colnames(Z) <- colnames(df)[4:(3+p)]
  ID <- df[,"ID"]
  extract.by.id <- function(id){
    id_index <- (ID==id)
    status_id <- status[id_index]
    time_id <- time[id_index]
    Z_id <- matrix(Z[id_index,],sum(id_index),p)
    time_H <- ifelse( 2 %in% status_id, min(time_id[status_id==2]),tau+1)
    time_D <- ifelse( 1 %in% status_id, time_id[status_id==1],tau+1)
    time_C <- ifelse( 0 %in% status_id, time_id[status_id==0],tau+1)
    Zi <- Z_id[1,1:p]
    return(c(ID = id, timeH = time_H, timeD = time_D,timeC = time_C, Z=Zi))
  }
  df2 <- as.matrix(Reduce(rbind.data.frame, lapply(unique(ID), function(x) extract.by.id(x))))
  colnames(df2) <- c("ID","timeH","timeD","timeC",colnames(Z))
  return(list(datap=df2,tau=tau))
  #return(df2)
}

lt <- function(x,y,ep=1e-12){
  return(x < y-ep)
}

le <- function(x,y,ep=1e-12){
  return(x <= y+ep)
}

eq <- function(x,y,ep=1e-12){
  return(abs(x-y) <= ep)
}

delta.data <- function(datap,tau){
  m<-ncol(datap)
  n<-nrow(datap)

  contrast.ij<-function(x){
    di<-datap[x[1],]
    dj<-datap[x[2],]

    Hi<-di[2]
    Di<-di[3]
    Ci<-di[4]
    Zi<-di[5:m]

    Hj<-dj[2]
    Dj<-dj[3]
    Cj<-dj[4]
    Zj<-dj[5:m]

    comp<-0

    if (lt(Hj,min(Hi,Di,Dj))&&le(Hj,Ci)){
      w <- 1
      t1 <- Hj
      if (lt(Di,Dj)&&le(Di,Cj)){
        t2 <- Di
        comp <- 1
      }else{
        t2 <- tau+1
        comp <- 2
      }
    }else{
      if ((le(Dj, min(Hi,Hj,Ci))||eq(Hi,Hj))&&lt(Dj,Di)){
        w <- 1
        t1 <- Dj
        t2 <- tau+1
        comp <- 1
      }else{
        if (lt(Hi,min(Hj,Dj,Di))&&le(Hi,Cj)){
          w <- -1
          t1 <- Hi
          if (lt(Dj,Di)&le(Dj,Ci)){
            t2 <- Dj
            comp <- 1
          }else{
            t2 <- tau+1
            comp <- 2
          }
        }else{
          if ((le(Di, min(Hi,Hj,Cj))||eq(Hi,Hj))&&lt(Di,Dj)){
            w <- -1
            t1 <- Di
            t2 <- tau+1
            comp <- 1
          }else{
            w <- 0
            t1 <- 0
            t2 <- 0
          }
        }
      }
    }

    Zd <- Zi-Zj
    return(c(w,t1,t2,Zd,comp))

  }
  value <- combn(1:n, 2, FUN = contrast.ij)
  index <- combn(1:n, 2)
  outcome <- as.data.frame(t(rbind(index, value)))
  colnames(outcome) <- c("IDi","IDj","w","t1","t2",colnames(datap)[5:m],"comp")
  return(outcome)
}

int.deriv.ij <- function(pair){
  p <- (length(pair)-1)/2
  if (eq(pair[p+1],0)){
    return(matrix(0,p,p))
  }else{
    weight <- pair[1:p]
    dm <- pair[(p+2):(2*p+1)]
  }
  return(-weight%*%t(dm))
}

M.fun <- function(beta,dataf,tau){
  p <- length(beta)
  m <- nrow(dataf)
  Zd <- as.matrix(dataf[,6:(5+p)],m,p)
  eta <- exp(Zd%*%beta)
  mu <- eta/(1+eta)
  dmu <- Zd*matrix(rep(eta*(1+eta)^{-2},by=p),m,p)
  dataM.beta <- cbind(dataf,mu,dmu)
  return(dataM.beta)
}


surv.Hfun <- function(rho,data,dataf,datap,tau){
  m <- nrow(dataf)
  p <- ncol(data)-3
  Zd <- as.matrix(dataf[,6:(5+p)],m,p)

  if (rho==0){
    t <- sort(unique(datap[,"timeD"]))
    t <- t[t<=tau]
    l <- length(t)
    H <- Zd
    survZ <- NULL
  }else{
    data.surv <- data[data[,"status"]!=2,]
    time <- data.surv[,"time"]
    status <- data.surv[,"status"]
    n <- nrow(data.surv)
    Z <- as.matrix(data.surv[,4:(3+p)],n,p)
    obj.cox <- survival::coxph(survival::Surv(time,status)~Z)
    obj.cox.sum <- summary(obj.cox)
    gamma <- obj.cox$coeff
    base <- basehaz(obj.cox)
    t <- base$t
    l <- length(t)
    Lambda <- base$hazard
    survZ <- exp(-exp(Z%*%gamma)%*%t(Lambda))
    H <- Zd
  }
  return(list(H=H,t=t,survZ=survZ))
}


UEE.fun <- function(beta,H,t,rho,dataf,tau,survZ=NULL){
  p <- length(beta)
  m <- nrow(dataf)
  dataM.beta <- M.fun(beta,dataf,tau)

  t1 <- dataM.beta[,"t1"]
  t2 <- dataM.beta[,"t2"]

  if (!eq(rho,0)){
    obj1 <- outer(t,t1,"<")
    obj2 <- outer(t,t2,"<")
    ind1 <- colSums(obj1)
    ind1[ind1==0] <- 1
    ind2 <- colSums(obj2)
    ind2[ind2==0] <- 1
    dataM.beta[,c("ind1","ind2","pair_order")] <- cbind(ind1,ind2,1:m)

    EE.ij <- function(x){
      wij <- as.numeric(x["w"])
      if (eq(wij,0)){
        return(rep(0,p))
      }else{
        muij <- x["mu"]
        indij1 <- as.numeric(x["ind1"])
        indij2 <- as.numeric(x["ind2"])
        po <- as.numeric(x["pair_order"])
        IDi <- as.integer(x["IDi"])
        IDj <- as.integer(x["IDj"])
        Hij1 <- H[po,]*survZ[IDi,indij1]*survZ[IDj,indij1]
        Hij2 <- H[po,]*survZ[IDi,indij2]*survZ[IDj,indij2]
        tij2 <- as.numeric(x["t2"])
        return(Hij1*abs(wij)*(eq(wij,1)-muij)-wij*(tij2<tau)*Hij2)
        }
      }

    EE <- as.matrix(t(apply(dataM.beta,1,EE.ij)),m,p)
    if (dim(EE)[1]==1){EE <- t(EE)}

    dEE.ij <- function(x){

      wij <- as.numeric(x["w"])
      if (eq(wij,0)){
        return(matrix(0,p,p))
      }else{
      dmuij <- as.numeric(x[(7+p):(6+2*p)])
      indij1 <- as.numeric(x["ind1"])
      po <- as.numeric(x["pair_order"])
      IDi <- as.integer(x["IDi"])
      IDj <- as.integer(x["IDj"])
      Hij1 <- H[po,]*survZ[IDi,indij1]*survZ[IDj,indij1]
      return(-matrix(Hij1,p,1)%*%matrix(dmuij,1,p))
      }
    }
    dEE <- array(apply(dataM.beta,1,dEE.ij),c(p,p,m))
  }else{
    w <- dataM.beta[,"w"]
    mu <- dataM.beta[,"mu"]
    Zd <- as.matrix(dataM.beta[,6:(5+p)],m,p)
    dmu <- as.matrix(dataM.beta[,(7+p):(6+2*p)],m,p)
    EE <- as.matrix(H*abs(w)*(eq(w,1)-mu)-w*(t2<tau)*H,m,p)
    dEE.input <- cbind(H,w,dmu)
    dEE <- array(apply(dEE.input,1,int.deriv.ij),c(p,p,m))
  }
  return(list(EE=EE,dEE=dEE))
}


NRfunPW <- function(beta,H,t,rho,dataf,tau,survZ=NULL,eps=1e-6,maxiter=50){
  error <- 1
  iter <- 1
  conv <- FALSE
  while(error>eps&&iter<maxiter){
    obj <- UEE.fun(beta,H,t,rho,dataf,tau,survZ=survZ)
    EE <- obj$EE
    dEE <- obj$dEE
    UEE <- colMeans(EE)
    dUEE <- apply(dEE,1:2,mean)
    inv.dUEE <- solve(dUEE)
    inc <- -inv.dUEE%*%UEE
    beta <- beta+inc
    error <- sum(abs(inc))
    iter <- iter+1
  }
  if(iter<maxiter)(conv <- TRUE)
  return(list(beta=beta,EE=EE,UEE=UEE,dEE=dEE,dUEE=dUEE,conv=conv,iter=iter))
}


otimes2 <- function(y){
  return(outer(y,y))
}


IF.fun <- function(beta,H,t,rho,dataf,tau,survZ){
    n <- dataf[nrow(dataf),2]
    obj <- UEE.fun(beta,H,t,rho,dataf,tau,survZ)
    EE <- obj$EE
    dEE <- obj$dEE
    UEE <- colMeans(EE)
    dUEE <- apply(dEE,1:2,mean)
    inv.dUEE <- solve(dUEE)
    p <- ncol(EE)
    IDi <- dataf[,"IDi"]
    IDj <- dataf[,"IDj"]
    IF.fun.id <- function(id){
      EE.id <- as.matrix(EE[IDi==id|IDj==id,],n-1,p)
      return(colMeans(EE.id))
    }

    IFbeta <- -2*inv.dUEE%*%matrix(sapply(1:n,IF.fun.id),p,n)
    S_array <- array(apply(IFbeta,2,otimes2),c(p,p,n))
    S <- apply(S_array,1:2,mean)*n/(n-p-1)
    Vscore <- dUEE%*%S%*%t(dUEE)
    return(list(IFbeta=IFbeta,S=S,Vscore=Vscore))
}


## Programs for stratified PW models

Hfun <- function(data,dataf,datap,tau,h){
  # Calculate the weight function H(t|Z_i,Z_j) = (Z_i-Z_j)*h_l(t|Z_i,Z_j)
  # Here we only consider h_l(t|Z_i,Z_j) = h_l. For example h_l=1 for all l=1,...,L
  m <- nrow(dataf)
  p <- ncol(data)-3
  Zd <- as.matrix(dataf[,6:(5+p)],m,p)
  t <- sort(unique(datap[, "timeD"]))
  t <- t[t <= tau]
  H <- Zd * h
  return(list(H=H,t=t))
}


s.UEE.fun <- function(beta,H,t,dataf,tau){
  # Calculate estimating equation and its derivative.
  p <- length(beta)
  m <- nrow(dataf)
  dataM.beta <- M.fun(beta,dataf,tau)

  t1 <- dataM.beta[,"t1"]
  t2 <- dataM.beta[,"t2"]

  w <- dataM.beta[,"w"]
  mu <- dataM.beta[,"mu"]
  Zd <- as.matrix(dataM.beta[,6:(5+p)],m,p)
  dmu <- as.matrix(dataM.beta[,(7+p):(6+2*p)],m,p)
  EE <- as.matrix(H*abs(w)*(eq(w,1)-mu)-w*(t2<tau)*H,m,p)
  dEE.input <- cbind(H,w,dmu)
  dEE <- array(apply(dEE.input,1,int.deriv.ij),c(p,p,m))

  return(list(EE=EE,dEE=dEE))
}

NRfunPWS <- function(beta,H,t,dataf,tau,eps=1e-6,maxiter=50){
  error <- 1
  iter <- 1
  conv <- FALSE
  ns <- length(tau) # number of strata
  p <- length(beta) # dimension of covariates
  while(error>eps && iter<maxiter){
    UEE <- rep(0, p)
    dUEE <- matrix(0,p,p)
    for(i in c(1:ns)){
      obj <- s.UEE.fun(beta,H[[i]],t[[i]],dataf[[i]],tau[[i]])
      EE <- obj$EE
      dEE <- obj$dEE
      UEE <- UEE + colMeans(EE)
      dUEE <- dUEE + apply(dEE,1:2,mean)
    }
    inv.dUEE <- solve(dUEE)
    inc <- -inv.dUEE%*%UEE
    beta <- beta+inc
    error <- sum(abs(inc))
    iter <- iter+1
  }
  if(iter < maxiter)(conv <- TRUE)
  return(list(beta=beta,EE=EE,UEE=UEE,dEE=dEE,dUEE=dUEE,conv=conv,iter=iter))
}


s.IF.fun <- function(beta,H,t,dataf,tau,fixedL){
  ns <- length(tau)
  p <- length(beta)
  EE <- NULL
  UEE <- rep(0, p)
  dUEE <- matrix(0,p,p)
  USL <- matrix(0,p,p)
  for(l in c(1:ns)){
    obj <- s.UEE.fun(beta,H[[l]],t[[l]],dataf[[l]],tau[[l]])
    EE[[l]] <- obj$EE
    dEE <- obj$dEE
    UEE <- UEE + colMeans(EE[[l]])
    dUEE <- dUEE + apply(dEE,1:2,mean)
    if(!fixedL){
      USL <- USL + otimes2(colMeans(EE[[l]]))
    }
  }
  inv.dUEE <- solve(dUEE)

  if(!fixedL){
    Var <- inv.dUEE %*% USL %*% inv.dUEE
    V <- USL
  }else{
    Var <- matrix(0,p,p)
    V <- matrix(0, p, p)
    for(l in c(1:ns)){
      nl <- dataf[[l]][nrow(dataf[[l]]),2]
      IDi <- dataf[[l]][,"IDi"]
      IDj <- dataf[[l]][,"IDj"]
      EEl <- EE[[l]]
      IF.fun.id <- function(id){
        EE.id <- as.matrix(EEl[IDi==id|IDj==id,],nl-1,p)
        return(colMeans(EE.id))
      }
      IFbeta <- -2*inv.dUEE %*% matrix(sapply(1:nl,IF.fun.id),p,nl)
      S_array <- array(apply(IFbeta,2,otimes2), c(p,p,nl))
      #S <- apply(S_array,1:2,mean)*nl/(nl-p-1)
      S <- apply(S_array,1:2,mean)
      V <- V + dUEE %*% S %*% t(dUEE)/nl
      Var <- Var + S/nl
    }
  }
  return(list(Var=Var, V=V))
}

#' Computes the standardized score processes
#' @description Computes the standarized score processes for the covariates.
#' @param obj an object of class pwreg.
#' @param t a vector containing times. If not specified, the function will use
#' all unique event times from the data.
#' @return An object of class \code{pwreg.score} consisting of \code{t:}
#' a vector of times; and \code{score:} a matrix whose rows are the standardized score processes
#' as a function of \code{t}.
#' @seealso \code{\link{pwreg}}, \code{\link{print.pwreg}}
#' @export
#' @keywords pwreg
#' @references Mao, L. and Wang, T. (2020). A class of proportional win-fractions
#' regression models for composite outcomes. Biometrics, 10.1111/biom.13382
#' @examples
#' library(WR)
#' head(non_ischemic)
#'
#' # Randomly sample 200 subjects from non_ischemic data
#' id_unique <-unique(non_ischemic$ID)
#' set.seed(2019)
#' id_sample <- sample(id_unique, 200)
#' non_ischemic_reduce <- non_ischemic[non_ischemic$ID %in% id_sample, ]
#'
#' # Use the reduced non_ischemic data for analysis
#' nr <- nrow(non_ischemic_reduce)
#' p <- ncol(non_ischemic_reduce)-3
#' ID <- non_ischemic_reduce[,"ID"]
#' time <- non_ischemic_reduce[,"time"]
#' status <- non_ischemic_reduce[,"status"]
#' Z <- as.matrix(non_ischemic_reduce[,4:(3+p)],nr,p)
#' pwreg.obj <- pwreg(time=time,status=status,Z=Z,ID=ID)
#' score.obj <- score.proc(pwreg.obj)
#' #plot the standardized score process for the first covariate
#' plot(score.obj, k = 1)
score.proc <- function(obj,t=NULL){
  if(!is.null(obj$strata)){
    #stop("score process can only be calculated for PW model without stratification")
    sv <- levels(obj$strata)
    beta <- obj$beta
    varnames <- obj$varnames
    p <- length(beta)
    V <- obj$V

    t <- unlist(obj$t)
    tau <- max(unlist(obj$tau))
    t <- c(sort(t), tau)
    t <- t[!duplicated(t)]
    l <- length(t)

    std.score <- matrix(0, p, l)
    for (i in c(1:length(sv))) {
      H.i <- obj$H[[i]]
      tau.i <- obj$tau[[i]]
      dataf.i <- obj$dataf[[i]]
      m <- nrow(dataf.i)
      dataM.beta <- M.fun(beta,dataf.i,tau.i)
      t1 <- dataM.beta[,"t1"]
      t2 <- dataM.beta[,"t2"]
      w <- dataM.beta[,"w"]
      mu <- dataM.beta[,"mu"]
      score.fun <- function(t0){
        if(t0 >= tau.i){t0 = tau.i}
        EE <- as.matrix((t1<t0)*H.i*abs(w)*(eq(w,1)-mu)-w*(t2<t0)*H.i,m,p)
        return(colMeans(EE))
      }
      scorep <- as.matrix(sapply(t,score.fun),p,l)
      if (ncol(scorep)==1){
        scorep <- t(scorep)
      }
      #zz <-  sqrt(solve(diag(diag(V),p)))%*%scorep
      std.score <- std.score + sqrt(solve(diag(diag(V),p)))%*%scorep
    }
  }

  if(is.null(obj$strata)){
    beta <- obj$beta
    H <- obj$H
    if (is.null(t)){
      t <- obj$t
    }
    varnames <- obj$varnames
    dataf <- obj$dataf
    tau <- obj$tau
    V <- obj$V
    p <- length(beta)
    m <- nrow(dataf)
    dataM.beta <- M.fun(beta,dataf,tau)
    t <- sort(c(t,tau))
    l <- length(t)
    t1 <- dataM.beta[,"t1"]
    t2 <- dataM.beta[,"t2"]
    w <- dataM.beta[,"w"]
    mu <- dataM.beta[,"mu"]
    score.fun <- function(t0){
      if(t0 >= tau){t0 = tau}
      EE <- as.matrix((t1<t0)*H*abs(w)*(eq(w,1)-mu)-w*(t2<t0)*(t0 <= tau)*H,m,p)
      return(colMeans(EE))
    }
    scorep <- as.matrix(sapply(t,score.fun),p,l)
    if (ncol(scorep)==1){
      scorep=t(scorep)
    }
    std.score <- sqrt(solve(diag(diag(V),p)))%*%scorep
  }

  rownames(std.score) <- varnames
  score.obj <- list(score=std.score,t=t)
  class(score.obj) <- "pwreg.score"
  return(score.obj)
}


#' Fit a standard proportional win-fractions (PW) regression model
#' @description Fit a standard proportional win-fractions (PW) regression model.
#' @param time a vector of event times.
#' @param status a vector of event type labels. 0: censoring, 1:death
#' and 2: non-fatal event.
#' @param Z a matrix or a vector of covariates.
#' @param ID a vector of unique subject-level identifiers.
#' @param rho a non-negative number as the power of the survival function used
#' in the weight. Default (\code{rho=0}) is recommended. If there is a `strata` argument,
#' then `rho` is ignored.
#' @param strata a vector of stratifying variable if a stratified model is desired.
#' @param fixedL logical variable indicating which variance estimator to be used. If `TRUE`,
#' the type I variance estimator (for a small number strata) is used; otherwise
#' the type II variance estimator (for a large number strata) is used.
#' @param eps precision for the convergence of Newton-Raphson algorithm.
#' @param maxiter maximum number of iterations allow for the Newton-Raphson
#' algorithm.
#' @return An object of class \code{pwreg} with the following components.
#' \code{beta}:a vector of estimated regression coefficients. \code{Var}:estimated
#' covariance matrix for \code{beta}. \code{conv:} boolean variable indicating
#' whether the algorithm converged within the maximum number of iterations.
#' @seealso \code{\link{score.proc}}, \code{\link{print.pwreg}}
#' @export
#' @import survival
#' @importFrom utils combn
#' @importFrom stats pchisq pnorm printCoefmat qnorm
#' @keywords pwreg
#' @references Mao, L. and Wang, T. (2020). A class of proportional win-fractions
#' regression models for composite outcomes. Biometrics, 10.1111/biom.13382
#' @references Wang, T. and Mao, L. (2021+). Stratified Proportional Win-fractions
#' Regression Analysis.
#' @examples
#' library(WR)
#' head(non_ischemic)
#' id_unique <-unique(non_ischemic$ID)
#'
#' # Randomly sample 200 subjects from non_ischemic data
#' set.seed(2019)
#' id_sample <- sample(id_unique, 200)
#' non_ischemic_reduce <- non_ischemic[non_ischemic$ID %in% id_sample, ]
#'
#' # Use the reduced non_ischemic data for analysis
#' nr <- nrow(non_ischemic_reduce)
#' p <- ncol(non_ischemic_reduce)-3
#' ID <- non_ischemic_reduce[,"ID"]
#' time <- non_ischemic_reduce[,"time"]
#' status <- non_ischemic_reduce[,"status"]
#' Z <- as.matrix(non_ischemic_reduce[,4:(3+p)],nr,p)
#' ## unstratified analysis
#' pwreg.obj <- pwreg(time=time,status=status,Z=Z,ID=ID)
#' print(pwreg.obj)
#' \dontrun{
#' ## stratified PW by sex
#' sex<-Z[,3]
#' ## take out sex from the covariate matrix
#' Z1<-Z[,-3]
#' pwreg.obj1 <- pwreg(time=time,status=status,Z=Z1,ID=ID,strata=sex)
#' print(pwreg.obj1)
#' }
 pwreg <- function(ID,time,status,Z,rho=0,strata=NULL,fixedL=TRUE,eps=1e-4, maxiter=50){

  ID.save <- ID


  if(is.character(ID) | is.factor(ID)){
    ID <- as.numeric(factor(ID,levels = unique(ID)))
  }

  if(is.null(strata)){
    data <- cbind(ID,time,status,Z)
    obj <- extract.times(data)
    datap <- obj$datap
    n <- nrow(datap)
    tau <- obj$tau
    p <- ncol(datap)-4
    outcome <- delta.data(datap,tau)
    q <- ncol(outcome)
    dataf <- outcome[,-q]
    comp <- outcome[,q]
    beta <- rep(0,p)
    H.obj <- surv.Hfun(rho,data,dataf,datap,tau)
    H <- H.obj$H
    t <- H.obj$t
    survZ <- H.obj$survZ
    obj <- NRfunPW(beta,H,t,rho,dataf,tau,survZ=survZ,eps=eps,maxiter=maxiter)

    beta <- obj$beta
    conv <- obj$conv
    iter <- obj$iter

    if(conv){
      IF.obj <- IF.fun(beta,H,t,rho,dataf,tau,survZ)
      IFbeta <- IF.obj$IFbeta
      S <- IF.obj$S
      Vscore <- IF.obj$Vscore
    }else{
      message("Newton-Raphson did not converge! \n You may try increasing maximum
      number of iterations by specifying 'maxiter='.")
    }
    obj <- list(beta=beta,Var=S/n,IFbeta=IFbeta,conv=conv,V=Vscore/n,n=n,
                t=t,comp=comp,dataf=dataf,H=H,tau=tau,strata=strata,
                varnames=colnames(data)[4:(3+p)],i=iter,call=match.call())
  }else{
    if(!is.factor(strata)){
      strata <- factor(strata)
    }

    data <- cbind(ID,time,status,Z)
    sv <- levels(strata)
    ns <- length(sv)
    p <- ncol(data) - 3

    tmp <- extract.times(cbind(data,strata))
    h <- as.numeric(table(tmp$datap[,'strata']))/nrow(tmp$datap)
    rm(tmp)

    datap <- list()
    dataf <- list()
    tau <- list()
    H <- list()
    t <- list()
    for(i in c(1:ns)){
      obj <- extract.times(data[strata==sv[i],])
      datap[[i]] <- obj$datap
      tau[[i]] <- obj$tau
      outcome <- delta.data(obj$datap, obj$tau)
      q <- ncol(outcome)
      dataf[[i]] <- outcome[,-q]
      H.obj <- Hfun(data,dataf=dataf[[i]],datap=datap[[i]],tau=tau[[i]],h=h[i])
      H[[i]] <- H.obj$H
      t[[i]] <- H.obj$t
    }

    beta <- rep(0,p)
    obj <- NRfunPWS(beta, H, t, dataf, tau, eps=eps, maxiter=maxiter)

    beta <- obj$beta
    conv <- obj$conv
    iter <- obj$iter

    if(conv){
      IF.obj <- s.IF.fun(beta=obj$beta,H,t,dataf,tau,fixedL)
      Var <- IF.obj$Var
      V <- IF.obj$V
    }else{
      message("Newton-Raphson did not converge! \n You may try increasing maximum
      number of iterations by specifying 'maxiter='.")
    }
    obj <- list(
      beta=beta, conv=conv, Var=Var, V=V, fixedL=fixedL, h=h, strata=strata,
      t = t, dataf=dataf, H=H, tau=tau, varnames = colnames(data)[4:(3+p)],
      i=iter, call=match.call()
    )
  }

  class(obj) <- "pwreg"
  return(obj)
}

