#' Print the results of the proportional win-fractions regression model
#' @description Print the results of the proportional win-fractions regression model.
#' @param x an object of class \code{pwreg}.
#' @param ... further arguments passed to or from other methods.
#' @return Print the results of \code{pwreg} object
#' @seealso \code{\link{pwreg}}
#' @export
#' @keywords pwreg
#' @examples
#' # see the example for pwreg
print.pwreg <- function(x,...){

  cat("Call:\n")
  print(x$call)
  cat("\n")
  strata <- x$strata
  beta <- x$beta
  var <- x$Var
  t <- x$t
  n <- x$n
  comp <- x$comp
  varnames <- x$varnames
  se <- sqrt(diag(var))
  pv <- 2*(1-pnorm(abs(beta/se)))

  if(is.null(strata)){
    cat("Proportional win-fractions regression models for priority-adjusted composite endpoint\n
    (Mao and Wang, 2020):\n\n")
    tn <- n*(n-1)/2
    cat("Total number of pairs:", tn,"\n")
    cat("Wins-losses on death: ",paste(sum(comp==1)," (",round(100*sum(comp==1)/tn,1),"%)",sep=""),"\n")
    cat("Wins-losses on non-fatal event: ",paste(sum(comp==2)," (",round(100*sum(comp==2)/tn,1),"%)",sep=""),"\n")
    cat("Indeterminate pairs",paste(sum(comp==0)," (",round(100*sum(comp==0)/tn,1),"%)",sep=""),"\n\n")
  }else{
    cat("Stratified proportional win-fractions regression analysis\n
    (Wang and Mao, 2021+):\n\n")
    cat("Total number of strata:", length(levels(strata)), "\n")
  }

  cat("Newton-Raphson algorithm converged in", x$i,"iterations.\n\n")
  Stat <- t(beta)%*%solve(var)%*%beta
  pval <- 1-pchisq(Stat,length(beta))
  cat("Overall test: chisq test with",length(beta),"degrees of freedom;","\n",
      "Wald statistic", round(Stat,1),"with p-value", pval,"\n\n")

  table <- cbind(
    Estimate = beta,
    StdErr = se,
    z.value = beta/se,
    p.value = pv)

  colnames(table) <- c("Estimate","se","z.value","p.value")
  rownames(table) <- varnames

  cat("Estimates for Regression parameters:\n")
  printCoefmat(table, P.values=TRUE, has.Pvalue=TRUE)
  cat("\n")
  cat("\n")

  za <- qnorm(0.975)
  MR <- cbind(exp(beta),exp(beta-za*se),exp(beta+za*se))
  colnames(MR) <- c("Win Ratio","95% lower CL","95% higher CL")
  rownames(MR) <- varnames
  cat("Point and interval estimates for the win ratios:\n")
  print(MR)
  cat("\n")
}

#' Plot the standardized score processes
#' @description Plot the standardized score processes.
#' @param x an object of class \code{pwreg.score}.
#' @param k A positive integer indicating the order of covariate to be plotted. For example, \code{k=3} requests the standardized score process for the third covariate in the covariate matrix \code{Z}.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param lty the line type. Default is 1.
#' @param frame.plot a logical variable indicating if a frame should be drawn in the 1D case.
#' @param add a logical variable indicating whether add to current plot?
#' @param ylim a vector indicating the range of y-axis. Default is (-3,3).
#' @param xlim a vector indicating the range of x-axis. Default is NULL.
#' @param lwd the line width, a positive number. Default is 1.
#' @param ... further arguments passed to or from other methods.
#' @return A plot of the standardized score process for object \code{pwreg.score}.
#' @seealso \code{\link{score.proc}}
#' @export
#' @importFrom graphics lines plot
#' @keywords pwreg
#' @examples
#' # see the example for score.proc
plot.pwreg.score <- function(
  x,k,xlab="Time",ylab="Standardized score",lty=1,frame.plot=TRUE,add=FALSE,
  ylim=c(-3,3),xlim=NULL,lwd=1,...){
  score <- x$score
  scorek <- score[k,]
  t <- x$t
  if (is.null(xlim)){
    xlim <- c(0,max(t))
  }

  if (add==FALSE){
  plot(t,scorek,type='l',
    lty=lty,xlab=xlab,ylab=ylab,ylim=ylim,xlim=xlim,lwd=lwd,frame.plot=frame.plot,
    main=rownames(score)[k],...)
  }else{
    lines(t,scorek,lty=lty,...)
  }
}


#' Print information on the content of the pwreg.score object
#' @description Print information on the content of the pwreg.score object
#' @param x A object of class pwreg.score.
#' @param ... further arguments passed to or from other methods.
#' @return Print the results of \code{pwreg.score} object.
#' @seealso \code{\link{score.proc}}
#' @export
#' @keywords pwreg
#' @examples
#' # see the example for score.proc
print.pwreg.score <- function(x,...){
cat("This object contains two components:\n")
 cat(" 't': an l-vector of times\n")
 cat(" 'score': a p-by-l matrix whose k'th row is the standardized score process for the k'th covariate
          as a function of t\n\n")
 cat("Use 'plot(object,k=k)' to plot the k'th score process.\n")
}


#####################Recurrent event WR ######################
#' Print the results of the two-sample recurrent-event win ratio analysis
#' @description Print the results of the two-sample recurrent-event win ratio analysis.
#' @param x an object of class \code{WRrec}.
#' @param ... further arguments passed to or from other methods.
#' @return Print the results of \code{WRrec} object.
#' @seealso \code{\link{WRrec}}
#' @export
#' @keywords WRrec
#' @examples
#' # see the example for WRrec
print.WRrec <- function(x,...){

  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(x$desc)
  n0<-x$desc[1,1]
  n1<-x$desc[2,1]
  N<-n0*n1


  ## print basic results for LWR (only for unstratified case) ##
  # cat("By last-event-assisted win ratio (LWR):\n\n")
  # Nwin<-sum(x$R.mat==1)
  # Wp<-Nwin/N
  # Nloss<-sum(x$R.mat==-1)
  # Lp<-Nloss/N
  # cat(paste0("Among ",n1," x ",n0," = ",N," pairs, "))
  # cat(paste0(Nwin, " (",round(100*Wp,1),"%) wins and ",Nloss, " (",round(100*Lp,1),"%) losses\n"))

  # Output table
  tab<-matrix(NA,1,4)
  colnames(tab)<-c("Win prob","Loss prob","WR (95% CI)*","p-value")
  rownames(tab)<-"LWR"
  za<-qnorm(0.975)
  beta<-x$log.WR
  se<-x$se
  theta<-x$theta


  tab[1,]<-c(paste0(round(100*theta[1],1),"%"),
       paste0(round(100*theta[2],1),"%"),
       paste0(round(exp(beta),2)," (",round(exp(beta-za*se),2),", ",round(exp(beta+za*se),2),")"),
       round(1-pchisq((beta/se)^2,1),4))

  cat("\n")

  nf<-!is.null(x$log.WR.naive)

  if (!nf){
  cat("Analysis of last-event-assisted WR (LWR):\n")
  }else{
    cat("Analysis of last-event-assisted WR (LWR; recommended), first-event-assisted WR (FWR), and naive WR (NWR):\n")
}


  if (!is.null(x$log.WR.naive)){
    beta.naive<-x$log.WR.naive
    se.naive<-x$se.naive
    theta.naive<-x$theta.naive

    beta.FI<-x$log.WR.FI
    se.FI<-x$se.FI
    theta.FI<-x$theta.FI

    tab<-rbind(tab,
               c(paste0(round(100*theta.FI[1],1),"%"),
                 paste0(round(100*theta.FI[2],1),"%"),
                 paste0(round(exp(beta.FI),2)," (",round(exp(beta.FI-za*se.FI),2),", ",
                        round(exp(beta.FI+za*se.FI),2),")"),
                 round(1-pchisq((beta.FI/se.FI)^2,1),4)),
               c(paste0(round(100*theta.naive[1],1),"%"),
               paste0(round(100*theta.naive[2],1),"%"),
               paste0(round(exp(beta.naive),2)," (",round(exp(beta.naive-za*se.naive),2),", ",
                      round(exp(beta.naive+za*se.naive),2),")"),
               round(1-pchisq((beta.naive/se.naive)^2,1),4))
               )
  rownames(tab)[2:3]<-c("FWR","NWR")
  }

  print(tab,quote=F)
  cat("-----\n")
  cat("*Note: The scale of WR should be interpreted with caution as it depends on \n")
  cat("censoring distribution without modeling assumptions.")
}




