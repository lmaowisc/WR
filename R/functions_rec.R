


# ID=1,2,3,...n; time; status= 0 (censoring), 1 (death), 2 (hospitalization); trt=1, 0
#
#
#
# Re-organize the data into two lists by arm: trt=1, 0
# scalers: ID, status
# (N_i+1)-vectors: time, (N_i is the number of hosps)
process.dat <- function(ID,time,status,trt){
   o=order(ID,time)
   ID=ID[o]
   time=time[o]
   status=status[o]
   trt=trt[o]


  time1=time[trt==1]
  status1=status[trt==1]
  ID1=ID[trt==1]

  time0=time[trt==0]
  status0=status[trt==0]
  ID0=ID[trt==0]


  uid1=unique(ID1)
  n1=length(uid1)

  uid0=unique(ID0)
  n0=length(uid0)



  dat.list1=vector("list", n1)
  dat.list0=vector("list", n0)

  names(dat.list1)=uid1
  names(dat.list0)=uid0

  for (i in 1:n1){
    id=uid1[i]
    ind=(ID1==id)
    dat.list1[[i]]=list(ID=id,time=time1[ind],status=tail(status1[ind],n=1))
  }

  for (i in 1:n0){
    id=uid0[i]
    ind=(ID0==id)
    dat.list0[[i]]=list(ID=id,time=time0[ind],status=tail(status0[ind],n=1))
  }


return(list(dat1=dat.list1,dat0=dat.list0))

}




#Define the win function
#Arguments: two lists sub.i and sub.j (ID,time,status)
#Each list is an element in the lists returned by "process.dat"
#Return: a list of three components, each from a win function
# 1 if sub.i wins sub.j; -1 if vice versa; 0 if otherwise
win.fun.Rec<-function(sub.i,sub.j){
  #minimum of follow-up time
  max.time.i=max(sub.i$time)
  max.time.j=max(sub.j$time)

  X=min(max.time.i,max.time.j)
  status.i=sub.i$status
  status.j=sub.j$status

  ind.i=(sub.i$time<=X)
  ind.j=(sub.j$time<=X)

  time.i=sub.i$time
  time.j=sub.j$time



  if (max.time.i>X){
    time.i=c(time.i[ind.i],X)
    status.i=0
  }
  if (max.time.j>X){
    time.j=c(time.j[ind.j],X)
    status.j=0
  }




  #naive comparison result
  naive=status.j-status.i
  if (naive==0){
    naive=(length(time.j)>length(time.i))-(length(time.i)>length(time.j))
  }
  # regular and first-time-imputed
  if (naive==0&&length(time.i)>1){
    FI=(time.j[1]<time.i[1])-(time.i[1]<time.j[1])
    p=length(time.i)
    R=(time.j[p-1]<time.i[p-1])-(time.i[p-1]<time.j[p-1])
  }else{
    FI=naive
    R=naive
  }
  result=list(naive=naive,FI=FI,R=R)
}



#Compute the win-loss indicators from the pariwise
# comparisons between the two arms.
# Arguments: two lists (treatment and control) returned by process.dat()
#   naive=F (default): regular only; naive=T: regular, naive, and first-time-imputed
# Return: an n1 x n0 matrix from each win function
win.loss.matrix=function(dat1,dat0,naive=FALSE){
  n1=length(dat1)
  n0=length(dat0)
  R.matrix=matrix(NA,n1,n0)

  names1=names(dat1)
  names0=names(dat0)


  rownames(R.matrix)=names1
  colnames(R.matrix)=names0

  if (naive){
    naive.matrix=matrix(NA,n1,n0)
    FI.matrix=matrix(NA,n1,n0)
    rownames(naive.matrix)=names1
    colnames(naive.matrix)=names0
    rownames(FI.matrix)=names1
    colnames(FI.matrix)=names0
  }else{
    naive.matrix=NULL
    FI.matrix=NULL
  }
  for (i in 1:n1){
    sub.i=dat1[[i]]
      for (j in 1:n0){
        sub.j=dat0[[j]]
        obj.ij=win.fun.Rec(sub.i,sub.j)
        R.matrix[i,j]=obj.ij$R
        if (naive){
          naive.matrix[i,j]=obj.ij$naive
          FI.matrix[i,j]=obj.ij$FI
        }
      }
  }

  return(list(R.mat=R.matrix,naive.mat=naive.matrix,FI.mat=FI.matrix))
}

############################Cross-comparisons Ends################



#Calculate the statistic and variance
#Argument: R.mat, naive.mat, or FI.mat returned by win.loss.matrix
#Output: theta= 2-vector of win-loss probabilities; Sigma=(2x2)-covariance matrix for theta
wr.stat.se=function(mat){

  n1=nrow(mat)
  n0=ncol(mat)


  win.mat=(mat==1)
  loss.mat=(mat==-1)
  #win-loss probabilities
  theta1=mean(win.mat)
  theta0=mean(loss.mat)
  theta=c(theta1,theta0)

  #influence functions for win-loss probabilities
  #of treatment and control arms
  w1=cbind(rowMeans(win.mat)-theta1,rowMeans(loss.mat)-theta0)#n1 x 2
  w0=cbind(colMeans(win.mat)-theta1,colMeans(loss.mat)-theta0)#n0 x 2
  #colMeans(w1): should be c(0,0)

  #variance of the win-loss estimators
  Sigma=n1^{-2}*t(w1)%*%w1+n0^{-2}*t(w0)%*%w0

  return(list(theta=theta,Sigma=Sigma))



}


##############################
#      Main function         #
##############################


#' Generalized win ratio tests
#'
#' Perform stratified two-sample test of possibly recurrent nonfatal
#' event and death using the recommended last-event assisted win ratio (LWR), and/or
#' naive win ratio (NWR) and first-event assisted win ratio (FWR) (Mao et al., 2022).
#' The LWR and FWR reduce to the standard win ratio of Pocock et al. (2012).
#'
#' @param ID A vector of unique patient identifiers.
#' @param time A numeric vector of event times.
#' @param status A vector of event type variable; 2 = recurrent event, 1 = death,
#' and 0 = censoring.
#' @param trt A vector of binary treatment indicators.
#' @param strata A vector of categorical variable for strata; Default is NULL,
#' which leads to unstratified analysis.
#' @param naive If TRUE, results for NWR and FWR will be provided in addition
#' to LWR; Default is FALSE, which gives LWR only.
#' @return An object of class \code{WRrec}, which contains the following
#' elements.
#' \item{theta}{A bivariate vector of win/loss fractions by LWR.}
#' \item{log.WR, se}{Log-win ratio estimate and its standard error by LWR.}
#' \item{pval}{\eqn{p}-value by the LWR test.}
#' \item{theta.naive}{A bivariate vector of win/loss fractions by NWR.}
#' \item{log.WR.naive, se.naive}{Log-win ratio estimate and its standard error by NWR.}
#' \item{theta.FI}{A bivariate vector of win/loss fractions by FWR.}
#' \item{log.WR.FI, se.FI}{Log-win ratio estimate and its standard error by FWR.}
#' \item{...}{}
#'
#' @references Mao, L., Kim, K. and Li, Y. (2022). On recurrent-event win ratio.
#' Statistical Methods in Medical Research, under review.
#' @references Pocock, S., Ariti, C., Collier, T., and Wang, D. (2012). The win ratio: a new approach
#' to the analysis of composite endpoints in clinical trials based on clinical priorities.
#'  European Heart Journal, 33, 176--182.
#' @examples
#' ## load the HF-ACTION trial data
#' library(WR)
#' head(hfaction_cpx9)
#' dat<-hfaction_cpx9
#' ## Comparing exercise training to usual care by LWR, FWR, and NWR
#' obj<-WRrec(ID=dat$patid,time=dat$time,status=dat$status,
#'           trt=dat$trt_ab,strata=dat$age60,naive=TRUE)
#' ## print the results
#' obj
#' @keywords WRrec
#' @importFrom stats median
#' @export
#' @aliases WRrec
#' @seealso \code{\link{print.WRrec}}.
WRrec=function(ID,time,status,trt,strata=NULL,naive=FALSE){
  n=length(unique(ID))

  #### descriptive statistics
  n1<-length(unique(ID[trt==1]))
  n0<-length(unique(ID[trt==0]))
  Nrec1<-sum(status[trt==1]==2)
  Nrec0<-sum(status[trt==0]==2)
  Ndeath1<-sum(status[trt==1]==1)
  Ndeath0<-sum(status[trt==0]==1)
  X<-time[!status==2]
  MedFU1<-median(X[trt[!status==2]==1])
  MedFU0<-median(X[trt[!status==2]==0])

  desc<-rbind(c(n0,Nrec0,Ndeath0,MedFU0),
              c(n1,Nrec1,Ndeath1,MedFU1)
              )
  colnames(desc)<-c("N","Rec. Event","Death","Med. Follow-up")
  rownames(desc)<-c("Control","Treatment")

  ## formal analysis ##
  if (!is.null(strata)){
    stratum.list=unique(strata)
    K=length(stratum.list)
  }else{
    strata=rep(1,length(ID))
    stratum.list=1
    K=1
  }

  weights=rep(NA,K)

  # Theta: (2xK)-matrix
  # the kth column of Theta contains the win-loss probabilities for the kth stratum
  Theta=matrix(NA,2,K)
  #Sigma: (2 x 2)-matrix
  Sigma=matrix(0,2,2)

  Theta.naive=matrix(NA,2,K)
  Sigma.naive=matrix(0,2,2)


  Theta.FI=matrix(NA,2,K)
  Sigma.FI=matrix(0,2,2)


  #compute the
  for (k in 1:K){

    ind=(strata==stratum.list[k])

    #compute the win-loss probabilities and their covariance
    #matrix for the kth stratum
    obj=process.dat(ID[ind],time[ind],status[ind],trt[ind])
    dat1=obj$dat1
    dat0=obj$dat0
    obj1=win.loss.matrix(dat1,dat0,naive=naive)
     R.mat=obj1$R.mat
     if (naive){
     naive.mat=obj1$naive.mat
     FI.mat=obj1$FI.mat
     }
     #weights=stratum-specific proportions of sample size
     weights[k]=sum(dim(R.mat))/n

     obj2=wr.stat.se(mat=R.mat)
     Theta[,k]=obj2$theta
     Sigma=Sigma+weights[k]^2*obj2$Sigma

     if (naive){
       obj2.naive=wr.stat.se(mat=naive.mat)
       Theta.naive[,k]=obj2.naive$theta
       Sigma.naive=Sigma.naive+weights[k]^2*obj2.naive$Sigma

       obj2.FI=wr.stat.se(mat=FI.mat)
       Theta.FI[,k]=obj2.FI$theta
       Sigma.FI=Sigma.FI+weights[k]^2*obj2.FI$Sigma
     }

  }

     #compute the overall log-win ratio and its
     #standard error
     theta=rowSums(Theta*rbind(weights,weights))
     #derivative of log-WR
     f=c(1/theta[1],-1/theta[2])
     se=sqrt(t(f)%*%Sigma%*%f)
     log.WR=log(theta[1]/theta[2])
      pval<-2*(1-pnorm(abs(log.WR/se)))

     if (naive){
       theta.naive=rowSums(Theta.naive*rbind(weights,weights))
       f.naive=c(1/theta.naive[1],-1/theta.naive[2])
       se.naive=sqrt(t(f.naive)%*%Sigma.naive%*%f.naive)
       log.WR.naive=log(theta.naive[1]/theta.naive[2])


       theta.FI=rowSums(Theta.FI*rbind(weights,weights))
       f.FI=c(1/theta.FI[1],-1/theta.FI[2])
       se.FI=sqrt(t(f.FI)%*%Sigma.FI%*%f.FI)
       log.WR.FI=log(theta.FI[1]/theta.FI[2])
     }else{
       log.WR.naive=NULL; se.naive=NULL;theta.naive=NULL;naive.mat=NULL
       log.WR.FI=NULL; se.FI=NULL;theta.FI=NULL;FI.mat=NULL
     }


     result<-list(log.WR=log.WR,se=se,pval=pval,theta=theta,R.mat=R.mat,
                  log.WR.naive=log.WR.naive,se.naive=se.naive,theta.naive=theta.naive,naive.mat=naive.mat,
                  log.WR.FI=log.WR.FI,se.FI=se.FI,theta.FI=theta.FI,FI.mat=FI.mat,
                  call=match.call(),desc=desc)
     class(result)<-"WRrec"

     return(result)


  }










