######################################################################################
#                                                                                    #
#                     Final Size estimation of small populations                     #
#                     Methods derived from Kroese and De Jong, 2001                  #
#                                                                                    #
#                     Copyright R-code: Egil A.J.Fischer                             #
#                                       e.a.j.fischer@uu.nl/egil@egilfischer.nl      #
######################################################################################

#######################################################################################################
# Functions to calculate the exact probability of a specific final size given inital conditions and R #
#######################################################################################################
library("compiler")
library("Rcpp")
#pSI is a function to describe the probability of being in a state S,I given R and S0 and I0
rm(pSI)
options(expressions=10000)
pSI <- function(s,i,R,s0,i0){
    n <- s0+i0
    if(s > s0 | i > n){return(0)}#limit impossible states
    if(s == s0 & i==i0){return(1)}#probability of being in the start state = 1
    if(i <= 1){return(pSI(s,i+1,R,s0,i0) *  n/(R*s + n))}  #State i = 0,1 can only be reached by recovery, because we start with at least one infected animal
    return(pSI(s + 1,i - 1,R,s0,i0) * (R * (s + 1)/(R * (s + 1) + n)) +   pSI(s,i+1,R,s0,i0) * n/(R*s +n)) #otherwise this state was reached from either s+1, i-1 by infection or from i +1 by recovery
}

pSI_comp <- cmpfun(pSI)

####Function to determine the probability of an observation x1, ..., xm final sizes given R#####
pFS<- function(R,x,s0,i0){
  #Itterate over the trials 
  pfs <- 1
  for(m in c(1:length(x))){
    pfs<-pfs*pSI_comp(s0[m] - x[m], 0, R, s0[m],i0[m])
  }
  return(pfs)
}

 pFS_comp <- cmpfun(pFS) #function(R,x,s0,i0){
#   #Itterate over the trials 
#   pfs <- 1
#   for(m in c(1:length(x))){
#     pfs<-pfs*pSI_comp(s0[m] - x[m], 0, R, s0[m],i0[m])
#   }
#  return(pfs)
#})

system.time(sapply(rlnorm(10,1,1),pFS_comp, x =c(1,2,3,4), s0 = c(50,10,10,10),i0=c(1,1,1,1)))
system.time(sapply(rlnorm(10,1,1),pFS, x =c(1,2,3,4), s0 = c(50,10,10,10),i0=c(1,1,1,1)))
####Function to determine the probability of an observation x1, ..., xm final sizes given R#####
### this function will first create a generation based table with states

FSdist.tab <- function(R,s0,i0){
  n = s0+i0;
  
}
# s0 = 7
# i0 = 7
# n = s0+i0
# R =3.5
# prec =function(s){1-s*R/(s*R + n)}
# gen.mat <- matrix(rep(0, (s0+i0+1)*(s0+1)),nrow= (s0+1))
# gen.mat[s0+1,i0+1]<- 1;
# #fill first
# for(j in c((s0+1):1)){gen.mat[s0+1,j-1]<- gen.mat[s0+1,j]*prec(s0) }
# gen.mat
# #next rows added by infection
# for(i in c(0:(s0-1))){
#     for(j in c((s0+i0+1):1)){
#           gen.mat[s0-i,j]<- gen.mat[s0-i,j]*prec(s0-i-1) + gen.mat[s0-(i-1),j-1]*(1-prec(s0-(i-1)-1))
#       }
#   }
# gen.mat

pFS<- function(R,x,s0,i0){
  #Itterate over the trials 
  pfs <- 1
  for(m in c(1:length(x))){
    pfs<-pfs*pSI(s0[m] - x[m], 0, R, s0[m],i0[m])
  }
  return(pfs)
}

pFS_comp <- cmpfun( function(R,x,s0,i0){
  #Itterate over the trials 
  pfs <- 1
  for(m in c(1:length(x))){
    pfs<-pfs*pSI_comp(s0[m] - x[m], 0, R, s0[m],i0[m])
  }
  return(pfs)
})

#return the complete final size distribution 
FSdist <- function(R,s0,i0){
  max.cases <- max(s0+i0)
  dist <- NULL
  for(m in c(1:length(s0))){
    dist<-rbind(dist,sapply(c(0:max.cases), FUN = pSI, i = 0, R = R, s0 = s0[m],i0 = i0[m]))
  }
  return(dist)
}

system.time(FSdist(1,c(5,5),c(1,1)))


FSdist.tab <- function(R,s0,i0){
  n = s0+i0;
  
}
# s0 = 7
# i0 = 7
# n = s0+i0
# R =3.5
# prec =function(s){1-s*R/(s*R + n)}
# gen.mat <- matrix(rep(0, (s0+i0+1)*(s0+1)),nrow= (s0+1))
# gen.mat[s0+1,i0+1]<- 1;
# #fill first
# for(j in c((s0+1):1)){gen.mat[s0+1,j-1]<- gen.mat[s0+1,j]*prec(s0) }
# gen.mat
# #next rows added by infection
# for(i in c(0:(s0-1))){
#     for(j in c((s0+i0+1):1)){
#           gen.mat[s0-i,j]<- gen.mat[s0-i,j]*prec(s0-i-1) + gen.mat[s0-(i-1),j-1]*(1-prec(s0-(i-1)-1))
#       }
#   }
# gen.mat

pFS<- function(R,x,s0,i0){
  #Itterate over the trials 
  pfs <- 1
  for(m in c(1:length(x))){
    pfs<-pfs*pSI(s0[m] - x[m], 0, R, s0[m],i0[m])
  }
  return(pfs)
}

pFS_comp <- cmpfun( function(R,x,s0,i0){
  #Itterate over the trials 
  pfs <- 1
  for(m in c(1:length(x))){
    pfs<-pfs*pSI_comp(s0[m] - x[m], 0, R, s0[m],i0[m])
  }
  return(pfs)
})

#return the complete final size distribution 
FSdist <- function(R,s0,i0){
  max.cases <- max(s0+i0)
  dist <- NULL
  for(m in c(1:length(s0))){
    dist<-rbind(dist,sapply(c(0:max.cases), FUN = pSI, i = 0, R = R, s0 = s0[m],i0 = i0[m]))
  }
  return(dist)
}

system.time(FSdist(1,c(5,5),c(1,1)))




#hand calculation compared to this one given R0 = 2 and x = 1
#2 R0/(2 R0 + N) * (1 - 1 R0/(1 R0 + N))^2
2 * 2/(2 * 2 + 3) * (1 -  2 / (2 + 3))^2

pFS(1,0,2,1)
pFS(2,1,2,1)
pFS(2,2,2,1)
pFS(2,0,2,1) + pFS(2,1,2,1) +pFS(2,2,2,1)






###Function to determine the probability of more extreme values given R####
# r = R, x = final number of cases, s0 = initial susceptibles, i0 = initial infectious, comp = the type of extreme 
rm(pExtremes)
pExtremes<-  function(r,x,s0,i0,comp = `<`){
  #create all possible outcomes of these transmission experiments
  #this means all possibilities between 0 and s0 contact infection (hence s0 + 1 options per trial)
  out <- matrix(ncol = length(s0),nrow = prod(s0+1))
  #repeat the outcome as many times as the previous trial possibilities
  for(k in c(1:length(s0)))
  {
    #check how often to repeat the same number given previous trials
    repetition <- ifelse(k > 1,prod(s0[1:k-1]+1),1)
    #put it in the matrix
    out[,k]<- matrix(sapply(c(0:s0[k]),FUN = function(x){rep(x,repetition)}), ncol = 1, nrow =prod(s0+1))[,1]
  }
  #select the extremes by selecting those for which the total number of cases
  #and the sum of x are given by the comparison "comp"  thus either <,>,<= or >= 
  #calculate for each extreme the probability of the Final Size under the hypothesis R = r
  #and sum all probabilities of these extreme outcomes
  #If only one extreme outcome possible: 
  if(is.null(dim(out[comp(apply(out,1,sum),sum(x)),]))){return(pFS(r,out[comp(apply(out,1,sum),sum(x)),],s0,i0))}
  #else
  return(c(sum(apply(out[comp(apply(out,1,sum),sum(x)),],1,function(ext){pFS(r,ext,s0,i0)}))))
}

pExtremes_comp<- cmpfun(function(r,x,s0,i0,comp = `<`){
  #create all possible outcomes of these transmission experiments
  #this means all possibilities between 0 and s0 contact infection (hence s0 + 1 options per trial)
  out <- matrix(ncol = length(s0),nrow = prod(s0+1))
  #repeat the outcome as many times as the previous trial possibilities
  for(k in c(1:length(s0)))
  {
    #check how often to repeat the same number given previous trials
    repetition <- ifelse(k > 1,prod(s0[1:k-1]+1),1)
    #put it in the matrix
    out[,k]<- matrix(sapply(c(0:s0[k]),FUN = function(x){rep(x,repetition)}), ncol = 1, nrow =prod(s0+1))[,1]
  }
  #select the extremes by selecting those for which the total number of cases
  #and the sum of x are given by the comparison "comp"  thus either <,>,<= or >= 
  #calculate for each extreme the probability of the Final Size under the hypothesis R = r
  #and sum all probabilities of these extreme outcomes
  #If only one extreme outcome possible: 
  if(is.null(dim(out[comp(apply(out,1,sum),sum(x)),]))){return(pFS_comp(r,out[comp(apply(out,1,sum),sum(x)),],s0,i0))}
  #else
  return(c(sum(apply(out[comp(apply(out,1,sum),sum(x)),],1,function(ext){pFS_comp(r,ext,s0,i0)}))))
})
##################################################################################################
#                                                                                                #
#                  Function to estimate R using the final size method                            #
#                  Includes confidence interval size 1- alpha and R >= 1 test                    #
#                                                                                                #
##################################################################################################
FinalSize<- function(x,s0,i0, alpha = 0.05, onesided = FALSE){
  res <- data.frame(point.est = as.numeric(1), ci.ll = as.numeric(1),ci.ul= as.numeric(1), pval = as.numeric(1))
  #determine the point estimate by optimization of the log-likelihood function
  res$point.est <- optimize(interval = c(0.0,25), f = function(R){-log(pFS(R,x,s0,i0))})$minimum
  #determine the confidence intervals
  #if one-sided is FALSE both sides, either only lower or upper limit of CI
  #lowerlimit is found for values of R for which the probability of extremes above the observations 
  res$ci.ll <- optimize(interval = c(0.0,25),f = function(R){( pExtremes(R,x,s0,i0,comp = `>=`) - alpha / (2 - onesided))^2})$minimum
  #upperlimit is found for values of R for which the probability of extremes below the observations
  res$ci.ul <- optimize(interval = c(0.0,25),f = function(R){( pExtremes(R,x,s0,i0,comp = `<=`) - alpha / (2 - onesided))^2})$minimum
  
  #probability of R >= 1 is found be calculating the probability to find an equal or less positive under the assumption R0 = 1
  res$pval = pExtremes(1,x,s0,i0,comp = `<=`)
  
  return(res)
}

FinalSize_comp <- cmpfun(function(x,s0,i0, alpha = 0.05, onesided = FALSE){
  res <- data.frame(point.est = as.numeric(1), ci.ll = as.numeric(1),ci.ul= as.numeric(1), pval = as.numeric(1))
  #determine the point estimate by optimization of the log-likelihood function
  res$point.est <- optimize(interval = c(0.0,25), f = function(R){-log(pFS_comp(R,x,s0,i0))})$minimum
  #determine the confidence intervals
  #if one-sided is FALSE both sides, either only lower or upper limit of CI
  #lowerlimit is found for values of R for which the probability of extremes above the observations 
  res$ci.ll <- optimize(interval = c(0.0,25),f = function(R){( pExtremes_comp(R,x,s0,i0,comp = `>=`) - alpha / (2 - onesided))^2})$minimum
  #upperlimit is found for values of R for which the probability of extremes below the observations
  res$ci.ul <- optimize(interval = c(0.0,25),f = function(R){( pExtremes_comp(R,x,s0,i0,comp = `<=`) - alpha / (2 - onesided))^2})$minimum
  
  #probability of R >= 1 is found be calculating the probability to find an equal or less positive under the assumption R0 = 1
  res$pval = pExtremes(1,x,s0,i0,comp = `<=`)
  
  return(res)
})
#check with paper of Kroese and De Jong, 2001
s0.KdJ <- c(5,5)
i0.KdJ<- c(5,4)
x.KdJ <- c(2,1)
start <- proc.time()
FinalSize(x.KdJ,s0.KdJ,i0.KdJ)
proc.time()-start

s0.KdJ <- c(5,5)
i0.KdJ<- c(5,5)
x.KdJ <- c(1,0)
start <- proc.time()
FinalSize(x.KdJ,s0.KdJ,i0.KdJ)
proc.time()-start

###test compiled functions
#uncompiled
s0.KdJ <- c(5,5)
i0.KdJ<- c(5,4)
x.KdJ <- c(2,1)
res.uncomp <-c()
for(i in c(1,10))
{ 
  start <- proc.time()
  for(j in 1:i){FinalSize(x.KdJ,s0.KdJ,i0.KdJ)}
  res.uncomp<- rbind(res.uncomp,(proc.time()-start)[3])
}

#uncompiled
s0.KdJ <- c(5,5)
i0.KdJ<- c(5,4)
x.KdJ <- c(2,1)
res.comp <-c()
for(i in c(1,10))
{ 
  start <- proc.time()
  for(j in 1:i){FinalSize_comp(x.KdJ,s0.KdJ,i0.KdJ)}
  res.comp<- rbind(res.comp,(proc.time()-start)[3])
}
res.comp
res.uncomp

#######################################################################################
#                                                                                     #
#                              Test if R of two treatments is the same                #
#                                                                                     #
#######################################################################################
rm(Test.TwoPops)
# parameter treat defines which colums are those of the treatment group
#this is a vector with booleans of the length of the number of experiments in x, s0, i0
Test.TwoPops<-  function(x,s0,i0,treat){
  #create all possible outcomes of these transmission experiments
  #this means all possibilities between 0 and s0 contact infection (hence s0 + 1 options per trial)
  out <- matrix(ncol = length(s0),nrow = prod(s0+1))
  #repeat the outcome as many times as the previous trial possibilities
  for(k in c(1:length(s0)))
  {
    #check how often to repeat the same number given previous trials
    repetition <- ifelse(k > 1,prod(s0[1:k-1]+1),1)
    #put it in the matrix
    out[,k]<- matrix(sapply(c(0:s0[k]),FUN = function(x){rep(x,repetition)}), ncol = 1, nrow =prod(s0+1))[,1]
  }
  
  #Select those possible outcomes for which the difference is <= sum(x[!treat]) - sum(x[treat])
  #and sum all probabilities of these same or extreme outcomes
  likeFun <- function(R){sum(apply(out[sapply(X = apply(out[,!treat],1,sum)-apply(out[,treat],1,sum), FUN = function(z){abs(z) >=  sum(x[!treat])-sum(x[treat])}),],1,function(ext){pFS(R,ext,s0,i0)}))}
  #find R with the maximum probability of occurence, this is the p-value for R1 = R2
  return(optimize(interval = c(0,25),f = likeFun, maximum = TRUE)$objective)
}

#check with paper of Kroese and De Jong, 2001
s0.KdJ.ct <- c(5,5,5,5)
i0.KdJ.ct <- c(5,5,5,5)
x.KdJ.ct <- c(2,1,5,5)
treat.KdJ.ct <- c(TRUE,TRUE,FALSE,FALSE)

Test.TwoPops(x = x.KdJ.ct,s0  = s0.KdJ.ct,i0 = i0.KdJ.ct, treat = treat.KdJ.ct)



