library(MASS)
library(mvtnorm)

# function for z test
z.test = function(a, mu=0, var=1){
  zeta = (mean(a) - mu) / (sqrt(var / length(a)))
  return(zeta)
}

########################################################################
empricalFDR_m <- function(rho,m,B,stat){
  
  #n.reject
  R0_t <- c()
  R_t <- c()
  fdr <- c()
  
  # create correlation mareix with equal correlation
  sigma <- diag(m)+matrix(rep(rho,m*m),nrow=m)-diag(rep(rho,m))
  
  # Simulation start
  for (i in 1:B) {
    # create multivariate normal observations
    df <- data.frame(mvrnorm(n = 200, mu=c(rep(0,0.9*m),rep(0.5,0.1*m)), sigma))
    
    # calculate z stats
    mul.z <- apply(df,2,z.test)
    
    # numberof false rejection
    #n.reject <- sum(mul.z>=qnorm(0.95))
    R0_t[i] <- sum(mul.z[1:(0.9*m)]>=qnorm(0.95))
    R_t[i] <- sum(mul.z>=qnorm(0.95))
    R_t[i] <- max(R_t[i],1)
    
    fdr[i] <- R0_t[i]/R_t[i]
    
  }
  
  # calculate empirical mean/var of # rejections
  mean.fdr <- mean(fdr)
  var.fdr <- var(fdr)
  
  if (stat=='mean'){
    return(mean.fdr)
  } else if (stat=='var'){
    return(var.fdr)
  }
 
  

}

#  mean


rho.list <- seq(0,1,by=0.05)
set.seed(100)
error.mean_m  <- sapply(rho.list, empricalFDR_m, m=100, B=1000,stat='mean')
plot(error.mean_m ~rho.list,xlab='Correlation', ylab='Empirical FDR' )





##########################################################################################
 #                                     E(v)/E(R)
##########################################################################################

empricalEvEr_m <- function(rho,m,B,stat){
  
  #n.reject
  R0_t <- c()
  R_t <- c()
  fdr <- c()
  
  # create correlation mareix with equal correlation
  sigma <- diag(m)+matrix(rep(rho,m*m),nrow=m)-diag(rep(rho,m))
  
  # Simulation start
  for (i in 1:B) {
    # create multivariate normal observations
    df <- data.frame(mvrnorm(n = 200, mu=c(rep(0,0.9*m),rep(0.1,0.1*m)), sigma))
    
    # calculate z stats
    mul.z <- apply(df,2,z.test)
    
    # numberof false rejection
    #n.reject <- sum(mul.z>=qnorm(0.95))
    R0_t[i] <- sum(mul.z[1:(0.9*m)]>=qnorm(0.95))
    R_t[i] <- sum(mul.z>=qnorm(0.95))
    R_t[i] <- max(R_t[i],1)
    
    fdr[i] <- R0_t[i]/R_t[i]
    
  }
  
  # calculate E(v)/E(R)
  mean.R0_t<- mean(R0_t)
  mean.R_t<- mean(R_t)
  eVeR <- mean.R0_t/mean.R_t
  
  return(eVeR)
#   var.fdr <- var(fdr)
#   
#   if (stat=='mean'){
#     return(mean.fdr)
#   } else if (stat=='var'){
#     return(var.fdr)
#   }
  
  
  
}

#  mean
rho.list <- seq(0,1,by=0.05)
#set.seed(100)
error.mean_m  <- sapply(rho.list, empricalEvEr_m, m=100, B=1000)
plot(error.mean_m ~rho.list,xlab='Correlation', ylab='Empirical E(v)/E(R) (mean)' )
