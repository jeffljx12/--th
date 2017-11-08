####################################################################################
#  
#               explore E(V/R) vs E(v)/E(R)
#
####################################################################################

library(MASS)
library(mvtnorm)


######################################## m=3 #####################################################

# function for z test
z.test = function(a, mu=0, var=1){
  zeta = (mean(a) - mu) / (sqrt(var / length(a)))
  return(zeta)
}



#############################  E(V/R) vs rho ############################################
# m=3
# function for empricalFDR
eVR_m3 <- function(rho,B,stat){
  
  # initilize storage vectors
  z1 <- c()
  z2 <- c()
  z3 <- c()
  R0_t <- c()
  R_t <- c()
  empricalFDR <- c()
  
  # create covariance matrix 
  mean <- c(0,0,0.5)
  sigma <- matrix(c(1,rho ,rho, 
                    rho,1 ,rho,
                    rho ,rho,1 ),3,3)
  
  # Simulation
  for (i in 1:B) {
    
    # create dataset: x1 is true null; x2 is false null
    df <- data.frame(mvrnorm(n = 200, mu=mean, sigma))
    
    # calculate z stats
    z1[i] <- z.test(df$X1)
    z2[i] <- z.test(df$X2)
    z3[i] <- z.test(df$X3)
    
    # numberof false rejection
    R0_t[i] <- ifelse(z1[i]>=qnorm(0.95),1,0) + ifelse(z2[i]>=qnorm(0.95),1,0)
    
    # numberof all rejection
    R_t[i] <- ifelse(z1[i]>=qnorm(0.95),1,0) + ifelse(z2[i]>=qnorm(0.95),1,0)+ifelse(z3[i]>=qnorm(0.95),1,0)
    R_t[i] <- max(R_t[i],1)
    
    # emprical fdr
    empricalFDR[i] <- R0_t[i]/R_t[i]
  }
  
  # numberof all rejection
  
  # calculate empirical mean/var of fdr
  mean.empricalFDR <- mean(empricalFDR)
  var.empricalFDR <- var(empricalFDR)
  
  if (stat=='mean'){
    return(mean.empricalFDR)
  } else if (stat=='var'){
    return(var.empricalFDR)
  }
}

##############  eVR vs rho #############
rho.list <- seq(0,1,by=0.05)
error.mean_m3  <- sapply(rho.list, eVR_m3 ,B=5000,stat='mean')
plot(error.mean_m3 ~rho.list,xlab='Correlation', ylab='Empirical FDR(mean)' )


##############  fdr.var vs rho ##############
rho.list <- seq(0,1,by=0.05)
error.var_m3 <- sapply(rho.list, eVR_m3 ,B=5000,stat='var')
plot(error.var_m3~rho.list,xlab='Correlation', ylab='Empirical FDR(var)' )











###############################################################################################################
#############################  E(V)/E(R) vs rho  vs rho ############################################
# m=3
# function for empricalFDR
eVeR_m3 <- function(rho,B,stat){
  
  # initilize storage vectors
  z1 <- c()
  z2 <- c()
  z3 <- c()
  R0_t <- c()
  R_t <- c()
  #empricalFDR <- c()
  
  # create covariance matrix 
  mean <- c(0,0,0.1)
  sigma <- matrix(c(1,rho ,rho, 
                    rho,1 ,rho,
                    rho ,rho,1 ),3,3)
  
  # Simulation
  for (i in 1:B) {
    
    # create dataset: x1 is true null; x2 is false null
    df <- data.frame(mvrnorm(n = 200, mu=mean, sigma))
    
    # calculate z stats
    z1[i] <- z.test(df$X1)
    z2[i] <- z.test(df$X2)
    z3[i] <- z.test(df$X3)
    
    # numberof false rejection
    R0_t[i] <- ifelse(z1[i]>=qnorm(0.95),1,0) + ifelse(z2[i]>=qnorm(0.95),1,0)
    
    # numberof all rejection
    R_t[i] <- ifelse(z1[i]>=qnorm(0.95),1,0) + ifelse(z2[i]>=qnorm(0.95),1,0)+ifelse(z3[i]>=qnorm(0.95),1,0)
    R_t[i] <- max(R_t[i],1)
    # emprical fdr
    #empricalFDR[i] <- R0_t[i]/max(R_t[i],1)
  }
  
  # numberof all rejection
  
  # calculate empirical mean/var of fdr
  mean.eVeR <- mean(R0_t)/mean(R_t)
  var.eVeR <- var(R0_t)/var(R_t)
  
  if (stat=='mean'){
    return(mean.eVeR)
  } else if (stat=='var'){
    return(var.eVeR )
  }
}

##############  eVeR vs rho #############
rho.list <- seq(0,1,by=0.05)
error.mean_m3  <- sapply(rho.list, eVeR_m3 ,B=5000,stat='mean')
plot(error.mean_m3 ~rho.list,xlab='Correlation', ylab='E(v)/E(R) (mean)' )


##############  fdr.var vs rho ##############
# rho.list <- seq(0,1,by=0.05)
# error.var_m3 <- sapply(rho.list, eVeR_m3 ,B=5000,stat='var')
# plot(error.var_m3~rho.list,xlab='Correlation', ylab='var(v)/var(R) ' )















##################################################################################################################
####################################################################################
#  
#               explore E(V/R) vs E(v)/E(R), m=100
#
####################################################################################

#############################  Empirical FDR m>3 ############################################

empricalFDR_m <- function(rho,m,B,stat){
  
  # number of total rejection
  n.f_reject <- c()
  # number of false rejection
  n.reject <- c()
 
  
  # create correlation mareix with equal correlation
  sigma <- diag(m)+matrix(rep(rho,m*m),nrow=m)-diag(rep(rho,m))
  
  # Simulation start
  for (i in 1:B) {
    # create multivariate normal observations
    df <- data.frame(mvrnorm(n = 200, mu=c(rep(0.5,0.1*m),rep(0,0.9*m)), sigma))
    
    # calculate z stats
    mul.z <- apply(df,2,z.test)
    
    # number of rejections
    n.reject <- sum(mul.z >= qnorm(0.95))
    n.reject <-max(n.reject,1)
    
    n.f_reject <- sum(mul.z[(0.1*m+1):m]>=qnorm(0.95))
    
    fdr <- n.f_reject/n.reject
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
error.mean_m  <- sapply(rho.list, empricalFDR_m, m=100, B=5000,stat='mean')
plot(error.mean_m ~rho.list,xlab='Correlation', ylab='Empirical fdr (mean)' )

#  var
rho.list <- seq(0,1,by=0.05)
error.var_m  <- sapply(rho.list, empricalFDR_m, m=100, B=5000,stat='var')
plot(error.var_m ~rho.list,xlab='Correlation', ylab='Empirical fdr (var)' )




###################################### E(v)/E(R) , m=100 ##################################
eVeR_m <- function(rho,m,B,stat){
  
  # number of total rejection
  n.f_reject <- c()
  
  # number of false rejection
  
  # create correlation mareix with equal correlation
  sigma <- diag(m)+matrix(rep(rho,m*m),nrow=m)-diag(rep(rho,m))
  
  # Simulation start
  for (i in 1:B) {
    # create multivariate normal observations
    df <- data.frame(mvrnorm(n = 200, mu=c(rep(0.5,0.1*m),rep(0,0.9*m)), sigma))
    
    # calculate z stats
    mul.z <- apply(df,2,z.test)
    
    # number of rejections
    n.reject <- sum(mul.z>=qnorm(0.95))
    n.reject <-max(n.reject,1)
    
    n.f_reject <- sum(mul.z[(0.1*m+1):m]>=qnorm(0.95))
    
    #fdr <- n.f_reject/n.reject
  }
  
  
  # calculate empirical mean/var of # rejections
  mean.eVeR <- mean(n.f_reject)/mean(n.reject)
  var.eVeR <- var(n.f_reject)/var(n.reject)


  if (stat=='mean'){
    return(mean.eVeR)
  } else if (stat=='var'){
    return(var.eVeR)
  }
  
}

#  mean
rho.list <- seq(0,1,by=0.05)
error.mean_m  <- sapply(rho.list, eVeR_m, m=100, B=5000,stat='mean')
plot(error.mean_m ~rho.list, xlab='Correlation', ylab='E(v)/E(R) (mean)' )




