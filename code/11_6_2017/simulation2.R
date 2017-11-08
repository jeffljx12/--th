library(MASS)
library(mvtnorm)


######################################## m=3 #####################################################


# function for z test
z.test = function(a, mu=0, var=1){
  zeta = (mean(a) - mu) / (sqrt(var / length(a)))
  return(zeta)
}


#############################  Empirical FDR vs rho ############################################
# m=3
# function for empricalFDR
empricalFDR_m3 <- function(rho,B,stat){
  
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
    
    # emprical fdr
    empricalFDR[i] <- R0_t[i]/max(R_t[i],1)
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

##############  fdr.mean vs rho #############
rho.list <- seq(0,1,by=0.05)
error.mean_m3  <- sapply(rho.list, empricalFDR_m3,B=5000,stat='mean')
plot(error.mean_m3 ~rho.list,xlab='Correlation', ylab='Empirical FDR(mean)' )


##############  fdr.var vs rho ##############
rho.list <- seq(0,1,by=0.05)
error.var_m3 <- sapply(rho.list, empricalFDR_m3,B=5000,stat='var')
plot(error.var_m3~rho.list,xlab='Correlation', ylab='Empirical FDR(var)' )


###############################################################################################################

#############################  Empirical FDR vs alpha level ############################################
# m=3
# function for empricalFDR
empricalFDR_m3.2 <- function(rho,B,stat, alpha){
  
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
    R0_t[i] <- ifelse(z1[i]>=qnorm(1-alpha),1,0) + ifelse(z2[i]>=qnorm(1-alpha),1,0)
    
    # numberof all rejection
    R_t[i] <- ifelse(z1[i]>=qnorm(1-alpha),1,0) + ifelse(z2[i]>=qnorm(1-alpha),1,0)+ifelse(z3[i]>=qnorm(1-alpha),1,0)
    
    # emprical fdr
    empricalFDR[i] <- R0_t[i]/max(R_t[i],1)
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

##############  fdr.mean vs alpha #############
alpha.list <- seq(0.05,0.95,by=0.05)
error.mean_m3  <- sapply(alpha.list , empricalFDR_m3.2,B=5000,stat='mean',rho=0.5)
plot(error.mean_m3 ~alpha.list,xlab='alpha level', ylab='Empirical FDR(mean), rho=0.5' )
abline(a=0,b=1,col="red", lty=2)

##############  fdr.var vs alpha #############
alpha.list <- seq(0.05,0.95,by=0.05)
error.var_m3  <- sapply(alpha.list , empricalFDR_m3.2,B=5000,stat='var',rho=0.5)
plot(error.mean_m3 ~alpha.list,xlab='alpha level', ylab='Empirical FDR(var), rho=0.5' )









#############################  Empirical # rejections vs rho ############################################
# m=3
# function for empricalFDR
empricalRt_m3 <- function(rho,B,stat){
  
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
    #R0_t[i] <- ifelse(z1[i]>=qnorm(0.95),1,0) + ifelse(z2[i]>=qnorm(0.95),1,0)
    
    # numberof all rejection
    R_t[i] <- ifelse(z1[i]>=qnorm(0.95),1,0) + ifelse(z2[i]>=qnorm(0.95),1,0)+ifelse(z3[i]>=qnorm(0.95),1,0)
    R_t[i] <- max(R_t[i],1)
    # emprical fdr
   # empricalFDR[i] <- R0_t[i]/max(R_t[i],1)
  }
  
  # numberof all rejection
  
  # calculate empirical mean/var of fdr
  mean.R_t <- mean(R_t)
  var.R_t <- var(R_t)
  
  if (stat=='mean'){
    return(mean.R_t)
  } else if (stat=='var'){
    return(var.R_t)
  }
}

##############  fdr.mean vs rho #############
rho.list <- seq(0,1,by=0.05)
error.mean_m3  <- sapply(rho.list, empricalRt_m3,B=5000,stat='mean')
plot(error.mean_m3 ~rho.list,xlab='Correlation', ylab='Empirical Rt(mean)' )

##############  fdr.var vs rho #############
rho.list <- seq(0,1,by=0.05)
error.var_m3  <- sapply(rho.list, empricalRt_m3,B=5000,stat='var')
plot(error.var_m3 ~rho.list,xlab='Correlation', ylab='Empirical Rt(var)' )



#############################  Empirical # rejections vs alpha ############################################
# m=3
# function for empricalFDR
empricalRt_m3.2 <- function(rho,B,stat,alpha){
  
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
    #R0_t[i] <- ifelse(z1[i]>=qnorm(1-alpha),1,0) + ifelse(z2[i]>=qnorm(1-alpha),1,0)
    
    # numberof all rejection
    R_t[i] <- ifelse(z1[i]>=qnorm(1-alpha),1,0) + ifelse(z2[i]>=qnorm(1-alpha),1,0)+ifelse(z3[i]>=qnorm(1-alpha),1,0)
    R_t[i] <- max(R_t[i],1)
    # emprical fdr
    # empricalFDR[i] <- R0_t[i]/max(R_t[i],1)
  }
  
  # numberof all rejection
  
  # calculate empirical mean/var of fdr
  mean.R_t <- mean(R_t)
  var.R_t <- var(R_t)
  
  if (stat=='mean'){
    return(mean.R_t)
  } else if (stat=='var'){
    return(var.R_t)
  }
}

##############  fdr.mean vs rho #############
alpha.list <- seq(0.05,0.95,by=0.05)
error.mean_m3  <- sapply(rho.list, empricalRt_m3.2,B=5000, rho=0.5, stat='mean')
plot(error.mean_m3 ~rho.list,xlab='Alpha Level', ylab='Empirical Rt(mean)' )

##############  fdr.var vs rho #############
rho.list <- seq(0,1,by=0.05)
error.var_m3  <- sapply(rho.list, empricalRt_m3,B=5000,stat='var')
plot(error.var_m3 ~rho.list,xlab='Alpha Level', ylab='Empirical Rt(var)' )

