library(MASS)
library(mvtnorm)

# function for z test
z.test = function(a, mu=0, var=1){
  zeta = (mean(a) - mu) / (sqrt(var / length(a)))
  return(zeta)
}





#############################  Empirical FDR m=3 ############################################
# m=3
# function for empricalFDR
empricalFDR_m3 <- function(rho,B,stat){
  
  # initilize storage vectors
  z1 <- c()
  z2 <- c()
  z3 <- c()
  R0_t <- c()
  
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
   #z3[i] <- z.test(df$X3)
    
    # numberof false rejection
    R0_t[i] <- ifelse(z1[i]>=qnorm(0.95),1,0) + ifelse(z2[i]>=qnorm(0.95),1,0)
  }
  
  
  
  # calculate empirical mean/var of # rejections
  mean.R0_t <- mean(R0_t)
  var.R0_t <- var(R0_t)
  
  if (stat=='mean'){
    return(mean.R0_t)
  } else if (stat=='var'){
    return(var.R0_t)
  }
}

##############  mean #############
rho.list <- seq(0,1,by=0.1)
error.mean_m3  <- sapply(rho.list, empricalFDR_m3,B=5000,stat='mean')
plot(error.mean_m3 ~rho.list,xlab='Correlation', ylab='Empirical R0' )
abline(h=0.10,col='red',lty=2)

##############  variance ##############
rho.list <- seq(0,1,by=0.1)
error.var_m3 <- sapply(rho.list, empricalFDR_m3,B=5000,stat='var')
plot(error.var_m3~rho.list,xlab='Correlation', ylab='Empirical R0' )

#  the number of average false rejections seems have no association with rho and mu3?

#############################  Exact FDR m=3 ############################################

######## function to calculate the expected number of false rejections based on rho########
R0_T <- function(rho) {
  
  mean <- c(0,0,0)
  sigma <- matrix(c(1,rho ,rho, 
                    rho,1 ,rho,
                    rho ,rho,1 ),3,3)
  
  prob1 <- pmvnorm(lower=c(qnorm(0.95), -Inf, -Inf), upper=c(Inf,qnorm(0.95),qnorm(0.95)), mean, sigma) *3 # prob of falsely reject one
  prob2 <- pmvnorm(lower=c(qnorm(0.95), qnorm(0.95), -Inf), upper=c(Inf,Inf,qnorm(0.95)), mean, sigma)*3  # prob of falsely reject two
  prob3 <- pmvnorm(lower=c(qnorm(0.95), qnorm(0.95), qnorm(0.95)), upper=c(Inf,Inf,Inf), mean, sigma)  # prob of falsely reject three
  E.R0_t <- prob1 + 2*prob2 + 3*prob3
  return(E.R0_t)
}

rho.list <- seq(0,1,by=0.1)
error <- sapply(rho.list, R0_T)*2/3
plot(error~rho.list ,xlab='Correlation', ylab='E(R0) *2/3')


############  calculate the variacne of R0t  ########
R0_T.var <- function(rho) {
  
  mean <- c(0,0,0)
  sigma <- matrix(c(1,rho ,rho, 
                    rho,1 ,rho,
                    rho ,rho,1 ),3,3)
  
  prob1 <- pmvnorm(lower=c(qnorm(0.95), -Inf, -Inf), upper=c(Inf,qnorm(0.95),qnorm(0.95)), mean, sigma) *3 # prob of falsely reject one
  prob2 <- pmvnorm(lower=c(qnorm(0.95), qnorm(0.95), -Inf), upper=c(Inf,Inf,qnorm(0.95)), mean, sigma)*3  # prob of falsely reject two
  prob3 <- pmvnorm(lower=c(qnorm(0.95), qnorm(0.95), qnorm(0.95)), upper=c(Inf,Inf,Inf), mean, sigma)  # prob of falsely reject three
 
   E.R0_t <- prob1 + 2*prob2 + 3*prob3
   E.R0_t.square <- prob1 + 2^2*prob2 + 3^2*prob3
   var <- E.R0_t.square - E.R0_t^2
  return(var)
}

rho.list <- seq(0,1,by=0.1)
error.var <- sapply(rho.list, R0_T.var)*4/9
plot(error.var~rho.list,xlab='Correlation', ylab='Var(R0) *4/9' )



#############################  Empirical FDR m>3 ############################################
rho= 0.9
m= 100

empricalFDR_m <- function(rho,m,B,stat){
  
  #n.reject
  n.f_reject <- c()
  
  # create correlation mareix with equal correlation
  sigma <- diag(m)+matrix(rep(rho,m*m),nrow=m)-diag(rep(rho,m))
  
  # Simulation start
  for (i in 1:B) {
    # create multivariate normal observations
    df <- data.frame(mvrnorm(n = 200, mu=c(rep(0.5,0.1*m),rep(0,0.9*m)), sigma))
    
    # calculate z stats
    mul.z <- apply(df,2,z.test)
    
    # numberof false rejection
    #n.reject <- sum(mul.z>=qnorm(0.95))
    n.f_reject[i] <- sum(mul.z[(0.1*m+1):m]>=qnorm(0.95))
  }
  
  # calculate empirical mean/var of # rejections
  mean.R0_t <- mean(n.f_reject)
  var.R0_t <- var(n.f_reject)
  
  if (stat=='mean'){
    return(mean.R0_t)
  } else if (stat=='var'){
    return(var.R0_t)
  }
  
}

#  mean
rho.list <- seq(0,1,by=0.1)
error.mean_m  <- sapply(rho.list, empricalFDR_m, m=100, B=5000,stat='mean')
plot(error.mean_m ~rho.list,xlab='Correlation', ylab='Empirical R0 (mean)' )
abline(h=0.10,col='red',lty=2)

# var
rho.list <- seq(0,1,by=0.1)
error.var_m <- sapply(rho.list, empricalFDR_m,m=100,B=5000,stat='var')
plot(error.var_m ~rho.list,xlab='Correlation', ylab='Empirical R0 (var)' )

