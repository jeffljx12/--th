

#############################  Empirical FDR  ############################################
library(MASS)

# function for z test
z.test = function(a, mu=0, var=1){
  zeta = (mean(a) - mu) / (sqrt(var / length(a)))
  return(zeta)
}

# function for empricalFDR
empricalFDR <- function(rho,B,stat){
  
   # initilize storage vectors
    z1 <- c()
    z2 <- c()
    R0_t <- c()
    
  # create covariance matrix 
   sigma <- matrix(c(1,rho ,rho ,1),2,2)

  # Simulation start 
    for (i in 1:B) {
  
       # create dataset: x1 is true null; x2 is false null
       df <- data.frame(mvrnorm(n = 200, mu=c(0,0.5), sigma))

       # calculate z stats
       z1[i] <- z.test(df$X1)
       #z2[i] <- z.test(df$X2)
       
       # numberof false rejection
       R0_t[i] <- ifelse(z1[i]>=qnorm(0.95),1,0)
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


##############  mean ##############
rho.list <- seq(0,1,by=0.1)
error.mean <- sapply(rho.list, empricalFDR,B=1000, stat='mean')
plot(error.mean ~rho.list,xlab='Correlation', ylab='Empirical R0' )
abline(h=0.05,col='red',lty=2)

##############  variance ##############
rho.list <- seq(0,1,by=0.1)
error.var <- sapply(rho.list, empricalFDR,B=5000,stat='var')
plot(error.var~rho.list,xlab='Correlation', ylab='Empirical R0' )



#############################  Exact FDR  ############################################

# calculate exact v(t)
library(mvtnorm)

######## function to calculate the expected number of false rejections based on rho########
R0_T <- function(rho) {
  mean <- c(0,0)
  sigma <- matrix(c(1,rho ,rho ,1),2,2)
  prob1 <- pmvnorm(lower=c(qnorm(0.95), -Inf), upper=c(Inf,qnorm(0.95)), mean, sigma)*2  # prob of falsely reject one 
  prob2 <- pmvnorm(lower=c(qnorm(0.95), qnorm(0.95)), upper=c(Inf,Inf), mean, sigma)  # prob of falsely reject both two
  E.R0_t <- prob1 + 2*prob2 
  return(E.R0_t)
}
rho.list <- seq(0,1,by=0.1)
error <- sapply(rho.list, R0_T)*1/2
plot(error~rho.list ,xlab='Correlation', ylab='E(R0) *1/2')

############  calculate the variacne of R0t  ########
R0_T.var <- function(rho) {
  mean <- c(0,0)
  sigma <- matrix(c(1,rho ,rho ,1),2,2)
  prob1 <- pmvnorm(lower=c(qnorm(0.95), -Inf), upper=c(Inf,qnorm(0.95)), mean, sigma)*2  # prob of falsely reject one 
  prob2 <- pmvnorm(lower=c(qnorm(0.95), qnorm(0.95)), upper=c(Inf,Inf), mean, sigma)  # prob of falsely reject both two
  E.R0_t <- prob1 + 2*prob2 
  E.R0_t.square <- prob1 + 2^2*prob2 
  
  var <- E.R0_t.square - E.R0_t^2
  return(var)
}

rho.list <- seq(0,1,by=0.1)
error.var <- sapply(rho.list, R0_T.var)*((1/2)^2)
plot(error.var~rho.list,xlab='Correlation', ylab='Var(R0) *1/4' )



