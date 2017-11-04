

# m >= 2
library(MASS)

# create correlation mareix with equal correlation
rho= 0.9
m= 100
sigma <- diag(m)+matrix(rep(rho,m*m),nrow=m)-diag(rep(rho,m))

# create multivariate normal observations
df <- data.frame(mvrnorm(n = 200, mu=c(rep(0.5,0.1*m),rep(0,0.9*m)), sigma))

# z statistics
mul.z <- apply(df,2,z.test)

#
n.reject <- sum(mul.z>=qnorm(0.95))
n.f_reject <- sum(mul.z[(0.1*m+1):m]>=qnorm(0.95))

n.reject
n.f_reject


###############################
# m=3
# function for empricalFDR
empricalFDR <- function(rho,B){
  
  # initilize z statieitcs
  z1 <- c()
  z2 <- c()
  z3 <- c()
  
  # create covariance matrix 
  mean <- c(0,0,5)
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
  }
  
  # one-sided test 
  reject_z1 <- ifelse(z1>=qnorm(0.95),1,0)
  reject_z2 <- ifelse(z2>=qnorm(0.95),1,0)
  reject_z3 <- ifelse(z3>=qnorm(0.95),1,0)
  
  # calculate empirical fdr
  #fdr.empirical = sum(reject_z1==1)/(sum(reject_z1==1)+sum(reject_z2==1))
  
  # calculate # rejections
  R = sum(reject_z1==1) + sum(reject_z2==1)+ sum(reject_z3==1)
  V = sum(reject_z1==1) + sum(reject_z2==1)
 # return (list(R,V))
  return(V/B)
}

empricalFDR(0,1000)
rho.list <- seq(0,1,by=0.1)
error.er1 <- sapply(rho.list, empricalFDR,B=5000)
plot(error.er1~rho.list )
#  the number of average false rejections seems have no association with rho and mu3?

###################################################
library(mvtnorm)


# function to calculate the expected number of false rejections based on rho
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
error <- sapply(rho.list, R0_T)
plot(error~rho.list )


