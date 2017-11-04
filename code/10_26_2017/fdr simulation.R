
# simulate 
library(MASS)

rho <- 1

# function for z test
z.test = function(a, mu=0, var=1){
  zeta = (mean(a) - mu) / (sqrt(var / length(a)))
  return(zeta)
}

empricalFDR <- function(rho){
sigma <- matrix(c(1,rho ,rho ,1),2,2)

# function for z test
z.test = function(a, mu=0, var=1){
  zeta = (mean(a) - mu) / (sqrt(var / length(a)))
  return(zeta)
}

# initilize z statieitcs
z1 <- c()
z2 <- c()

# Simulation
B=1000
for (i in 1:B) {
  
# create dataset: x1 is true null; x2 is false null
df <- data.frame(mvrnorm(n = 1000, mu=c(0,0.5), sigma))

# calculate z stats
z1[i] <- z.test(df$X1)
z2[i] <- z.test(df$X2)

}

# rejection
reject_z1 <- ifelse(z1>=qnorm(0.95),1,0)
reject_z2 <- ifelse(z2>=qnorm(0.95),1,0)

# calculate empirical fdr
#fdr.empirical = sum(reject_z1==1)/(sum(reject_z1==1)+sum(reject_z2==1))

# calculate # false rejections
return (sum(reject_z1==1))

}

rho.list <- seq(-1,1,by=0.1)
error.er <- sapply(rho.list, empricalFDR)
plot(error.er~rho.list )


#######################################

# calculate exact v(t)
library(mvtnorm)
mean <- c(0,0)  # under the null
lower <- c(-Inf, -Inf)
upper <- c(-qnorm(0.95), -qnorm(0.95))
rho <- 1
sigma <- matrix(c(1,rho ,rho ,1),2,2)
pmvnorm(lower, upper, mean, sigma)

mean <- c(0,0)  # under the null
lower <- c(1.96, 1.96)
upper <- c(Inf, Inf)
rho <- 1
sigma <- matrix(c(1,rho ,rho ,1),2,2)
pmvnorm(lower, upper, mean, sigma)

mean <- c(0,0)  # under the null
lower <- c(1.96, 1.96)
upper <- c(Inf, Inf)
rho <- seq()
sigma <- matrix(c(1,rho ,rho ,1),2,2)
pmvnorm(lower, upper, mean, sigma)



rho.list <- seq(-1,1,by=0.1)

# function to calculate the expected number of false rejections based on rho
R0_T <- function(rho) {
  mean <- c(0,0)
  sigma <- matrix(c(1,rho ,rho ,1),2,2)
  prob1 <- pmvnorm(lower=c(qnorm(0.95), -Inf), upper=c(Inf,qnorm(0.95)), mean, sigma) * 2 # prob of falsely reject either one 
  prob2 <- pmvnorm(lower=c(qnorm(0.95), qnorm(0.95)), upper=c(Inf,Inf), mean, sigma) * 2 # prob of falsely reject both two
  E.R0_t <- prob1 + 2*prob2 
  return(E.R0_t)
}
rho.list <- seq(-1,1,by=0.1)
error <- sapply(rho.list, R0_T)
plot(error~rho.list )

rho=0
mean <- c(0,0)
sigma <- matrix(c(1,rho ,rho ,1),2,2)
prob1 <- pmvnorm(lower=c(qnorm(0.95), -Inf), upper=c(Inf,qnorm(0.95)), mean, sigma) * 2 # prob of falsely reject either one 
prob2 <- pmvnorm(lower=c(qnorm(0.95), qnorm(0.95)), upper=c(Inf,Inf), mean, sigma) * 2 # prob of falsely reject both two
E.R0_t <- prob1 + 2*prob2 

