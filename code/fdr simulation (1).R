
# simulation: Empirical FDR
library(MASS)

# function for z test
z.test = function(a, mu=0, var=1){
  zeta = (mean(a) - mu) / (sqrt(var / length(a)))
  return(zeta)
}

# function for empricalFDR
empricalFDR <- function(rho,B){
  
# initilize z statieitcs
z1 <- c()
z2 <- c()

# create covariance matrix 
sigma <- matrix(c(1,rho ,rho ,1),2,2)

# Simulation
for (i in 1:B) {
  
# create dataset: x1 is true null; x2 is false null
df <- data.frame(mvrnorm(n = 200, mu=c(0,0.2), sigma))

# calculate z stats
z1[i] <- z.test(df$X1)
z2[i] <- z.test(df$X2)

}

# one-sided test 
reject_z1 <- ifelse(z1>=qnorm(0.95),1,0)
reject_z2 <- ifelse(z2>=qnorm(0.95),1,0)

# calculate empirical fdr
#fdr.empirical = sum(reject_z1==1)/(sum(reject_z1==1)+sum(reject_z2==1))

# calculate # rejections
R = sum(reject_z1==1) + sum(reject_z2==1)
V = sum(reject_z1==1)
# return (list(R,V))
return(V/R/B)
}


rho.list <- seq(0,1,by=0.1)
error.er <- sapply(rho.list, empricalFDR,B=5000)
plot(error.er~rho.list )

# empricalFDR.all <- replicate(500,empricalFDR(0.9,1000))
# empricalR.all <- empricalFDR.all[seq(1,by=2,len=500)]
# empricalV.all <- empricalFDR.all[seq(2,by=2,len=500)]
# mean(unlist(empricalR.all))
# sd(unlist(empricalR.all))
# mean(unlist(empricalV.all))
# sd(unlist(empricalV.all))
# 
# empricalFDR2.all <- replicate(500,empricalFDR(0.9,1000))
# empricalR2.all <- empricalFDR2.all[seq(1,by=2,len=500)]
# empricalV2.all <- empricalFDR2.all[seq(2,by=2,len=500)]
# mean(unlist(empricalR2.all))
# sd(unlist(empricalR2.all))
# mean(unlist(empricalV2.all))
# sd(unlist(empricalV2.all))

#######################################

# calculate exact v(t)
library(mvtnorm)


# function to calculate the expected number of false rejections based on rho
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




