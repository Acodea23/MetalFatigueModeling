library(R2jags)
pc <- read.csv(paper_clip_data.csv)
head(pc)

y <- pc$Num_bends
temp <- pc$Temp
gauge <- pc$Gauge
tester <- as.numeric(pc$Block=="M")

PoisModel <- "model {
  for(i in 1:length(y)){
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta[1] + beta[2]*temp[i] + beta[3]*gauge[i] + beta[4]*tester[i]
    
  }
  beta[1] ~ dnorm(0,1/100)
  beta[2] ~ dnorm(0,1/100)
  beta[3] ~ dnorm(0,1/100)
  beta[4] ~ dnorm(0,1/100)
  PPDK ~ dpois(exp(beta[1] + beta[2]*75 + beta[3]*.8 + beta[4]*0))
  PPDM ~ dpois(exp(beta[1] + beta[2]*75 + beta[3]*.8 + beta[4]*1))
  PPDhot ~ dpois(exp(beta[1] + beta[2]*74 + beta[3]*1 + beta[4]*1))
  PPDcold ~ dpois(exp(beta[1] + beta[2]*42 + beta[3]*1 + beta[4]*1))
  PPD.8 ~ dpois(exp(beta[1] + beta[2]*75 + beta[3]*.1 + beta[4]*1))
  PPD1 ~ dpois(exp(beta[1] + beta[2]*75 + beta[3]*1 + beta[4]*1))
}
"
Pois.sim <- jags(
  data=c('y','temp','gauge', 'tester'),
  parameters.to.save=c('beta', 'PPDhot', 'PPDcold', 'PPDK', 'PPDM', 'PPD.8', 'PPD1'),
  model.file=textConnection(PoisModel),
  n.iter=22000,
  n.burnin=2000,
  n.chains=5,
  n.thin=1
) 

beta0 <- Pois.sim$BUGSoutput$sims.matrix[,"beta[1]"]
beta1 <- Pois.sim$BUGSoutput$sims.matrix[,"beta[2]"]
beta2 <- Pois.sim$BUGSoutput$sims.matrix[,"beta[3]"]
beta3 <- Pois.sim$BUGSoutput$sims.matrix[,"beta[4]"]
PPDK <- Pois.sim$BUGSoutput$sims.matrix[,"PPDK"]
PPDM <- Pois.sim$BUGSoutput$sims.matrix[,"PPDM"]
PPDhot <- Pois.sim$BUGSoutput$sims.matrix[,"PPDhot"]
PPDcold <- Pois.sim$BUGSoutput$sims.matrix[,"PPDcold"]
PPD.8 <- Pois.sim$BUGSoutput$sims.matrix[,"PPD.8"]
PPD1 <- Pois.sim$BUGSoutput$sims.matrix[,"PPD1"]

cat("# of posterior samples:", length(beta0),"/n")
coda::effectiveSize(cbind(beta0, beta1, beta2, beta3))
par(mfrow=c(2,2))
acf(beta0)
acf(beta1)
acf(beta2)
acf(beta3)
par(mfrow=c(1,1))

# trace plots to check convergence
plot.ts(cbind(beta0, beta1, beta2, beta3))
gelman.diag(Pois.sim$BUGSoutput)

# check model fit
GoF_Test <- function(fitted_quantiles) {
  n <- length(fitted_quantiles)
  K <- round((n)^(0.4))
  mK <- table(cut(fitted_quantiles,(0:K)/K))
  np <- n/K
  RB <- sum(((mK-np)^2)/np)
  return(1-pchisq(RB,K-1))
}
GoF <- matrix(NA,ncol=length(temp),nrow=length(beta0))
for (i in 1:length(beta0)) {
  for (j in 1:length(temp)) {
    vals <- ppois(c(y[j]-1,y[j]),exp(beta0[i]+beta1[i]*(temp[j])+beta2[i]*gauge[j]+beta3[i]*tester[j]))
    GoF[i,j] <- runif(1,vals[1],vals[2])
  }
}
GoF_Summary <- apply(GoF,1,GoF_Test)
hist(GoF_Summary,xlim=c(0,1))
mean(GoF_Summary < 0.05)

### Posterior Inference

hist(PPDM-PPDK,main = "Matthew vs Katelyn PPD")
mean(PPDM-PPDK)
quantile(PPDM-PPDK, c(.05, .95))

hist(PPDhot-PPDcold)
quantile(PPDhot-PPDcold, c(.05, .95))

hist(PPD1-PPD.8)
quantile(PPD1-PPD.8, c(.05, .95))


hist(beta0)
mean(beta0)
quantile(beta0,c(.05,.95))

hist(beta1)
mean(beta1)
quantile(beta1,c(.05,.95))



hist(beta2)
mean(beta2)
quantile(beta2,c(.05,.95))


hist(beta3)
mean(beta3)
quantile(beta3,c(.05,.95))

par(mfrow=c(2,2))
hist(beta0)
hist(beta1)
hist(beta2)
hist(beta3)





