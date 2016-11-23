### Data Analysis Project

# Loading packages
library(R2jags)
library(lattice)
library(ggplot2)

# Useful functions
source("AddBurnin.R")
logit = function(x) log(x/(1-x)) # creates logit function
ilogit = function(x) exp(x)/(1+exp(x)) # creates inverse logit function

# Data
namesd = list(c("Brener","Buszman","Chieffo","Makikallio","Palmerini","Sanmartin","Serruys","Seung","Silvestri","White","Wu"),c("Stent","CABG")) # labels for data table
nn = matrix(c(71,174,52,53,107,142,49,238,86,103,49,241,357,348,516,512,186,218,36,41,70,80), nrow=11, ncol=2, byrow=T) # number of people in each treatment group in each study
yy = matrix(c(4,1,3,2,12,3,15,19,18,5,11,10,4,9,25,11,20,15,17,25,4,5), nrow=11, ncol=2, byrow=F) # idem, but for the number of deaths after one year of the treatment
dimnames(nn) = dimnames(yy) = namesd
pp = yy/nn # mortality rates
round(pp,2) 
round(apply(pp,2,mean),2)

# Defining priors


# Jags model
sink("model.txt")
cat("
  model{
    for(i in 1:N){
      for(j in 1:2){
        y[i,j] ~ dbin(pie[i,j],n[i,j])
        logit(pie[i,j]) <- mu[i] + (2*j - 3) * delta[i]/2
        }
      mu[i] ~ dnorm(mu0, tau.mu)
      delta[i] ~ dnorm(delta0, tau.delta)
      }
    mu0 ~ dnorm(0, 1/3.36)
    tau.mu ~ dunif(0, 1/3.36)
    delta0 ~ dnorm(0, 1/2.17)
    tau.delta ~ dunif(0, 1/2.17)
    }    
  ", fill = TRUE)
sink()

# Jags input info
data = list(n=nn, y=yy, N=nrow(nn))
inits = rep(
  list(list(
    pie = matrix(rep(0.5,22), ncol= 2, nrow=11, byrow=TRUE),
    mu = rep(0,11), delta = rep(0,11),
    mu0 = 0, tau.mu = 0.001, delta0 = 0, tau.delta = 0.001))
  ,3)
params = c("pie[1:11,1:2]","mu[1:11]","delta[1:11]","mu0","delta0","tau.mu","tau.delta")

# Jags output
model.out = jags(data=data, inits= NULL, parameter=params, "model.txt", n.chains=3, n.iter=101000, n.burnin=0, n.thin=20, DIC=F)
outparams = dimnames(model.out$BUGSoutput$sims.array)[[3]] # all parameter names estimated by jags

# Take the first 1000 runs as burnin
Output = AddBurnin(model.out$BUGSoutput$sims.array, burnin=1000,n.thin=1)

# Checking for MCMC convergence - 
pdf("acf.pdf",paper="letter")
par(mfrow=c(2,3))
for(i in outparams){
  acf(model.out$BUGSoutput$sims.array[1:5000, 1, i], lag.max= 160, main=i)
  }
dev.off()
par(mfrow=c(1,1))

# Time series
cols = rainbow(3,alpha=0.7)
pdf("timeSeries.pdf",paper="letter")
par(mfrow=c(4,1))
for(i in outparams){
  plot(model.out$BUGSoutput$sims.array[1:1000, 1, i], type="l", col=cols[1], main=i, ylab="", xlab="Iteration")
  lines(model.out$BUGSoutput$sims.array[1:1000, 2, i], type="l", col=cols[2])
  lines(model.out$BUGSoutput$sims.array[1:1000, 3, i], type="l", col=cols[3])  
  }
dev.off()
par(mfrow=c(1,1))

# Results
Output$Burnin.Summary

# Stent vs. CABG (pie[i])

# Stent vs. CABG (pie0)
posts = Output$Burnin.sims.matrix
pp0Stent = ilogit(posts[,"mu0"]-posts[,"delta0"]) # posterior distribution of pie 0 for Stent
pp0CABG = ilogit(posts[,"mu0"]+posts[,"delta0"]) # posterior distribution of pie 0 for CABG
plot(density(pp0Stent),xlab=expression(paste(pi[0])),main="")
lines(density(pp0CABG),lty=2)
legend("topright",legend=c("Stent","CABG"),lty=1:2,bty="n")






