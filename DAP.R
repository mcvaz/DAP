### Data Analysis Project

# Loading packages
library(R2jags)
library(lattice)
library(ggplot2)

# Useful functions
load("AddBurnin.RData")
logit = function(x) log(x/(1-x)) # creates logit function
ilogit = function(x) exp(x)/(1+exp(x)) # creates inverse logit function

# Data
nn = matrix(c(71,174,52,53,107,142,49,238,86,103,49,241,357,348,516,512,186,218,36,41,70,80), nrow=11, ncol=2, byrow=TRUE) # number of people in each treatment group in each study
yy = matrix(c(4,1,3,2,12,3,15,19,18,5,11,10,4,9,25,11,20,15,17,25,4,5), nrow=11, ncol=2, byrow=TRUE) # idem, but for the number of deaths after one year of the treatment
yy/nn

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




