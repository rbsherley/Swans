library(jagsUI)
library(lme4)
library(boot)
library(effects)

# read in the data
dat <- read.csv(file="StandardisedWhooperData.csv", header=T)
names(dat)

# Set up some of the variables so they will work in JAGS
dat$SWID.n <- as.numeric(as.factor(dat$SWID))
dat$Season.n <- dat$Season-(min(dat$Season)-1)
dat$WWTNEW.n <- as.numeric(as.factor(dat$WWTNEW))  ## WWTBOTH = 1, WWTIN = 2, WWTNON = 3, WWTOUT = 4, 
dat$Age1 <- dat$Age.Z   # Z score standardised age
dat$Age2 <- (dat$Age.Z*dat$Age.Z)
dat$PBS[1:2] <- 0
dat$PBS.Z[1:2] <- -0.4518592

# Set up some objects for use in prediction
pred.Age <- seq(from=min(dat$Age1),to=max(dat$Age1),length.out=27)
pred.pbs <-c(seq(from=min(na.omit(dat$PBS)),to=max(na.omit(dat$PBS)),length.out=25),max(na.omit(dat$PBS)),max(na.omit(dat$PBS)),max(na.omit(dat$PBS)))
pred.pbs1 <- pred.pbs[2:28]

# This is the minimally adequate model from the frequentist model selection:
global<-glmer(CygnetsIII~
                WWTNEW+Age.Z+I(Age.Z*Age.Z)+PBS.Z+WWTNEW*I(Age.Z*Age.Z)+ 
                (1|Season)+(1|SWID),data=dat,family=binomial,control=glmerControl(optimizer="bobyqa"))
summary(global)


############################################################################
##### JAGS Model
############################################################################

cat(
  "
  model{
  # priors
    b2 ~ dnorm(0, 0.001) # prior for the linear Age effect
    b4 ~ dnorm(0, 0.001) # prior for the PBS (breeding experience) effect 
  
 for(i in 1:M.WWT){
    b1[i] ~ dnorm(0, 0.001)  # priors for the site effect (one level per WWT code)
    b3[i] ~ dnorm(0, 0.001)  # prior for the age-site interaction (one level per WWT code)
 }

# Define the season random effect (M.s = Number of winter seasons)
 for (j in 1:M.s) {
    u.s[j] ~ dnorm(0,tau.u.s)
 }
# Variance for the season random effect
    tau.u.s <- 1/pow(sd.u.s,2)
    sd.u.s ~ dunif(0,5)

# Define the individual swan random effect (M.id = Number of swans)
 for (k in 1:M.id) {
    u.id[k] ~ dnorm(0,tau.u.id)
 }
  # Variance for the individual swan random effect
    tau.u.id <- 1/pow(sd.u.id,2)
    sd.u.id ~ dunif(0,5)
  
  # likelihood
 for (l in 1:n.y){
    y[l] ~ dbern(mu[l])
    logit(mu[l]) <- b1[WWT[l]] + (b2 * age[l]) +  (b3[WWT[l]] * age2[l]) + (b4 * pbs[l]) + u.s[season[l]] + u.id[swan[l]]
    }


# Derived parameters:
for(m in 1:length(pred.Age)){
Pred.pbs[m] ~ dpois(mean.pbs) I(0, pred.pbs[m])
WWTBoth.Pred[m] <- b1[1] + (b2 * pred.Age[m]) +  (b3[1] * (pred.Age[m]*pred.Age[m])) + (b4 * Pred.pbs[m])
WWTIn.Pred[m] <- b1[2] + (b2 * pred.Age[m]) +  (b3[2] * (pred.Age[m]*pred.Age[m])) + (b4 * Pred.pbs[m])
WWTNon.Pred[m] <- b1[3] + (b2 * pred.Age[m]) +  (b3[3] * (pred.Age[m]*pred.Age[m])) + (b4 * Pred.pbs[m])
WWTOut.Pred[m] <- b1[4] + (b2 * pred.Age[m]) +  (b3[4] * (pred.Age[m]*pred.Age[m])) + (b4 * Pred.pbs[m])
}
}  # model
  ", file="swans_model.jags"
)

# bundle the data for input into JAGS
jd <- list(y=dat$CygnetsIII, age=dat$Age1,
          age2=dat$Age2,pbs=dat$PBS,
          n.y = length(dat$CygnetsIII),   
          M.s=max(dat$Season.n), M.id=max(dat$SWID.n),
          WWT=dat$WWTNEW.n, M.WWT=max(dat$WWTNEW.n),
          season=dat$Season.n, swan=dat$SWID.n,
          pred.Age=pred.Age,pred.pbs=pred.pbs1,
          mean.pbs=mean(na.omit(dat$PBS[dat$PBS>0])))

# Parameters monitored
params <- c("b1", "b2", "b3", "b4","sd.u.id","sd.u.s","WWTBoth.Pred","WWTIn.Pred","WWTNon.Pred","WWTOut.Pred","Pred.pbs")

# MCMC settings: Use ni = 550k, nb = 50k and nt = 5 for inference.
ni <- 550000
nt <- 5
nb <- 50000
nc <- 3

# Call JAGS from R 
mod1 <- jags(jd, parallel=T, inits=NULL, params, "swans_model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# Summarize posteriors
print(mod1, digits = 4)
# Save JAGS model object
save(mod1, file="BreedingSuccessModel.RData")
# Save JAGS model output as csv file
write.csv(mod1$summary, file = "Swans_breeding.csv")

## EOF
