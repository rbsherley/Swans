setwd("/Volumes/GoogleDrive/My Drive/Main_File_Store/Papers/My Papers/In review/PNAS - Whooper Swans/Final_code_and_data")

library(jagsUI)
library(lme4)
library(boot)
library(effects)

dat <- read.csv(file="StandardisedWhooperData.csv", header=T)
names(dat)

dat$SWID.n <- as.numeric(as.factor(dat$SWID))
dat$Season.n <- dat$Season-(min(dat$Season)-1)
dat$WWTNEW.n <- as.numeric(as.factor(dat$WWTNEW))  ## WWTBOTH = 1, WWTIN = 2, WWTNON = 3, WWTOUT = 4, 
dat$Age1 <- dat$Age.Z
dat$Age2 <- (dat$Age.Z*dat$Age.Z)
dat$PBS[1:2] <- 0
dat$PBS.Z[1:2] <- -0.4518592

pred.Age <- seq(from=min(dat$Age1),to=max(dat$Age1),length.out=27)
pred.pbs <-c(seq(from=min(na.omit(dat$PBS)),to=max(na.omit(dat$PBS)),length.out=25),max(na.omit(dat$PBS)),max(na.omit(dat$PBS)),max(na.omit(dat$PBS)))
pred.pbs1 <- pred.pbs[2:28]

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
#params <- c("b1", "b2", "b3", "b4", "b5", "tau.u.id","tau.u.s")

# MCMC settings: Use ni = 550k, nb = 50k and nt = 5 for inference.
ni <- 55000#0
nt <- 5
nb <- 5000#0
nc <- 3

# Call JAGS from R 
mod1 <- jags(jd, parallel=T, inits=NULL, params, "swans_model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# Summarize posteriors
print(mod1, digits = 4)

write.csv(mod1$summary, file = "Swans_breeding_19Mar.csv")

## Extract model fitted proportions and CIs ##
Fit<- effect(term="WWTNEW*I(Age.Z*Age.Z)", mod=global, se = T, confidence.level=0.95, xlevels=30)

FP<- inv.logit(Fit$fit)
UP<- inv.logit(Fit$upper)
LP<- inv.logit(Fit$lower) 
#Extract and plot every 4th number from Fit$fit
#Fit$x displays all the predictions with the factor level
quartz()
par(mfrow=c(1,2))

plot(FP[seq(3,length(FP),4)],type="n",xlab="Age",ylab="Breeding Probability",main="Interaction between WWT Site Switching and Age")
lines(FP[seq(1,length(FP),4)],type="l",lwd="3",col="blue")             #BOTH
lines(FP[seq(2,length(FP),4)],type="l",lwd="3",col="red")             #IN
lines(FP[seq(3,length(FP),4)],type="l",lwd="3",col="green")       #NON   
lines(FP[seq(4,length(FP),4)],type="l",lwd="3",col="purple")  #OUT
boxplot(FP~Fit$x[,1],col="lightblue")

#Include credibility intervals
###### Comparing JAGS and lme4 #######
quartz(width=6,height=6)
par(mar=c(4,4,1.2,1),mfrow=c(2,2))
#WWTBOTH
plot(FP[seq(1,length(FP),4)],type="n",xlab="Age",ylab="Breeding Probability",main="WWT BOTH", ylim=c(0,0.4))
lines(FP[seq(1,length(FP),4)],type="l",lwd="3",col="blue") 
lines(UP[seq(1,length(FP),4)],type="l",lwd="1",col="blue") 
lines(LP[seq(1,length(FP),4)],type="l",lwd="1",col="blue") 
lines(inv.logit(mod1$mean$WWTBoth.Pred))
lines(inv.logit(mod1$q2.5$WWTBoth.Pred))
lines(inv.logit(mod1$q97.5$WWTBoth.Pred))

#WWTIN
plot(FP[seq(2,length(FP),4)],type="n",xlab="Age",ylab="Breeding Probability",main="WWT IN", ylim=c(0,0.4))
lines(FP[seq(2,length(FP),4)],type="l",lwd="3",col="red") 
lines(UP[seq(2,length(FP),4)],type="l",lwd="1",col="red") 
lines(LP[seq(2,length(FP),4)],type="l",lwd="1",col="red")
lines(inv.logit(mod1$mean$WWTIn.Pred))
lines(inv.logit(mod1$q2.5$WWTIn.Pred))
lines(inv.logit(mod1$q97.5$WWTIn.Pred))

#WWTNON
plot(FP[seq(3,length(FP),4)],type="n",xlab="Age",ylab="Breeding Probability",main="WWT NON", ylim=c(0,0.4))
lines(FP[seq(3,length(FP),4)],type="l",lwd="3",col="green") 
lines(UP[seq(3,length(FP),4)],type="l",lwd="1",col="green") 
lines(LP[seq(3,length(FP),4)],type="l",lwd="1",col="green") 
lines(inv.logit(mod1$mean$WWTNon.Pred))
lines(inv.logit(mod1$q2.5$WWTNon.Pred))
lines(inv.logit(mod1$q97.5$WWTNon.Pred))

#WWTOUT
plot(FP[seq(4,length(FP),4)],type="n",xlab="Age",ylab="Breeding Probability",main="WWT OUT", ylim=c(0,0.4))
lines(FP[seq(4,length(FP),4)],type="l",lwd="3",col="purple") 
lines(UP[seq(4,length(FP),4)],type="l",lwd="1",col="purple") 
lines(LP[seq(4,length(FP),4)],type="l",lwd="1",col="purple") 
lines(inv.logit(mod1$mean$WWTOut.Pred))
lines(inv.logit(mod1$q2.5$WWTOut.Pred))
lines(inv.logit(mod1$q97.5$WWTOut.Pred))

quartz.save(file='JAGSvsR.pdf',type='pdf')
dev.off()
################################

###### Comparing JAGS and lme4 #######
quartz(width=6,height=6)
par(mar=c(4,4,1.2,1),mfrow=c(2,2))
#WWTBOTH
plot(inv.logit(mod1$mean$WWTBoth.Pred),type="n",xlab="Age",ylab="Breeding Probability",main="WWT BOTH", ylim=c(0,1))
lines(inv.logit(mod1$mean$WWTBoth.Pred),type="l",lwd="3",col="blue")
lines(inv.logit(mod1$q2.5$WWTBoth.Pred),type="l",lwd="1",col="blue")
lines(inv.logit(mod1$q97.5$WWTBoth.Pred),type="l",lwd="1",col="blue")

#WWTIN
plot(inv.logit(mod1$mean$WWTIn.Pred),type="n",xlab="Age",ylab="Breeding Probability",main="WWT IN", ylim=c(0,1))
lines(inv.logit(mod1$mean$WWTIn.Pred),type="l",lwd="3",col="red")
lines(inv.logit(mod1$q2.5$WWTIn.Pred),type="l",lwd="1",col="red")
lines(inv.logit(mod1$q97.5$WWTIn.Pred),type="l",lwd="1",col="red")

#WWTNON
plot(inv.logit(mod1$mean$WWTNon.Pred),type="n",xlab="Age",ylab="Breeding Probability",main="WWT NON", ylim=c(0,1))
lines(inv.logit(mod1$mean$WWTNon.Pred),type="l",lwd="3",col="green")
lines(inv.logit(mod1$q2.5$WWTNon.Pred),type="l",lwd="1",col="green")
lines(inv.logit(mod1$q97.5$WWTNon.Pred),type="l",lwd="1",col="green")

#WWTOUT
plot(inv.logit(mod1$mean$WWTOut.Pred),type="n",xlab="Age",ylab="Breeding Probability",main="WWT OUT", ylim=c(0,1))
lines(inv.logit(mod1$mean$WWTOut.Pred),type="l",lwd="3",col="purple")
lines(inv.logit(mod1$q2.5$WWTOut.Pred),type="l",lwd="1",col="purple")
lines(inv.logit(mod1$q97.5$WWTOut.Pred),type="l",lwd="1",col="purple")

quartz.save(file='JAGS_final_output_13Mar.pdf',type='pdf')
dev.off()

## EOF

### not used:
test <- rbind(mod1$sims.list$WWTBoth.Pred,mod1$sims.list$WWTBoth.Pred)
test1 <-rbind(mod1$sims.list$WWTNon.Pred,mod1$sims.list$WWTOut.Pred)
WWTBOTH.IN <- WWTNON.OUT <- 0 
for(i in 1:30){
  WWTBOTH.IN[i]=mean(test[,i])
  WWTNON.OUT[i]=mean(test1[,i])
}


plot(inv.logit(mod1$mean$b1)+inv.logit(mod1$mean$b3),ylim=c(0.3,0.6))
segments(seq(1,4,1),(inv.logit(mod1$q2.5$b1)+inv.logit(mod1$q2.5$b3)),seq(1,4,1),(inv.logit(mod1$q97.5$b1)+inv.logit(mod1$q97.5$b3)),col="grey60") 
lab()



quartz(width=9,height=3)
par(mar=c(4,4,0.2,0.2),mfrow=c(1,3))
for(i in 2:4){
temp <- inv.logit(mod1$sims.list$b3[,i])-inv.logit(mod1$sims.list$b3[,1])
plot(density(temp),main='',xlab='',ylab='',cex.axis=1.2)
abline(v=quantile(temp,0.5))
abline(v=quantile(temp,0.025),lty=2)
abline(v=quantile(temp,0.975),lty=2)
abline(v=0,col='red')
mtext(side=1,"Coefficient estimate",line=2.5,cex=1)
mtext(side=2,"Density",line=2.5,cex=1)
mtext(c("a)","b)","c)")[i-1], at = (max(density(temp)$y)*0.98),side=2, line=2.5,cex=1.1,las=1,font=1)
}
quartz.save(file='bs_interaction_effect.pdf',type='pdf')
dev.off()

