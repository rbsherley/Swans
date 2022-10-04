#setwd("~/Google Drive/Main File Store/Papers, Reports, Theses etc/My Papers/In progress/Whooper Swans")

setwd("C:/NMSU/Projects/WhooperSwans")

require(jagsUI)
#require(denstrip)
#require(plotrix)
#require(R2ucare)
#require(unmarked)


## Capture years are 1979 to 2009
CH1 <- as.matrix(read.csv2('Whooper_MS_Age_Encounter_5.csv',header=F,sep=","))
CH <- matrix(data=rep(0,(dim(CH1)[1]*dim(CH1)[2])),nrow = dim(CH1)[1], ncol = dim(CH1)[2])#matrix
for(i in 1:dim(CH)[1]) {CH[i,] <- CH1[i,]}

CH.UCARE <- CH
CH.UCARE[CH==4] <- 1
CH.UCARE[CH==5] <- 2
CH.UCARE[CH==6] <- 3
CH.freq <- rep(1,nrow(CH))

temp <- group_data(CH.UCARE,CH.freq)
CH.freq <- temp[,32]
CH.UCARE <- temp[,1:31]

test3Gsr(CH.UCARE,CH.freq) # transience
test3Gsm(CH.UCARE,CH.freq)
test3Gwbwa(CH.UCARE,CH.freq) # memory
testMitec(CH.UCARE,CH.freq) # trap dependence
testMltec(CH.UCARE,CH.freq) 
overall_JMV(CH.UCARE,CH.freq)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first) ### Tell the model when the time of first capture was

# Recode CH matrix:
# 1 = seen as juvenile @ WWT, 2 = seen as immature @ WWT, 3 = seen as adult @ WWT 
# 4 = seen as juvenile @ non-WWT, 5 = seen as immature @ non-WWT, 6 = seen as adult @ non-WWT 
# 7 = not seen

rCH <- CH
rCH[rCH==0] <- 7


################################################################################################
###### Multistate model with three age-classes and movements rates:
# Multistate model that accounts for trap-dependence (Adapted from Kery & Schaub (2012))
# Trap-dependence modeled in the adults
################################################################################################

# 
sink("swans_M0_no_i_td.jags")
cat("
    model{
    
    #********************************************************
    # Parameters:
    #********************************************************
    # phi: survival probs. age and site-specific
    # phi.jW: juvenile survival @ WWT sites
    # phi.1W: immature survival @ WWT sites
    # phi.aW: adult survival @ WWT sites
    # phi.jN: juvenile survival @ non-WWT sites
    # phi.1N: immature survival @ non-WWT sites
    # phi.aN: adult survival @ non-WWT sites
    #*****************************************************
    # p: recapture probs. age and site-specific
    # p.1W: immature recapture @ WWT sites 
    # p.aW: adult recapture @ WWT sites, given captured at t-1
    # p.aWu: adult recapture @ WWT sites, given not captured at t-1
    # p.1N: immature recapture @ non-WWT sites
    # p.aN: adult recapture @ non-WWT sites, given captured at t-1
    # p.aNu: adult recapture @ non-WWT sites, given not captured at t-1
    #*****************************************************
    # psi: movement probs. age and site-specific
    # psi.jW: juvenile movement from WWT to non-WWT sites
    # psi.jN: juvenile movement from non-WWT sites to WWT sites
    # psi.1W: immature movement from WWT to non-WWT sites
    # psi.1N: immature movement from non-WWT sites to WWT sites
    # psi.aW: adult movement from WWT to non-WWT sites
    # psi.aN: adult movement from non-WWT sites to WWT sites
    #*****************************************************
    # 11 true states (S):
    # 1 = Juvenile alive @ WWT, 
    # 2 = Immature alive, seen at t-1 @ WWT, 
    # 3 = Adult alive, seen at t-1 @ WWT, 
    # 4 = Juvenile alive @ non-WWT, 
    # 5 = Immature alive, seen at t-1 @ non-WWT, 
    # 6 = Adult alive, seen at t-1 @ non-WWT 
    # 7 = Dead
    # 8 = Immature alive, not seen at t-1 @ WWT, 
    # 9 = Adult alive, not seen at t-1 @ WWT 
    # 10 = Immature alive, not seen at t-1 @ non-WWT, 
    # 11 = Adult alive, not seen at t-1 @ non-WWT, 

    #*****************************************************
    # 7 observation states (O):
    # 1 = seen as juvenile @ WWT, 
    # 2 = seen as immature @ WWT, 
    # 3 = seen as adult @ WWT 
    # 4 = seen as juvenile @ non-WWT, 
    # 5 = seen as immature @ non-WWT, 
    # 6 = seen as adult @ non-WWT 
    # 7 = not seen
    #*****************************************************
    
    #******************************
    # Priors and constraints
    #******************************

 for(t in 1:(n.occasions-1)){
    
    # Survival
    phi.jW[t] <- mean.phi[1]
    phi.1W[t] <- mean.phi[2]
    phi.aW[t] <- mean.phi[3] 
    
    phi.jN[t] <- mean.phi[4]
    phi.1N[t] <- mean.phi[5] 
    phi.aN[t] <- mean.phi[6]
    
    # Recapture
    p.1W[t] <- mean.p[1] 
    p.aW[t] <- mean.p[2] 
    p.aWu[t] <- mean.p[3] 

    p.1N[t] <- mean.p[4]
    p.aN[t] <- mean.p[5]
    p.aNu[t] <- mean.p[6]

    # Movement between WWT and Non-WWT (juv(1), 1yr(2), Adu(a))
    
    psi.jW[t] <- mean.psi[1]
    psi.jN[t] <- mean.psi[2]

    psi.1W[t] <- mean.psi[3]
    psi.1N[t] <- mean.psi[4] 

    psi.aW[t] <- mean.psi[5]
    psi.aN[t] <- mean.psi[6]
  } #t

  for (k in 1:6){
    mean.phi[k] ~ dunif(0,1)  # priors for mean state-specific survival
    mean.psi[k] ~ dunif(0,1) # priors for mean state-specific transitions
  }#k
  for (k in 1:6){
    mean.p[k] ~ dunif(0,1) # priors for mean state-specific recapture
  }#k

    #**************************   
    # Define state-transition and observation matrices
    #**************************
    # Define probabilities of state S(t+1) given S(t)
   for(t in 1:(n.occasions-1)){  # loop over time
    # First index = states at time t-1, last index = states at time t
    ps[1,t,1] <- 0
    ps[1,t,2] <- phi.jW[t]*(1-psi.jW[t])*p.1W[t] 
    ps[1,t,3] <- 0
    ps[1,t,4] <- 0
    ps[1,t,5] <- phi.jW[t]*psi.jW[t]*p.1N[t]
    ps[1,t,6] <- 0 
    ps[1,t,7] <- 1-phi.jW[t]
    ps[1,t,8] <- phi.jW[t]*(1-psi.jW[t])*(1-p.1W[t])
    ps[1,t,9] <- 0
    ps[1,t,10] <- phi.jW[t]*psi.jW[t]*(1-p.1N[t])
    ps[1,t,11] <- 0

    ps[2,t,1] <- 0
    ps[2,t,2] <- 0
    ps[2,t,3] <- phi.1W[t]*(1-psi.1W[t])*p.aW[t] 
    ps[2,t,4] <- 0
    ps[2,t,5] <- 0
    ps[2,t,6] <- phi.1W[t]*psi.1W[t]*p.aN[t] 
    ps[2,t,7] <- 1-phi.1W[t]
    ps[2,t,8] <- 0 
    ps[2,t,9] <- phi.1W[t]*(1-psi.1W[t])*(1-p.aW[t]) 
    ps[2,t,10] <- 0
    ps[2,t,11] <- phi.1W[t]*psi.1W[t]*(1-p.aN[t])
       
    ps[3,t,1] <- 0
    ps[3,t,2] <- 0
    ps[3,t,3] <- phi.aW[t]*(1-psi.aW[t])*p.aW[t] 
    ps[3,t,4] <- 0
    ps[3,t,5] <- 0
    ps[3,t,6] <- phi.aW[t]*psi.aW[t]*p.aN[t]
    ps[3,t,7] <- 1-phi.aW[t] 
    ps[3,t,8] <- 0 
    ps[3,t,9] <- phi.aW[t]*(1-psi.aW[t])*(1-p.aW[t]) 
    ps[3,t,10] <- 0
    ps[3,t,11] <- phi.aW[t]*psi.aW[t]*(1-p.aN[t])

    ps[4,t,1] <- 0
    ps[4,t,2] <- phi.jN[t]*psi.jN[t]*p.1W[t]
    ps[4,t,3] <- 0 
    ps[4,t,4] <- 0
    ps[4,t,5] <- phi.jN[t]*(1-psi.jN[t])*p.1N[t] 
    ps[4,t,6] <- 0
    ps[4,t,7] <- 1-phi.jN[t]
    ps[4,t,8] <- phi.jN[t]*psi.jN[t]*(1-p.1W[t])
    ps[4,t,9] <- 0
    ps[4,t,10] <- phi.jN[t]*(1-psi.jN[t])*(1-p.1N[t])
    ps[4,t,11] <- 0 

    ps[5,t,1] <- 0
    ps[5,t,2] <- 0
    ps[5,t,3] <- phi.1N[t]*psi.1N[t]*p.aW[t]
    ps[5,t,4] <- 0
    ps[5,t,5] <- 0
    ps[5,t,6] <- phi.1N[t]*(1-psi.1N[t])*p.aN[t]
    ps[5,t,7] <- 1-phi.1N[t]
    ps[5,t,8] <- 0
    ps[5,t,9] <- phi.1N[t]*psi.1N[t]*(1-p.aW[t]) 
    ps[5,t,10] <- 0 
    ps[5,t,11] <- phi.1N[t]*(1-psi.1N[t])*(1-p.aN[t]) 

    ps[6,t,1] <- 0
    ps[6,t,2] <- 0
    ps[6,t,3] <- phi.aN[t]*psi.aN[t]*p.aW[t] 
    ps[6,t,4] <- 0
    ps[6,t,5] <- 0
    ps[6,t,6] <- phi.aN[t]*(1-psi.aN[t])*p.aN[t]  
    ps[6,t,7] <- 1-phi.aN[t]
    ps[6,t,8] <- 0
    ps[6,t,9] <- phi.aN[t]*psi.aN[t]*(1-p.aW[t]) 
    ps[6,t,10] <- 0 
    ps[6,t,11] <- phi.aN[t]*(1-psi.aN[t])*(1-p.aN[t]) 

    ps[7,t,1] <- 0
    ps[7,t,2] <- 0       
    ps[7,t,3] <- 0     
    ps[7,t,4] <- 0
    ps[7,t,5] <- 0     
    ps[7,t,6] <- 0 
    ps[7,t,7] <- 1
    ps[7,t,8] <- 0       
    ps[7,t,9] <- 0     
    ps[7,t,10] <- 0
    ps[7,t,11] <- 0     

    ps[8,t,1] <- 0
    ps[8,t,2] <- 0 
    ps[8,t,3] <- phi.1W[t]*(1-psi.1W[t])*p.aWu[t]
    ps[8,t,4] <- 0 
    ps[8,t,5] <- 0 
    ps[8,t,6] <- phi.1W[t]*psi.1W[t]*p.aNu[t]
    ps[8,t,7] <- 1-phi.1W[t] 
    ps[8,t,8] <- 0
    ps[8,t,9] <- phi.1W[t]*(1-psi.1W[t])*(1-p.aWu[t])
    ps[8,t,10] <- 0
    ps[8,t,11] <- phi.1W[t]*psi.1W[t]*(1-p.aNu[t])

    ps[9,t,1] <- 0
    ps[9,t,2] <- 0 
    ps[9,t,3] <- phi.aW[t]*(1-psi.aW[t])*p.aWu[t] 
    ps[9,t,4] <- 0 
    ps[9,t,5] <- 0
    ps[9,t,6] <- phi.aW[t]*psi.aW[t]*p.aNu[t]
    ps[9,t,7] <- 1-phi.aW[t] 
    ps[9,t,8] <- 0
    ps[9,t,9] <- phi.aW[t]*(1-psi.aW[t])*(1-p.aWu[t]) 
    ps[9,t,10] <- 0
    ps[9,t,11] <- phi.aW[t]*psi.aW[t]*(1-p.aNu[t])
    
    ps[10,t,1] <- 0
    ps[10,t,2] <- 0
    ps[10,t,3] <- phi.1N[t]*psi.1N[t]*p.aWu[t] 
    ps[10,t,4] <- 0
    ps[10,t,5] <- 0
    ps[10,t,6] <- phi.1N[t]*(1-psi.1N[t])*p.aNu[t]
    ps[10,t,7] <- 1-phi.1N[t] 
    ps[10,t,8] <- 0
    ps[10,t,9] <- phi.1N[t]*psi.1N[t]*(1-p.aWu[t]) 
    ps[10,t,10] <- 0 
    ps[10,t,11] <- phi.1N[t]*(1-psi.1N[t])*(1-p.aNu[t]) 

    ps[11,t,1] <- 0
    ps[11,t,2] <- 0
    ps[11,t,3] <- phi.aN[t]*psi.aN[t]*p.aWu[t] 
    ps[11,t,4] <- 0
    ps[11,t,5] <- 0
    ps[11,t,6] <- phi.aN[t]*(1-psi.aN[t])*p.aNu[t] 
    ps[11,t,7] <- 1-phi.aN[t]
    ps[11,t,8] <- 0
    ps[11,t,9] <- phi.aN[t]*psi.aN[t]*(1-p.aWu[t]) 
    ps[11,t,10] <- 0 
    ps[11,t,11] <- phi.aN[t]*(1-psi.aN[t])*(1-p.aNu[t]) 
    
    # Define probabilities of O(t) given S(t)
    # First index = states at time t, last index = observations at time t
    
    po[1,t,1] <- 0
    po[1,t,2] <- 0
    po[1,t,3] <- 0
    po[1,t,4] <- 0
    po[1,t,5] <- 0
    po[1,t,6] <- 0
    po[1,t,7] <- 1

    po[2,t,1] <- 0
    po[2,t,2] <- 1
    po[2,t,3] <- 0
    po[2,t,4] <- 0
    po[2,t,5] <- 0
    po[2,t,6] <- 0
    po[2,t,7] <- 0

    po[3,t,1] <- 0
    po[3,t,2] <- 0
    po[3,t,3] <- 1
    po[3,t,4] <- 0
    po[3,t,5] <- 0
    po[3,t,6] <- 0
    po[3,t,7] <- 0

    po[4,t,1] <- 0
    po[4,t,2] <- 0
    po[4,t,3] <- 0
    po[4,t,4] <- 0
    po[4,t,5] <- 0
    po[4,t,6] <- 0
    po[4,t,7] <- 1

    po[5,t,1] <- 0
    po[5,t,2] <- 0
    po[5,t,3] <- 0
    po[5,t,4] <- 0
    po[5,t,5] <- 1
    po[5,t,6] <- 0
    po[5,t,7] <- 0

    po[6,t,1] <- 0
    po[6,t,2] <- 0
    po[6,t,3] <- 0
    po[6,t,4] <- 0
    po[6,t,5] <- 0
    po[6,t,6] <- 1
    po[6,t,7] <- 0

    po[7,t,1] <- 0
    po[7,t,2] <- 0
    po[7,t,3] <- 0
    po[7,t,4] <- 0
    po[7,t,5] <- 0
    po[7,t,6] <- 0
    po[7,t,7] <- 1
    
    po[8,t,1] <- 0
    po[8,t,2] <- 0
    po[8,t,3] <- 0
    po[8,t,4] <- 0
    po[8,t,5] <- 0
    po[8,t,6] <- 0
    po[8,t,7] <- 1

    po[9,t,1] <- 0
    po[9,t,2] <- 0
    po[9,t,3] <- 0
    po[9,t,4] <- 0
    po[9,t,5] <- 0
    po[9,t,6] <- 0
    po[9,t,7] <- 1

    po[10,t,1] <- 0
    po[10,t,2] <- 0
    po[10,t,3] <- 0
    po[10,t,4] <- 0
    po[10,t,5] <- 0
    po[10,t,6] <- 0
    po[10,t,7] <- 1

    po[11,t,1] <- 0
    po[11,t,2] <- 0
    po[11,t,3] <- 0
    po[11,t,4] <- 0
    po[11,t,5] <- 0
    po[11,t,6] <- 0
    po[11,t,7] <- 1
  } #t
  
    #***********************************
    # Define the state-space likelihood
    #***********************************
 for(i in 1:nind){
     z[i,f[i]] <- y[i,f[i]]
  for(t in (f[i]+1):n.occasions){  # loop over time
    # State equation: draw S(t) given S(t-1)
     z[i,t] ~ dcat(ps[z[i,t-1],t-1,])
    # Observation equation: draw O(t) given S(t)
     y[i,t] ~ dcat(po[z[i,t],t-1,])
  } #t
 } #i  
    
    #**************************   
    # Derived fidelity rates
    #**************************    
    fid.jWW <- (1-mean.psi[1])
    fid.jNN<- (1-mean.psi[2])
    fid.1WW <- (1-mean.psi[3])
    fid.1NN<- (1-mean.psi[4])
    fid.aWW <- (1-mean.psi[5])
    fid.aNN<- (1-mean.psi[6])

    }#model
    ",fill=TRUE)
sink()

# Function to create known latent states z

known.state.ms <- function(ms, notseen){
  #notseen: label for 'not seen'
  state <- ms
  state[state==notseen] <- NA
  for (i in 1:dim(ms)[1]){
    m <- min(which(!is.na(state[i,])))
    state[i,m] <- NA
  }
  return(state)
}   

# Function to create a matrix of initial values for unknown z
age.init <- function(ch,f){
  age <- array(NA, dim=dim(ch))
  
  for (i in 1:nrow(ch)){
    
    # WWT-juvenile  
    if(ch[i,f[i]]==1) {
      age[i,f[i]] <- 1
      if(f[i]==ncol(ch)) next  
      age[i,(f[i]+1)] <- sample(c(8,10),1)
      
      if((f[i]+1)==ncol(ch)) next   
      for (t in (f[i]+2):ncol(ch)){
        age[i,t] <- sample(c(9,11),1)
      }#t  
    } 
    
    # WWT-immature 
    if(ch[i,f[i]]==2){
      age[i,f[i]] <- 2 
      if(f[i]==ncol(ch)) next
      for (t in (f[i]+1):ncol(ch)){
        age[i,t] <- sample(c(9,11),1)
      }
    } 
    
    #WWT-adult
    if(ch[i,f[i]]==3){
      age[i,f[i]] <- 3 
      if(f[i]==ncol(ch)) next
      for (t in (f[i]+1):ncol(ch)){
        age[i,t] <- sample(c(9,11),1)
      }
    }
    
    # Non-WWT-juvenile
    if(ch[i,f[i]]== 4) {
      age[i,f[i]] <- 4
      if(f[i]==ncol(ch)) next  
      age[i,(f[i]+1)] <- sample(c(8,10),1)
      
      if((f[i]+1)==ncol(ch)) next   
      for (t in (f[i]+2):ncol(ch)){
        age[i,t] <- sample(c(9,11),1)
      }#t  
    }
    
    # Non-WWT-immature 
    if(ch[i,f[i]]==5){
      age[i,f[i]] <- 5 
      if(f[i]==ncol(ch)) next
      for (t in (f[i]+1):ncol(ch)){
        age[i,t] <- sample(c(9,11),1)
      }
    }
    
    #Non-WWT-adult
    if(ch[i,f[i]]==6){
      age[i,f[i]] <- 6 
      if(f[i]==ncol(ch)) next
      for (t in (f[i]+1):ncol(ch)){
        age[i,t] <- sample(c(9,11),1)
      }
    }
  }#i
  ini <- array(ch, dim=dim(ch))
  for (i in 1:nrow(ch)){
    for (t in f[i]:ncol(ch)){
      if(ch[i,t]==7) ini[i,t] <- age[i,t]
      if(ch[i,t]!=7) ini[i,t] <- NA
      ini[i,1:f[i]] <- NA
    }
  }

  return(ini)
}   #function


# Initial values
inits <- function(){list(z = age.init(rCH,f),
                         mean.phi = runif(6, 0, 1), 
                         mean.psi = runif(6, 0, 1),
                         mean.p = runif(6, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p", "mean.psi",
                "fid.jWW","fid.jNN","fid.1WW",
                "fid.1NN","fid.aWW","fid.aNN")

# MCMC settings -  Computing time = 2.138683 days !!!!
ni <- 15000
nt <- 10
nb <- 5000
nc <- 3

# Bundle data
jags.data <- list(y = rCH, f = f, nind = dim(rCH)[1], n.occasions = dim(rCH)[2], z= known.state.ms(rCH,7)) 

# Call JAGS from R 

start <-Sys.time() 
start
ws.ms.1.td <- jags(jags.data, inits, parameters, "swans_M0_no_i_td.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)
end <-Sys.time() 
end
end-start

# Summarize posteriors
print(ws.ms.1.td, digits = 3)

save(ws.ms.1.td, file="ws_ms_1_td.RData")

write.csv(ws.ms.1.td$summary, file = "WA_S(a)_M(a)_E(a)_td.csv")

#************************
# Figures
#************************

pdf(file="fig1-Survival.pdf")
op <- layout(matrix(c(1:6), nrow=2, byrow=T))
xl <- c("Ju.phi.WWT","Im.phi.WWT","Ad.phi.WWT","Ju.phi.non-WWT","Im.phi.non-WWT","Ad.phi.non-WWT")
for(k in 1:6)
{
  hist(ws.ms.1.td$sims.list$mean.phi[,k],main="",xlab=xl[k]) 
  abline(v=ws.ms.1.td$mean$mean.phi[k],col="red",lty=2,lwd=2)
}
par(op)
dev.off()

pdf(file="fig1-Transition.pdf")
op <- layout(matrix(c(1:6), nrow=2, byrow=F))
xl <- c("Ju.psi.WWT","Ju.psi.non-WWT","Im.psi.WWT","Im.psi.non-WWT","Ad.psi.WWT","Ad.psi.non-WWT")
for(k in 1:6)
{
  hist(ws.ms.1.td$sims.list$mean.psi[,k],main="",xlab=xl[k]) 
  abline(v=ws.ms.1.td$mean$mean.psi[k],col="red",lty=2,lwd=2)
}
par(op)
dev.off()

pdf(file="fig1-Recapture.pdf")
op <- layout(matrix(c(1:6), nrow=2, byrow=T))
xl <- c("Im.p.WWT","Ad.TrapA.p.WWT","Ad.TrapU.p.WWT","Im.p.non-WWT","Ad.TrapA.p.non-WWT","Ad.TrapU.p.non-WWT")
for(k in 1:6)
{
  hist(ws.ms.1.td$sims.list$mean.p[,k],main="", xlab=xl[k]) 
  abline(v=ws.ms.1.td$mean$mean.p[k],col="red",lty=2,lwd=2)
}
par(op)
dev.off()

pdf(file="fig1-Fidelity.pdf")
op <- layout(matrix(c(1:6), nrow=2, byrow=T))
hist(ws.ms.1.td$sims.list$fid.jWW,main="",xlab="fid.jWW")
abline(v=ws.ms.1.td$mean$fid.jWW,col="red",lty=2,lwd=2)
hist(ws.ms.1.td$sims.list$fid.1WW,main="",xlab="fid.1WW")
abline(v=ws.ms.1.td$mean$fid.1WW,col="red",lty=2,lwd=2)
hist(ws.ms.1.td$sims.list$fid.aWW,main="",xlab="fid.aWW")
abline(v=ws.ms.1.td$mean$fid.aWW,col="red",lty=2,lwd=2)
hist(ws.ms.1.td$sims.list$fid.jNN,main="",xlab="fid.jNN")
abline(v=ws.ms.1.td$mean$fid.jNN,col="red",lty=2,lwd=2)
hist(ws.ms.1.td$sims.list$fid.1NN,main="",xlab="fid.1NN")
abline(v=ws.ms.1.td$mean$fid.1NN,col="red",lty=2,lwd=2)
hist(ws.ms.1.td$sims.list$fid.aNN,main="",xlab="fid.aNN")
abline(v=ws.ms.1.td$mean$fid.aNN,col="red",lty=2,lwd=2)
par(op)
dev.off()

##################################################################################################

# Plots 
Image <- function(){
  par(mfrow=c(3,1),mar=c(3, 3, 0.2, 1.5),mgp = c(2, 0.6, 0))
  plotCI(x=1:6,y=ws.ms.1$mean$mean.phi[1:6],ui=ws.ms.1$q97.5$mean.phi[1:6],li=ws.ms.1$q2.5$mean.phi[i],xaxt='n',xlim=c(1,6),ylim=c(0.55,1.05),xlab='Site and age-class',ylab='Survival',main='',pch=16, col='white',lty=2)
  for(i in 1:6){
  denstrip(density(ws.ms.1$sims.list$mean.phi[,i])$x,density(ws.ms.1$sims.list$mean.phi[,i])$y, at=i, horiz=FALSE,colmax='grey20',
           mticks=c(ws.ms.1$mean$mean.phi[i]),mlen=1.5,mwd=1.5,mcol='black',
           ticks=c(ws.ms.1$q2.5$mean.phi[i],ws.ms.1$q97.5$mean.phi[i]),
           tcol=c('grey70'),tlen=1, twd=1.5,scale=0.8, colmin = 'grey95')
  }#i
  axis(1, at=1:6, labels=c("Ju. WWT","Im. WWT","Ad. WWT","Ju. non-WWT","Im. non-WWT","Ad. non-WWT"))
  text("Survival",x=3.5,y=1.04,cex=1.1)
  points(x=1,y=1.04,pch=1,col='black',cex=2.3)
  text("A",x=1,y=1.04,cex=0.8)
  
  plotCI(x=1:4,y=ws.ms.1$mean$mean.p[1:4],ui=ws.ms.1$q97.5$mean.p[1:4],li=ws.ms.1$q2.5$mean.p[i],xaxt='n',xlim=c(1,4),ylim=c(0.45,1.05),xlab='Site and age-class',ylab='Encounter rate',main='',pch=16, col='white',lty=2)
  for(i in 1:4){
    denstrip(density(ws.ms.1$sims.list$mean.p[,i])$x,density(ws.ms.1$sims.list$mean.p[,i])$y, at=i, horiz=FALSE,colmax='grey20',
             mticks=c(ws.ms.1$mean$mean.p[i]),mlen=1.5,mwd=1.5,mcol='black',
             ticks=c(ws.ms.1$q2.5$mean.p[i],ws.ms.1$q97.5$mean.p[i]),
             tcol=c('grey70'),tlen=1, twd=1.5,scale=0.8, colmin = 'grey95')
  }#i
  
  axis(1, at=1:4, labels=c("Im. WWT","Ad. WWT","Im. non-WWT","Ad. non-WWT"))
  text("Encounter",x=2.5,y=1.04,cex=1.1)
  points(x=1,y=1.04,pch=1,col='black',cex=2.3)
  text("B",x=1,y=1.04,cex=0.8)
  
  plotCI(x=1:6,y=ws.ms.1$mean$mean.psi[1:6],ui=ws.ms.1$q97.5$mean.psi[1:6],li=ws.ms.1$q2.5$mean.psi[i],xaxt='n',xlim=c(1,6),ylim=c(0,0.55),xlab='Site and age-class',ylab='Movement rate',main='',pch=16, col='white',lty=2)
  for(i in 1:6){
    denstrip(density(ws.ms.1$sims.list$mean.psi[,i])$x,density(ws.ms.1$sims.list$mean.psi[,i])$y, at=i, horiz=FALSE,colmax='grey20',
             mticks=c(ws.ms.1$mean$mean.psi[i]),mlen=1.5,mwd=1.5,mcol='black',
             ticks=c(ws.ms.1$q2.5$mean.psi[i],ws.ms.1$q97.5$mean.psi[i]),
             tcol=c('grey70'),tlen=1, twd=1.5,scale=0.8, colmin = 'grey95')
  }#i
  
  axis(1, at=1:6, labels=c("Ju. W-N","Ju. N-W", "Im. W-N","Im. N-W", "Ad. W-N","Ad. N-W"))
  text("Movement",x=3.5,y=.54,cex=1.1)
  points(x=1,y=0.54,pch=1,col='black',cex=2.3)
  text("C",x=1,y=0.54,cex=0.8)
  
}# image
pdf(file = "W.Swans.pdf", width = 5, height = 8, onefile = TRUE,pagecentre=TRUE,colormodel="rgb")
Image()
dev.off()
jpeg("W.Swans.jpg",width=5,height=8,units='in',res=350)
Image()
dev.off()
  