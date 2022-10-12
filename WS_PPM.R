require(jagsUI)
require(denstrip)
require(plotrix)
require(boot)
# #####################################################
# Demographic model for Whooper Swans in JAGS
# #####################################################

# Set max age class for JAGS:
maxage <-20

# Set number of years of projections into the future (post 2014 data):
pYr <- 16

# Observed population counts
ObsPop <- read.csv("PopData.csv")
TP <- (ObsPop$TP)/2
WWTP <- (ObsPop$WWTP)/2
Yr <- seq(from=min(ObsPop$YearStart),to=(max(ObsPop$YearStart)+pYr),by=1)


###############################################
### RUN THIS COMMENTED SECTION ONLY IF READING IN OUTPUTS
### FROM YOUR OWN RUNS OF THE CMR AND PRODUCTIVITY
### MODELS. SKIP THIS SECTION IF USING THE POSTERIORS
### PROVIDED IN "SwansAllParamsFinal.RData"
###############################################

#load("ws_ms_1_td.RData") # load in the CMR (survival and transition) model posteriors
#load("BreedingSuccessModel.RData") # load in the breeding probability model posteriors

### SET VECTORS WITH DEMOGRAPHIC PARAMETRES ###
# #Survival rates

#SjuvWWT <- ws.ms.1.td$sims.list$mean.phi[,1] #Juvenile survival inside WWT
#SjuvNN <- ws.ms.1.td$sims.list$mean.phi[,4] #Juvenile survival outside WWT

#S1yWWT <- ws.ms.1.td$sims.list$mean.phi[,2] #1-year old survival inside WWT
#S1yNN <- ws.ms.1.td$sims.list$mean.phi[,5] #1-year old survival outside WWT
#SadWWT <- ws.ms.1.td$sims.list$mean.phi[,3] #Adult survival inside WWT
#SadNN <- ws.ms.1.td$sims.list$mean.phi[,6] #Adult survival outside WWT

# #Transitions

#TjuvWWT_1yNN <- ws.ms.1.td$sims.list$mean.psi[,1] #Juvenile WWT to 1-year old non-WWT 
#TjuvWWT_1yWWT <- ws.ms.1.td$sims.list$fid.jWW #Juvenile WWT to 1-year old WWT
#TjuvNN_1yWWT <- ws.ms.1.td$sims.list$mean.psi[,2] #Juvenile non-WWT to 1-year old WWT
#TjuvNN_1yNN <- ws.ms.1.td$sims.list$fid.jNN #Juvenile non-WWT to 1-year old non-WWT

#T1yWWT_adNN <- ws.ms.1.td$sims.list$mean.psi[,3] # First-year WWT to adult non-WWT
#T1yWWT_adWWT <- ws.ms.1.td$sims.list$fid.1WW # First-year WWT to adult WWT
# T1yNN_adWWT <-  ws.ms.1.td$sims.list$mean.psi[,4] # First-year non-WWT to adult WWT
# T1yNN_adNN <- ws.ms.1.td$sims.list$fid.1NN # First-year non-WWT to adult non-WWT
# 
# TadWWT_adNN <- ws.ms.1.td$sims.list$mean.psi[,5] #Adult WWT to adult non-WWT 
# TadWWT_adWWT <- ws.ms.1.td$sims.list$fid.aWW #Adult WWT to adult WWT 
# TadNN_adWWT <- ws.ms.1.td$sims.list$mean.psi[,6] #Adult non-WWT to adult WWT
# TadNN_adNN <- ws.ms.1.td$sims.list$fid.aNN #Adult non-WWT to adult non-WWT
# 
# #Breeding probability
# 
# PWWTBoth <- exp(mod1$sims.list$WWTBoth.Pred)/(1+exp(mod1$sims.list$WWTBoth.Pred))
# PWWTIn <- exp(mod1$sims.list$WWTIn.Pred)/(1+exp(mod1$sims.list$WWTIn.Pred))
# PNNBoth <- exp(mod1$sims.list$WWTNon.Pred)/(1+exp(mod1$sims.list$WWTNon.Pred))
# PWWTOut <- exp(mod1$sims.list$WWTOut.Pred)/(1+exp(mod1$sims.list$WWTOut.Pred))


# Setting up the parameters (mean and sd of breeding probabilities) for priors in the JAGS model (these are used in jd when passing data to JAGS below)
P.m <- P.sd <- P.L <- P.U <- PW.m <- PW.sd <- PW.L <- PW.U <- PIn.m <- PIn.sd <- PIn.U <- PIn.L <- PWIn.m <- PWIn.sd <- PWIn.L <- PWIn.U <- 0
for(i in 1:maxage){
  P.m[i] <- mean(PNNBoth[,i])
  P.sd[i] <- sd(PNNBoth[,i])
  P.L[i] <- quantile(PNNBoth[,i],0.025)
  P.U[i] <- quantile(PNNBoth[,i],0.975)
  
  PW.m[i] <- mean(PWWTBoth[,i])
  PW.sd[i] <- sd(PWWTBoth[,i])
  PW.L[i] <- quantile(PWWTBoth[,i],0.025)
  PW.U[i] <- quantile(PWWTBoth[,i],0.975)
  
  PIn.m[i] <- mean(PWWTOut[,i])
  PIn.sd[i] <- sd(PWWTOut[,i])
  PIn.L[i] <- quantile(PWWTOut[,i],0.025)
  PIn.U[i] <- quantile(PWWTOut[,i],0.975)
  
  PWIn.m[i] <- mean(PWWTIn[,i])
  PWIn.sd[i] <- sd(PWWTIn[,i])
  PWIn.L[i] <- quantile(PWWTIn[,i],0.025)
  PWIn.U[i] <- quantile(PWWTIn[,i],0.975)
}

#rm(mod1)
#rm(ws.ms.1.td)

# Function to translate mean and sd of rates into parameters of beta distributions (used in jd when passing data to JAGS below)
beta.parameters <- function(x, sd.x){
  u <- x*(1-x)/sd.x^2-1
  alpha <- x*u
  beta <- (1-x)*u
  return(list(alpha = alpha, beta = beta))
}

#Productivity
BSData <- read.csv("BSData.csv")
Prod <- BSData$BS

### parameters to specify a gamma prior for productivity assuming the parameterization gamma(shape,rate): shape = mean^2/sd^2, rate = mean/sd^2 (Kruschke 2015)
mp <- mean(Prod)
mpsd <- sd(Prod) # # derive the variance based on the central limit theorm
mp.shape <- ((mp)^2/(mpsd)^2)
mp.rate <- mp/(mpsd)^2

# ----------------------------------------
#  Non-WWT Population:
# ----------------------------------------
# Deterministic LMM for comparison
sa <- mean(SadNN)
si <- mean(S1yNN)
sj <- mean(SjuvNN)

P <- colMeans(PNNBoth)            # proportion of adults breeding
R  <- mean(Prod) 

mat <- matrix(c(0,((P[2]*R))/2*sa,((P[3]*R))/2*sa,((P[4]*R))/2*sa,((P[5]*R))/2*sa,((P[6]*R))/2*sa,((P[7]*R))/2*sa,((P[8]*R))/2*sa,((P[9]*R))/2*sa,((P[10]*R))/2*sa,((P[11]*R))/2*sa,((P[12]*R))/2*sa,((P[13]*R))/2*sa,((P[14]*R))/2*sa,((P[15]*R))/2*sa,((P[16]*R))/2*sa,((P[17]*R))/2*sa,((P[18]*R))/2*sa,((P[19]*R))/2*sa,((P[20]*R))/2*sa,
                sj,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,si,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,sa,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,sa,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,sa,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,sa,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,sa,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,sa,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,sa,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,sa,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,sa,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,0,sa,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,0,0,sa,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,0,0,0,sa,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,sa,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,sa,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,sa,0,0,0,
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,sa,0,0,
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,sa,sa), nrow = 20, ncol = 20, byrow = T)
Roots <- eigen(mat)
lambda <- Re(Roots$values[1])
w <- Roots$vectors[, 1]
V <- Conj(solve(Roots$vectors))
v <- Re(t(V)[, 1])
senmat <- Re(v %*% t(w))
emat <- senmat * mat/lambda

lambda                 # population growth rate
senmat                 # sensitivity matrix
emat                  # elasticity matrix

# stable age distribution
sage <- as.numeric(w/sum(w))       
op <- par(mfrow=c(1,2))
barplot(sage, names.arg = c(0:18, "ad"), las=1, ylab = "Proportion in population", xlab = "Age class")    
barplot(c(sum(sage[1]), sum(sage[2:maxage])), names.arg = c("juv", "ad"), las=1, ylab = "Proportion in population", xlab = "Age class")    
par(op)

# Stage-specific population size in first year:
SSPop <- round((sage/sum(sage[1:maxage]))*TP[1],0)



# ----------------------------------------
#  WWT Population:
# ----------------------------------------
# Deterministic LMM for comparison
saW <- mean(SadWWT)
siW <- mean(S1yWWT)
sjW <- mean(SjuvWWT)

PW <- colMeans(PWWTBoth)            # proportion of adults breeding
RW  <- mean(Prod) 

matW <- matrix(c(0,((PW[2]*RW))/2*saW,((PW[3]*RW))/2*saW,((PW[4]*RW))/2*saW,((PW[5]*RW))/2*saW,((PW[6]*RW))/2*saW,((PW[7]*RW))/2*saW,((PW[8]*RW))/2*saW,((PW[9]*RW))/2*saW,((PW[10]*RW))/2*saW,((PW[11]*RW))/2*saW,((PW[12]*RW))/2*saW,((PW[13]*RW))/2*saW,((PW[14]*RW))/2*saW,((PW[15]*RW))/2*saW,((PW[16]*RW))/2*saW,((PW[17]*RW))/2*saW,((PW[18]*RW))/2*saW,((PW[19]*RW))/2*saW,((PW[20]*RW))/2*saW,
                 sjW,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,siW,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,saW,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,saW,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,saW,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,saW,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,saW,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,saW,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,saW,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,saW,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,saW,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,0,saW,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,0,0,saW,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,saW,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,saW,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,saW,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,saW,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,saW,0,0,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,saW,saW), nrow = 20, ncol = 20, byrow = T)
RootsW <- eigen(matW)
lambdaW <- Re(RootsW$values[1])
wW <- RootsW$vectors[, 1]
VW <- Conj(solve(RootsW$vectors))
vW <- Re(t(VW)[, 1])
senmatW <- Re(vW %*% t(wW))
ematW <- senmatW * matW/lambdaW

lambdaW                 # population growth rate
senmatW                 # sensitivity matrix
ematW                   # elasticity matrix

# stable age distribution
sageW <- as.numeric(wW/sum(wW))       
op <- par(mfrow=c(1,2))
barplot(sageW, names.arg = c(0:18, "ad"), las=1, ylab = "Proportion in population", xlab = "Age class")    
barplot(c(sum(sageW[1]), sum(sageW[2:maxage])), names.arg = c("juv", "ad"), las=1, ylab = "Proportion in population", xlab = "Age class")    
par(op)

# Stage-specific population size in first year:
SSPopW <- round((sageW/sum(sageW[1:maxage]))*WWTP[1],0)
#sageW1 <- (ppm$mean$NW[,T]/sum(ppm$mean$NW[,T]))
#SSPopW <- round((sageW1/sum(sageW1[1:maxage]))*WWTP[1],0)


# ----------------------------------------
#  Full JAGS model allowing transitions:
# ----------------------------------------

sink("modSwans.jags")
cat("
    model { 
    ##########
    # Priors #
    ##########
    # Reproduction (offspring per pair)
    R ~ dgamma(mp.shape,mp.rate) I(minp, maxp)
    
    # Survival rates - Non-WWT population: 
    sa ~ dbeta(alpha.sa, beta.sa)
    si ~ dbeta(alpha.si, beta.si)
    sj ~ dbeta(alpha.sj, beta.sj)
    
    # Survival rates - WWT population:
    saW ~ dbeta(alpha.saW, beta.saW)
    siW ~ dbeta(alpha.siW, beta.siW)
    sjW ~ dbeta(alpha.sjW, beta.sjW)
    
    # Transition rates - Non-WWT to WWT:
    aTNW ~ dbeta(alpha.aTNW, beta.aTNW)
    iTNW ~ dbeta(alpha.iTNW, beta.iTNW)
    jTNW ~ dbeta(alpha.jTNW, beta.jTNW)   
    
    # Transition rates - WWT to Non-WWT:
    aTWN ~ dbeta(alpha.aTWN, beta.aTWN)
    iTWN ~ dbeta(alpha.iTWN, beta.iTWN)
    jTWN ~ dbeta(alpha.jTWN, beta.jTWN)
    
    # Priors for age-specific probability of breeding:
    for(j in 1:MAR){
    
    # non-WWT population:
    P[j] ~ dbeta(alpha.P[j], beta.P[j])
    PIn[j] ~ dbeta(alpha.PIn[j], beta.PIn[j])
    
    # WWT population:
    PW[j] ~ dbeta(alpha.PW[j], beta.PW[j])
    PWIn[j] ~ dbeta(alpha.PWIn[j], beta.PWIn[j])
    }
    
    ###########################
    # Model for initial state #
    ###########################
    for(i in 1:MAR){
    N[i,1] <- SSPop[i]
    NW[i,1] <- SSPopW[i]
    }#i
    
    # Loop over time
    for (t in 1:(T-1)){
    
    ####################
    # Population model #
    #################### 
    
    # Non-WWT population:
    N[1,t+1] ~ dpois(
    (N[2,t]) * (si * (1-iTNW) * R/2 * P[2])+
    (NW[2,t]) * (siW * iTWN * R/2 * PIn[2])+
    (N[3,t]) * (sa * (1-aTNW) * R/2 * P[3])+
    (NW[3,t]) * (saW * aTWN * R/2 * PIn[3])+
    (N[4,t]) * (sa * (1-aTNW) * R/2 * P[4])+
    (NW[4,t]) * (saW * aTWN * R/2 * PIn[4])+
    (N[5,t]) * (sa * (1-aTNW) * R/2 * P[5])+
    (NW[5,t]) * (saW * aTWN * R/2 * PIn[5])+
    (N[6,t]) * (sa * (1-aTNW) * R/2 * P[6])+
    (NW[6,t]) * (saW * aTWN * R/2 * PIn[6])+
    (N[7,t]) * (sa * (1-aTNW) * R/2 * P[7])+
    (NW[7,t]) * (saW * aTWN * R/2 * PIn[7])+
    (N[8,t]) * (sa * (1-aTNW) * R/2 * P[8])+
    (NW[8,t]) * (saW * aTWN * R/2 * PIn[8])+
    (N[9,t]) * (sa * (1-aTNW) * R/2 * P[9])+
    (NW[9,t]) * (saW * aTWN * R/2 * PIn[9])+
    (N[10,t]) * (sa * (1-aTNW) * R/2 * P[10])+
    (NW[10,t]) * (saW * aTWN * R/2 * PIn[10])+
    (N[11,t]) * (sa * (1-aTNW) * R/2 * P[11])+
    (NW[11,t]) * (saW * aTWN * R/2 * PIn[11])+
    (N[12,t]) * (sa * (1-aTNW) * R/2 * P[12])+
    (NW[12,t]) * (saW * aTWN * R/2 * PIn[12])+
    (N[13,t]) * (sa * (1-aTNW) * R/2 * P[13])+
    (NW[13,t]) * (saW * aTWN * R/2 * PIn[13])+
    (N[14,t]) * (sa * (1-aTNW) * R/2 * P[14])+
    (NW[14,t]) * (saW * aTWN * R/2 * PIn[14])+
    (N[15,t]) * (sa * (1-aTNW) * R/2 * P[15])+
    (NW[15,t]) * (saW * aTWN * R/2 * PIn[15])+
    (N[16,t]) * (sa * (1-aTNW) * R/2 * P[16])+
    (NW[16,t]) * (saW * aTWN * R/2 * PIn[16])+
    (N[17,t]) * (sa * (1-aTNW) * R/2 * P[17])+
    (NW[17,t]) * (saW * aTWN * R/2 * PIn[17])+
    (N[18,t]) * (sa * (1-aTNW) * R/2 * P[18])+
    (NW[18,t]) * (saW * aTWN * R/2 * PIn[18])+
    (N[19,t]) * (sa * (1-aTNW) * R/2 * P[19])+
    (NW[19,t]) * (saW * aTWN * R/2 * PIn[19])+
    (N[20,t]) * (sa * (1-aTNW) * R/2 * P[20])+
    (NW[20,t]) * (saW * aTWN * R/2 * PIn[20]))

    N[2,t+1] ~ dpois((sj * N[1,t] * (1-jTNW)) + (NW[1,t] * jTWN * sjW))
    N[3,t+1] ~ dpois((si * N[2,t] * (1-iTNW)) + (NW[2,t] * iTWN * siW))
    N[4,t+1] ~ dpois((sa * N[3,t] * (1-aTNW)) + (NW[3,t] * aTWN * saW))
    N[5,t+1] ~ dpois((sa * N[4,t] * (1-aTNW)) + (NW[4,t] * aTWN * saW))
    N[6,t+1] ~ dpois((sa * N[5,t] * (1-aTNW)) + (NW[5,t] * aTWN * saW))
    N[7,t+1] ~ dpois((sa * N[6,t] * (1-aTNW)) + (NW[6,t] * aTWN * saW))
    N[8,t+1] ~ dpois((sa * N[7,t] * (1-aTNW)) + (NW[7,t] * aTWN * saW))
    N[9,t+1] ~ dpois((sa * N[8,t] * (1-aTNW)) + (NW[8,t] * aTWN * saW))
    N[10,t+1] ~ dpois((sa * N[9,t] * (1-aTNW)) + (NW[9,t] * aTWN * saW))
    N[11,t+1] ~ dpois((sa * N[10,t] * (1-aTNW)) + (NW[10,t] * aTWN * saW)) 
    N[12,t+1] ~ dpois((sa * N[11,t] * (1-aTNW)) + (NW[11,t] * aTWN * saW))
    N[13,t+1] ~ dpois((sa * N[12,t] * (1-aTNW)) + (NW[12,t] * aTWN * saW))
    N[14,t+1] ~ dpois((sa * N[13,t] * (1-aTNW)) + (NW[13,t] * aTWN * saW))
    N[15,t+1] ~ dpois((sa * N[14,t] * (1-aTNW)) + (NW[14,t] * aTWN * saW))
    N[16,t+1] ~ dpois((sa * N[15,t] * (1-aTNW)) + (NW[15,t] * aTWN * saW))
    N[17,t+1] ~ dpois((sa * N[16,t] * (1-aTNW)) + (NW[16,t] * aTWN * saW))
    N[18,t+1] ~ dpois((sa * N[17,t] * (1-aTNW)) + (NW[17,t] * aTWN * saW))
    N[19,t+1] ~ dpois((sa * N[18,t] * (1-aTNW)) + (NW[18,t] * aTWN * saW))
    N[20,t+1] ~ dpois((sa * N[19,t] * (1-aTNW)) + (NW[19,t] * aTWN * saW)+(sa * N[20,t] * (1-aTNW)) + (NW[20,t] * aTWN * saW))
    
    # WWT population:
    NW[1,t+1] ~ dpois(
    (NW[2,t]) * (siW * (1-iTWN) * R/2 * PW[2])+
    (N[2,t]) * (si * iTNW * R/2 * PWIn[2])+
    (NW[3,t]) * (saW * (1-aTWN) * R/2 * PW[3])+
    (N[3,t]) * (sa * aTNW * R/2 * PWIn[3])+
    (NW[4,t]) * (saW * (1-aTWN) * R/2 * PW[4])+
    (N[4,t]) * (sa * aTNW * R/2 * PWIn[4])+
    (NW[5,t]) * (saW * (1-aTWN) * R/2 * PW[5])+
    (N[5,t]) * (sa * aTNW * R/2 * PWIn[5])+
    (NW[6,t]) * (saW * (1-aTWN) * R/2 * PW[6])+
    (N[6,t]) * (sa * aTNW * R/2 * PWIn[6])+
    (NW[7,t]) * (saW * (1-aTWN) * R/2 * PW[7])+
    (N[7,t]) * (sa * aTNW * R/2 * PWIn[7])+
    (NW[8,t]) * (saW * (1-aTWN) * R/2 * PW[8])+
    (N[8,t]) * (sa * aTNW * R/2 * PWIn[8])+
    (NW[9,t]) * (saW * (1-aTWN) * R/2 * PW[9])+
    (N[9,t]) * (sa * aTNW * R/2 * PWIn[9])+
    (NW[10,t]) * (saW * (1-aTWN) * R/2 * PW[10])+
    (N[10,t]) * (sa * aTNW * R/2 * PWIn[10])+
    (NW[11,t]) * (saW * (1-aTWN) * R/2 * PW[11])+
    (N[11,t]) * (sa * aTNW * R/2 * PWIn[11])+
    (NW[12,t]) * (saW * (1-aTWN) * R/2 * PW[12])+
    (N[12,t]) * (sa * aTNW * R/2 * PWIn[12])+
    (NW[13,t]) * (saW * (1-aTWN) * R/2 * PW[13])+
    (N[13,t]) * (sa * aTNW * R/2 * PWIn[13])+
    (NW[14,t]) * (saW * (1-aTWN) * R/2 * PW[14])+
    (N[14,t]) * (sa * aTNW * R/2 * PWIn[14])+
    (NW[15,t]) * (saW * (1-aTWN) * R/2 * PW[15])+
    (N[15,t]) * (sa * aTNW * R/2 * PWIn[15])+
    (NW[16,t]) * (saW * (1-aTWN) * R/2 * PW[16])+
    (N[16,t]) * (sa * aTNW * R/2 * PWIn[16])+
    (NW[17,t]) * (saW * (1-aTWN) * R/2 * PW[17])+
    (N[17,t]) * (sa * aTNW * R/2 * PWIn[17])+
    (NW[18,t]) * (saW * (1-aTWN) * R/2 * PW[18])+
    (N[18,t]) * (sa * aTNW * R/2 * PWIn[18])+
    (NW[19,t]) * (saW * (1-aTWN) * R/2 * PW[19])+
    (N[19,t]) * (sa * aTNW * R/2 * PWIn[19])+
    (NW[20,t]) * (saW * (1-aTWN) * R/2 * PW[20])+
    (N[20,t]) * (sa * aTNW * R/2 * PWIn[20]))

    NW[2,t+1] ~ dpois((sjW * NW[1,t] * (1-jTWN)) + (N[1,t] * jTNW * sj))
    NW[3,t+1] ~ dpois((siW * NW[2,t] * (1-iTWN)) + (N[2,t] * iTNW * si))
    NW[4,t+1] ~ dpois((saW * NW[3,t] * (1-aTWN)) + (N[3,t] * aTNW * sa))
    NW[5,t+1] ~ dpois((saW * NW[4,t] * (1-aTWN)) + (N[4,t] * aTNW * sa))
    NW[6,t+1] ~ dpois((saW * NW[5,t] * (1-aTWN)) + (N[5,t] * aTNW * sa))
    NW[7,t+1] ~ dpois((saW * NW[6,t] * (1-aTWN)) + (N[6,t] * aTNW * sa))
    NW[8,t+1] ~ dpois((saW * NW[7,t] * (1-aTWN)) + (N[7,t] * aTNW * sa))
    NW[9,t+1] ~ dpois((saW * NW[8,t] * (1-aTWN)) + (N[8,t] * aTNW * sa))
    NW[10,t+1] ~ dpois((saW * NW[9,t] * (1-aTWN)) + (N[9,t] * aTNW * sa))
    NW[11,t+1] ~ dpois((saW * NW[10,t] * (1-aTWN)) + (N[10,t] * aTNW * sa))
    NW[12,t+1] ~ dpois((saW * NW[11,t] * (1-aTWN)) + (N[11,t] * aTNW * sa))
    NW[13,t+1] ~ dpois((saW * NW[12,t] * (1-aTWN)) + (N[12,t] * aTNW * sa))
    NW[14,t+1] ~ dpois((saW * NW[13,t] * (1-aTWN)) + (N[13,t] * aTNW * sa))
    NW[15,t+1] ~ dpois((saW * NW[14,t] * (1-aTWN)) + (N[14,t] * aTNW * sa))
    NW[16,t+1] ~ dpois((saW * NW[15,t] * (1-aTWN)) + (N[15,t] * aTNW * sa))
    NW[17,t+1] ~ dpois((saW * NW[16,t] * (1-aTWN)) + (N[16,t] * aTNW * sa))
    NW[18,t+1] ~ dpois((saW * NW[17,t] * (1-aTWN)) + (N[17,t] * aTNW * sa))
    NW[19,t+1] ~ dpois((saW * NW[18,t] * (1-aTWN)) + (N[18,t] * aTNW * sa))
    NW[20,t+1] ~ dpois(((saW * NW[19,t] * (1-aTWN)) + (N[19,t] * aTNW * sa)) + ((saW * NW[20,t] * (1-aTWN)) + (N[20,t] * aTNW * sa)))
    
    ########################################
    # Calculation of population quantities #
    ######################################## 
    # Annual growth rate
    pop.change[t] <- sum(N[1:MAR,t+1])/sum(N[1:MAR,t])
    pop.changeW[t] <- sum(NW[1:MAR,t+1])/sum(NW[1:MAR,t])

    pop.numbers[t] <- sum(N[1:MAR,t])
    pop.numbersW[t] <- sum(NW[1:MAR,t])

    } #T

    lambda <- pop.change[T-1]
    mean.lambda <- mean(pop.change) #mean(pop.change[u:T-1])
    lambdaW <- pop.changeW[T-1]
    mean.lambdaW <- mean(pop.changeW) #mean(pop.changeW[u:T-1])
    
    } # model
    ",fill = TRUE)
sink()

# Bundle data
T <- length(Yr)
u <- 2   # Burn-in for the calculation of the stochastic population growth rate (T > u)

jd <- list(T = T, u=u,
           # Fecundity
           mp.shape = mp.shape, mp.rate = mp.rate, # mean number of observed offspring
           maxp = max(Prod), # max number of observed offspring
           minp = min(Prod), # max number of observed offspring
           # Starting population at stable-stage distributions
           SSPopW = SSPopW, SSPop = SSPop, MAR = maxage, 
           # Probability of breeding
           alpha.P = beta.parameters(P.m, P.sd)$alpha, beta.P = beta.parameters(P.m, P.sd)$beta,
           alpha.PW = beta.parameters(PW.m, PW.sd)$alpha, beta.PW = beta.parameters(PW.m, PW.sd)$beta,
           alpha.PIn = beta.parameters(PIn.m, PIn.sd)$alpha, beta.PIn = beta.parameters(PIn.m, PIn.sd)$beta,
           alpha.PWIn = beta.parameters(PWIn.m, PWIn.sd)$alpha, beta.PWIn = beta.parameters(PWIn.m, PWIn.sd)$beta,
           # survival rates
           alpha.sa = beta.parameters(mean(SadNN), sd(SadNN))$alpha, beta.sa = beta.parameters(mean(SadNN), sd(SadNN))$beta,
           alpha.si = beta.parameters(mean(S1yNN), sd(S1yNN))$alpha, beta.si = beta.parameters(mean(S1yNN), sd(S1yNN))$beta,
           alpha.sj = beta.parameters(mean(SjuvNN), sd(SjuvNN))$alpha, beta.sj = beta.parameters(mean(SjuvNN), sd(SjuvNN))$beta,
           alpha.saW = beta.parameters(mean(SadWWT), sd(SadWWT))$alpha, beta.saW = beta.parameters(mean(SadWWT), sd(SadWWT))$beta,
           alpha.siW = beta.parameters(mean(S1yWWT), sd(S1yWWT))$alpha, beta.siW = beta.parameters(mean(S1yWWT), sd(S1yWWT))$beta,
           alpha.sjW = beta.parameters(mean(SjuvWWT), sd(SjuvWWT))$alpha, beta.sjW = beta.parameters(mean(SjuvWWT), sd(SjuvWWT))$beta,
           # Transition rates
           alpha.aTNW = beta.parameters(mean(TadNN_adWWT), sd(TadNN_adWWT))$alpha, beta.aTNW = beta.parameters(mean(TadNN_adWWT), sd(TadNN_adWWT))$beta,
           alpha.iTNW = beta.parameters(mean(T1yNN_adWWT), sd(T1yNN_adWWT))$alpha, beta.iTNW = beta.parameters(mean(T1yNN_adWWT), sd(T1yNN_adWWT))$beta,
           alpha.jTNW = beta.parameters(mean(TjuvNN_1yWWT), sd(TjuvNN_1yWWT))$alpha, beta.jTNW = beta.parameters(mean(TjuvNN_1yWWT), sd(TjuvNN_1yWWT))$beta,
           alpha.aTWN = beta.parameters(mean(TadWWT_adNN), sd(TadWWT_adNN))$alpha, beta.aTWN = beta.parameters(mean(TadWWT_adNN), sd(TadWWT_adNN))$beta,
           alpha.iTWN = beta.parameters(mean(T1yWWT_adNN), sd(T1yWWT_adNN))$alpha, beta.iTWN = beta.parameters(mean(T1yWWT_adNN), sd(T1yWWT_adNN))$beta,
           alpha.jTWN = beta.parameters(mean(TjuvWWT_1yNN), sd(TjuvWWT_1yNN))$alpha, beta.jTWN = beta.parameters(mean(TjuvWWT_1yNN), sd(TjuvWWT_1yNN))$beta)


# Initial values
Ninit <- matrix(5000, nrow = maxage, ncol = T)
Ninit[,1] <- NA
NinitW<- matrix(50, nrow = maxage, ncol = T)
NinitW[,1] <- NA
inits <- function(){list(N = Ninit, NW = NinitW)} 

# Parameters monitored
parameters <- c("lambda","mean.lambda","lambdaW","mean.lambdaW",
                "pop.change","pop.changeW","pop.numbers","pop.numbersW", "N", "NW") #"N",

# MCMC settings
ni <- 500#00
nt <- 1#0
nb <- 100#00
nc <- 3

# Call JAGS from R
ppm <- jags(jd, inits, parameters, "modSwans.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC = F, parallel=T)

print(ppm, 5)


####### Model with non-WWT only (no transitions) ######

# JAGS model:
sink("modSwansNN.jags")
cat("
    model { 
    ##########
    # Priors #
    ##########
    # Reproduction (offspring per pair)
    R ~ dgamma(mp.shape,mp.rate) I(minp, maxp)
    
    # Survival rates - Non-WWT population: 
    sa ~ dbeta(alpha.sa, beta.sa)
    si ~ dbeta(alpha.si, beta.si)
    sj ~ dbeta(alpha.sj, beta.sj)
    
    
    # Priors for age-specific probability of breeding:
    for(j in 1:MAR){
    
    # non-WWT population:
    P[j] ~ dbeta(alpha.P[j], beta.P[j])
    
    }
    
    ###########################
    # Model for initial state #
    ###########################
    for(i in 1:MAR){
    N[i,1] <- SSPop[i]
    }#i
    
    # Loop over time
    for (t in 1:(T-1)){
    
    ####################
    # Population model #
    #################### 
    
    # Non-WWT population:
    N[1,t+1] ~ dpois(
    (N[2,t]) * (si * R/2 * P[2])+
    (N[3,t]) * (sa * R/2 * P[3])+
    (N[4,t]) * (sa * R/2 * P[4])+
    (N[5,t]) * (sa * R/2 * P[5])+
    (N[6,t]) * (sa * R/2 * P[6])+
    (N[7,t]) * (sa * R/2 * P[7])+
    (N[8,t]) * (sa * R/2 * P[8])+
    (N[9,t]) * (sa * R/2 * P[9])+
    (N[10,t]) * (sa * R/2 * P[10])+
    (N[11,t]) * (sa * R/2 * P[11])+
    (N[12,t]) * (sa * R/2 * P[12])+
    (N[13,t]) * (sa * R/2 * P[13])+
    (N[14,t]) * (sa * R/2 * P[14])+
    (N[15,t]) * (sa * R/2 * P[15])+
    (N[16,t]) * (sa * R/2 * P[16])+
    (N[17,t]) * (sa * R/2 * P[17])+
    (N[18,t]) * (sa * R/2 * P[18])+
    (N[19,t]) * (sa * R/2 * P[19])+
    (N[20,t]) * (sa * R/2 * P[20]))
    
    N[2,t+1] ~ dbin(sj, N[1,t])
    N[3,t+1] ~ dbin(si, N[2,t])
    N[4,t+1] ~ dbin(sa, N[3,t])
    N[5,t+1] ~ dbin(sa, N[4,t])
    N[6,t+1] ~ dbin(sa, N[5,t])
    N[7,t+1] ~ dbin(sa, N[6,t])
    N[8,t+1] ~ dbin(sa, N[7,t])
    N[9,t+1] ~ dbin(sa, N[8,t])
    N[10,t+1] ~ dbin(sa, N[9,t])
    N[11,t+1] ~ dbin(sa, N[10,t]) 
    N[12,t+1] ~ dbin(sa, N[11,t])
    N[13,t+1] ~ dbin(sa, N[12,t])
    N[14,t+1] ~ dbin(sa, N[13,t])
    N[15,t+1] ~ dbin(sa, N[14,t])
    N[16,t+1] ~ dbin(sa, N[15,t])
    N[17,t+1] ~ dbin(sa, N[16,t])
    N[18,t+1] ~ dbin(sa, N[17,t])
    N[19,t+1] ~ dbin(sa, N[18,t])
    N[20,t+1] ~ dbin(sa, (N[19,t]+N[20,t]))
    
    #######################################
    # Calculation of population quantities #
    ####################################### 
    # Annual growth rate
    pop.change[t] <- sum(N[1:MAR,t+1])/sum(N[1:MAR,t])
    pop.numbers[t] <- sum(N[1:MAR,t])
    
    } #T
    lambda <- pop.change[T-1]
    mean.lambda <- mean(pop.change)
    
    } # model
    ",fill = TRUE)
sink()

jd.nt <- list(T = T, u=u,
           # Fecundity
           mp.shape = mp.shape, mp.rate = mp.rate, # mean number of observed offspring
           maxp = max(Prod), # max number of observed offspring
           minp = min(Prod), # max number of observed offspring
           
           # Starting population at stable-stage distributions
           SSPop = SSPop, MAR = maxage, 
           
           # Probability of breeding
           alpha.P = beta.parameters(P.m, P.sd)$alpha, beta.P = beta.parameters(P.m, P.sd)$beta,
          
           # survival rates
           alpha.sa = beta.parameters(mean(SadNN), sd(SadNN))$alpha, beta.sa = beta.parameters(mean(SadNN), sd(SadNN))$beta,
           alpha.si = beta.parameters(mean(S1yNN), sd(S1yNN))$alpha, beta.si = beta.parameters(mean(S1yNN), sd(S1yNN))$beta,
           alpha.sj = beta.parameters(mean(SjuvNN), sd(SjuvNN))$alpha, beta.sj = beta.parameters(mean(SjuvNN), sd(SjuvNN))$beta)

# Initial values
Ninit <- matrix(5000, nrow = maxage, ncol = T)
Ninit[,1] <- NA
inits.nt <- function(){list(N = Ninit)} 

# Parameters monitored
parameters.nt <- c("lambda","mean.lambda",
                "pop.change", "N", "pop.numbers")  #"N",


# Call JAGS from R
ppm.nt <- jags(jd.nt, inits.nt, parameters.nt, parallel = T, "modSwansNN.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC = F)

print(ppm.nt, 5)

# 
tempN<-tempW<-tempN.nt<-0
for(i in 1:length(ppm$mean$pop.numbers)-1){
  tempN[i]<-ppm$mean$pop.numbers[i+1]/ppm$mean$pop.numbers[i]
  tempW[i]<-ppm$mean$pop.numbersW[i+1]/ppm$mean$pop.numbersW[i]
  tempN.nt[i]<-ppm.nt$mean$pop.numbers[i+1]/ppm.nt$mean$pop.numbers[i]
}
mean(tempN)
mean(tempW)
mean(tempN.nt)

### Calulated population trajectories and 95% B CI: 
Years <- seq(from=1,to=length(Yr),by=1)
pairsN <-pairsN.nt <-pairsW <- pairsW <-pairsWL<-pairsWU <- pairsLMM <- matrix(0, nrow = max(Years), ncol = 1, byrow = T)
pairsN[1] <- pairsN.nt[1] <- pairsLMM[1] <- sum(SSPop[1:maxage])
pairsW[1] <- sum(SSPopW[1:maxage])

for (i in 2:max(Years)){
  #pairsN.nt[i] <- sum(SSPop[1:maxage]) * ppm.nt$mean$mean.lambda^(i-1)
  pairsN[i] <- sum(SSPop[1:maxage]) * SSM.l.N^(i-1)
  pairsW[i] <- sum(SSPopW[1:maxage]) * SSM.l.W^(i-1)
  pairsLMM[i] <- sum(SSPop[1:maxage]) * lambda^(i-1)
}

modpairsN<-modpairsNL<-modpairsNU <- modpairsN.nt<-modpairsNL.nt<-modpairsNU.nt <- modpairsW <-modpairsWL<-modpairsWU <- matrix(0, nrow = max(Years), ncol = 1, byrow = T)
modpairsN[1] <- modpairsNL[1] <- modpairsNU[1] <- sum(SSPop[1:maxage])
modpairsN.nt[1] <- modpairsNL.nt[1] <- modpairsNU.nt[1] <- sum(SSPop[1:maxage])
modpairsW[1] <- modpairsWL[1] <- modpairsWU[1] <- sum(SSPopW[1:maxage])
for (i in 2:max(Years)){
  modpairsW[i] <- sum(ppm$mean$NW[1:maxage,i])
  modpairsWU[i] <- sum(ppm$q97.5$NW[1:maxage,i])
  modpairsWL[i] <- sum(ppm$q2.5$NW[1:maxage,i])
  
  modpairsN[i] <- sum(ppm$mean$N[1:maxage,i])
  modpairsNU[i] <- sum(ppm$q97.5$N[1:maxage,i])
  modpairsNL[i] <- sum(ppm$q2.5$N[1:maxage,i])
  
  modpairsN.nt[i] <- sum(ppm.nt$mean$N[1:maxage,i])
  modpairsNU.nt[i] <- sum(ppm.nt$q97.5$N[1:maxage,i])
  modpairsNL.nt[i] <- sum(ppm.nt$q2.5$N[1:maxage,i])
}


######## Lambdas #######

round(lambda,3) # LMM deterministic Non-WWT
round(ppm$mean$mean.lambda,3) # JAGS mean stochastic non-WWT 

round(lambdaW,3) # LMM deterministic WWT
round(ppm$mean$mean.lambdaW,3) # JAGS mean stochastic WWT 


######## Diagnostic plots #######

# quartz(height=6,width=8)
# par(mfrow=c(2,5),mar=c(4,4,0.2,0.2))
# plot(density(ppm$sims.list$sa),main='')
# lines(density(SadNN),col='red')
# 
# plot(density(ppm$sims.list$si),main='')
# lines(density(S1yNN),col='red')
# 
# plot(density(ppm$sims.list$sj),main='')
# lines(density(SjuvNN),col='red')
# 
# plot(density(ppm$sims.list$saW),main='')
# lines(density(SadWWT),col='red')
# 
# plot(density(ppm$sims.list$siW),main='')
# lines(density(S1yWWT),col='red')
# 
# plot(density(ppm$sims.list$sjW),main='')
# lines(density(SjuvWWT),col='red')
# 
# plot(density(ppm$sims.list$R),main='')
# lines(density(Prod),col='red')
# 
# plot(density(ppm$sims.list$jTNW),main='')
# lines(density(TjuvNN_1yWWT),col='red')
# 
# plot(density(ppm$sims.list$jTWN),main='')
# lines(density(TjuvWWT_1yNN),col='red')
# 
# plot(density(ppm$sims.list$iTNW),main='')
# lines(density(T1yNN_adWWT),col='red')
# 
# plot(density(ppm$sims.list$iTWN),main='')
# lines(density(T1yWWT_adNN),col='red')

# par(mfrow=c(3,5),mar=c(4,4,0.2,0.2))
# for(i in 1:15){
#   plot(density(ppm$sims.list$PW[,i]),main='')
#   lines(density(PWWTBoth[,i]),col='red')
# }
# 
# for(i in 1:15){
#   plot(density(ppm$sims.list$P[,i]),main='')
#   lines(density(PNNBoth[,i]),col='red')
# }
