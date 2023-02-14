### Whooper Swan Breeding Analysis
### Analysis done August 2017
### Write up started Jan 2018
### Further analysis Jan 2018
### Write up May 2020
### R.Inger

library(MASS)
library(lme4)
library(ez)
library(MuMIn)
library(arm)
library(chron)
library(lmerTest)
library(effects)
library(boot)
library(ggplot2)
library(beepr)
rm(list=ls())
setwd("")#Set this is working directory

data<- read.table("StandardisedWhooperData.csv",sep=",",header=TRUE)
attach(data)

#STANDARDISE
Age.Z<-scale(Age,center=TRUE,scale=TRUE)
PBS.Z<-scale(PBS,center=TRUE,scale=TRUE)
PopBSAll.Z<-scale(PopBSAll,center=TRUE,scale=TRUE)
PrevSeaMaxDUr.Z<-scale(PrevSeaMaxDUr,center=TRUE,scale=TRUE)
PrevSSites.Z<-scale(PrevSSites,center=TRUE,scale=TRUE)
PrevSeaJulian.Z<-scale(PrevSeaJulian,center=TRUE,scale=TRUE)
dataII<-cbind(data,Age.Z,PBS.Z,PopBSAll.Z,PrevSeaMaxDUr,PrevSSites.Z,PrevSeaJulian.Z)

#write.table(data,"StandardisedWhooperData.csv",sep = ",",col.names = TRUE)

#Global model
global<-glmer(CygnetsIII~
                Age.Z+                     # Age
                I(Age.Z*Age.Z)+            # Age squared
                
                WWTNEW+                    #WWT Site Use Previous 2 winters
                PBS.Z+                     # Previous Breeding success
                PopBSAll.Z+                # Population breeding success
                PrevSSites+                # Number of sites in previous season
                PrevSeaMaxDUr+             # Previous season max duration
                
                WWTNEW*Age.Z+              # WWT * Age
                WWTNEW*I(Age.Z*Age.Z)+ 
                PBS.Z*I(Age.Z*Age.Z)+               # PBS * Age
                PopBSAll.Z*Age.Z+
                PrevSSites.Z+Age.Z+        # Number of sites in previous season * Age
                PrevSeaMaxDUr.Z*Age.Z+     # Previous season max duration * Age
                PBS.Z*PopBSAll.Z+
                
                (1|Season)+(1|SWID)+(1|Ind),data=data,family=binomial,control=glmerControl(optimizer="bobyqa")) #Random Season + Individual + Observation level
#fixed-effect model matrix is rank deficient so dropping 2 columns / coefficients
#Warning messages:
# 1: In optwrap(optimizer, devfun, start, rho$lower, control = control,  :
#                 convergence code 1 from bobyqa: bobyqa -- maximum number of function evaluations exceeded
# 2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                                 Model failed to converge with max|grad| = 0.130438 (tol = 0.001, component 1)
# 3: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                                                 Model is nearly unidentifiable: very large eigenvalue
#                                               - Rescale variables?
#SUMMARY
summary(global)
r.squaredGLMM(global) #R2m = 0.422 R2c = 0.552 theoretical
relgrad <- with(global@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) #0.037

##########################################################
## Global model - WWTNEW*Age.Z - PBS & AGE Interactions  #
##########################################################
#PBS & Age interactions removed as one (PBS) is partially derived from the other (Age)

global2<-glmer(CygnetsIII~
                Age.Z+                     # Age
                I(Age.Z*Age.Z)+            # Age squared
                
                WWTNEW+                    #WWT Site Use Previous 2 winters
                PBS.Z+                     # Previous Breeding success
                PopBSAll.Z+                # Population breeding success
                PrevSSites+                # Number of sites in previous season
                PrevSeaMaxDUr+             # Previous season max duration
                
                #WWTNEW*Age.Z+              # WWT * Age
                WWTNEW*I(Age.Z*Age.Z)+ 
                #PBS.Z*I(Age.Z*Age.Z)+               # PBS * Age
                #PBS*Age.Z+
                PopBSAll.Z*Age.Z+
                PrevSSites.Z+Age.Z+        # Number of sites in previous season * Age
                PrevSeaMaxDUr.Z*Age.Z+     # Previous season max duration * Age
                PBS.Z*PopBSAll.Z+
                
                (1|Season)+(1|SWID)+(1|Ind),data=data,family=binomial,control=glmerControl(optimizer="bobyqa")) #Random Season + Individual + Observation level
# fixed-effect model matrix is rank deficient so dropping 2 columns / coefficients
# Warning messages:
#   1: In optwrap(optimizer, devfun, start, rho$lower, control = control,  :
#                   convergence code 1 from bobyqa: bobyqa -- maximum number of function evaluations exceeded
#                 2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                                   Model failed to converge with max|grad| = 0.257169 (tol = 0.001, component 1)
#                                 3: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                                                   Model is nearly unidentifiable: very large eigenvalue
#                                                 - Rescale variables?
#SUMMARY
summary(global2)
r.squaredGLMM(global2) #R2m = 0.4229 R2c = 0.5445
relgrad <- with(global2@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) #0.05

##########################################################
## Global model - WWTNEW*Age.Z - PBS & AGE Interactions
# - PrevSSites & Interactions
# - This stopped one of the data deficient warnings
# - Also main effect was highest p value
##########################################################


global3<-glmer(CygnetsIII~
                Age.Z+                     # Age
                I(Age.Z*Age.Z)+            # Age squared
                
                WWTNEW+                    #WWT Site Use Previous 2 winters
                PBS.Z+                     # Previous Breeding success
                PopBSAll.Z+                # Population breeding success
                #PrevSSites+                # Number of sites in previous season
                PrevSeaMaxDUr+             # Previous season max duration
                
                #WWTNEW*Age.Z+              # WWT * Age
                WWTNEW*I(Age.Z*Age.Z)+ 
                #PBS.Z*I(Age.Z*Age.Z)+               # PBS * Age
                #PBS*Age.Z+
                PopBSAll.Z*Age.Z+
                #PrevSSites.Z+Age.Z+        # Number of sites in previous season * Age
                PrevSeaMaxDUr.Z*Age.Z+     # Previous season max duration * Age
                PBS.Z*PopBSAll.Z+
                
                (1|Season)+(1|SWID)+(1|Ind),data=data,family=binomial,control=glmerControl(optimizer="bobyqa")) #Random Season + Individual + Observation level

# fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
# Warning messages:
#   1: In optwrap(optimizer, devfun, start, rho$lower, control = control,  :
#                   convergence code 1 from bobyqa: bobyqa -- maximum number of function evaluations exceeded
#                 2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                                   Model failed to converge with max|grad| = 0.466559 (tol = 0.001, component 1)
#                                 3: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                                                   Model is nearly unidentifiable: very large eigenvalue
#                                                 - Rescale variables?
#SUMMARY
summary(global3)
r.squaredGLMM(global3) #R2m = 0.4234 R2c = 0.5439
relgrad <- with(global3@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) #0.04


##########################################################
## Global model - WWTNEW*Age.Z - PBS & AGE Interactions
# - PrevSSites & Interactions
# - This stopped one of the data deficient warnings
# - Also main effect was highest p value
# - PopBSAll & Interactions
##########################################################


global4<-glmer(CygnetsIII~
                Age.Z+                     # Age
                I(Age.Z*Age.Z)+            # Age squared
                
                WWTNEW+                    #WWT Site Use Previous 2 winters
                PBS.Z+                     # Previous Breeding success
                #PopBSAll.Z+                # Population breeding success
                #PrevSSites+                # Number of sites in previous season
                PrevSeaMaxDUr+             # Previous season max duration
                
                #WWTNEW*Age.Z+              # WWT * Age
                WWTNEW*I(Age.Z*Age.Z)+ 
                #PBS.Z*I(Age.Z*Age.Z)+               # PBS * Age
                #PBS*Age.Z+
                #PopBSAll.Z*Age.Z+
                #PrevSSites.Z+Age.Z+        # Number of sites in previous season * Age
                PrevSeaMaxDUr.Z*Age.Z+     # Previous season max duration * Age
                #PBS.Z*PopBSAll.Z+
                
                (1|Season)+(1|SWID)+(1|Ind),data=data,family=binomial,control=glmerControl(optimizer="bobyqa")) #Random Season + Individual + Observation level
# fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
# Warning messages:
#   1: In optwrap(optimizer, devfun, start, rho$lower, control = control,  :
#                   convergence code 1 from bobyqa: bobyqa -- maximum number of function evaluations exceeded
#                 2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                                   Model failed to converge with max|grad| = 0.0569004 (tol = 0.001, component 1)
#                                 3: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                                                   Model is nearly unidentifiable: very large eigenvalue
#                                                 - Rescale variables?

#SUMMARY
summary(global4)
r.squaredGLMM(global4) #R2m = 0.4095 R2c = 0.5442
relgrad <- with(global4@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) #0.02

##FINAL MODEL
##########################################################
## Global model - WWTNEW*Age.Z - PBS & AGE Interactions
# - PrevSSites & Interactions
# - This stopped one of the data deficient warnings
# - Also main effect was highest p value
# - PopBSAll & Interactions
# - PrevSeaMaxDUr & Interactions
# - Gets rid of data deficient warnings
# - Also gets rid of all on convergence warnings
##########################################################


global5<-glmer(CygnetsIII~
                Age.Z+                     # Age
                I(Age.Z*Age.Z)+            # Age squared
                
                WWTNEW+                    #WWT Site Use Previous 2 winters
                PBS.Z+                     # Previous Breeding success
                #PopBSAll.Z+                # Population breeding success
                #PrevSSites+                # Number of sites in previous season
                #PrevSeaMaxDUr+             # Previous season max duration
                
                #WWTNEW*Age.Z+              # WWT * Age
                WWTNEW*I(Age.Z*Age.Z)+ 
                #PBS.Z*I(Age.Z*Age.Z)+               # PBS * Age
                #PBS*Age.Z+
                #PopBSAll.Z*Age.Z+
                #PrevSSites.Z+Age.Z+        # Number of sites in previous season * Age
                #PrevSeaMaxDUr.Z*Age.Z+     # Previous season max duration * Age
                #PBS.Z*PopBSAll.Z+
                
                (1|Season)+(1|SWID)+(1|Ind),data=data,family=binomial,control=glmerControl(optimizer="bobyqa")) #Random Season + Individual + Observation level
beep(2)

#Table in supplementary materials
global5<-glmer(CygnetsIII~
                Age.Z+                     # Age
                I(Age.Z*Age.Z)+            # Age squared
                WWTNEW+                    # WWT Site Use Previous 2 winters
                PBS.Z+                     # Previous Breeding success
                PopBSAll                   # Population Level Breeding Success 
                WWTNEW*I(Age.Z*Age.Z)+ 
                (1|Season)+(1|SWID)+(1|Ind),data=data,family=binomial,control=glmerControl(optimizer="bobyqa")) #Random Season + Individual + Observation level
summary(global5)
r.squaredGLMM(global5) #R2m = 0.419 R2c = 0.539
relgrad <- with(global5@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # 1.171156e-05

# Gives singularity warning, therefore remove PopBSAll as lower significance

#This is the model used to produce the Bayesian Model
global6<-glmer(CygnetsIII~
                Age.Z+                     # Age
                I(Age.Z*Age.Z)+            # Age squared
                WWTNEW+                    # WWT Site Use Previous 2 winters
                PBS.Z+                     # Previous Breeding Success
                WWTNEW*I(Age.Z*Age.Z)+ 
                (1|Season)+(1|SWID)+(1|Ind),data=data,family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(global6)
r.squaredGLMM(global6) #R2m = 0.405 R2c = 0.5417
relgrad <- with(global6@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) #0.00005
#NO ERRORS!

#SUMMARY
summary(global5)
r.squaredGLMM(global5) #R2m = 0.404 R2c = 0.5417
relgrad <- with(global5@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) #0.00005

## Extract model fitted proportions and CIs ##
Fit<- effect(term="WWTNEW*I(Age.Z*Age.Z)", mod=global2, se = T, confidence.level=0.95, xlevels=30)
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
quartz()
par(mfrow=c(2,2))
#WWTBOTH
plot(FP[seq(1,length(FP),4)],type="n",xlab="Age",ylab="Breeding Probability",main="Breeding Probability WWT BOTH", ylim=c(0,0.4))
lines(FP[seq(1,length(FP),4)],type="l",lwd="3",col="blue") 
lines(UP[seq(1,length(FP),4)],type="l",lwd="1",col="blue") 
lines(LP[seq(1,length(FP),4)],type="l",lwd="1",col="blue") 
#WWTIN
plot(FP[seq(2,length(FP),4)],type="n",xlab="Age",ylab="Breeding Probability",main="Breeding Probability WWT IN", ylim=c(0,0.4))
lines(FP[seq(2,length(FP),4)],type="l",lwd="3",col="red") 
lines(UP[seq(2,length(FP),4)],type="l",lwd="1",col="red") 
lines(LP[seq(2,length(FP),4)],type="l",lwd="1",col="red") 
#WWTNON
plot(FP[seq(3,length(FP),4)],type="n",xlab="Age",ylab="Breeding Probability",main="Breeding Probability WWT NON", ylim=c(0,0.4))
lines(FP[seq(3,length(FP),4)],type="l",lwd="3",col="green") 
lines(UP[seq(3,length(FP),4)],type="l",lwd="1",col="green") 
lines(LP[seq(3,length(FP),4)],type="l",lwd="1",col="green") 
#WWTOUT
plot(FP[seq(4,length(FP),4)],type="n",xlab="Age",ylab="Breeding Probability",main="Breeding Probability WWT OUT", ylim=c(0,0.4))
lines(FP[seq(4,length(FP),4)],type="l",lwd="3",col="purple") 
lines(UP[seq(4,length(FP),4)],type="l",lwd="1",col="purple") 
lines(LP[seq(4,length(FP),4)],type="l",lwd="1",col="purple") 

