## Including capelin and cod biomass as covariates as covariates on fecundity 
## Martin Biuw and Tor Arne Øigård
## Modified from:
## Fitting state-space models to seal populations with scarce data
## Tor Arne Øigård and Hans J. Skaug
## ICES Journal of Marine Science, 2014
## Contact:  Martin Biuw (martin.biuw@hi.no) or Tor Arne Øigård (tor.arne.oigard@nr.no)

## Assume the Rstudio project opens in the SAMFH directory, so no change in wd needed
##setwd('C:/Users/a5406/Documents/Populasjonsmodellering/rSPAMM/rSPAMM-master')
source('./R/read.ICES.R')

library(TMB)

#####################################
#Some functions needed below
ilogit1 <- function(x){
  return(1/(1+exp(-x)))
}

logit <- function(x){
  return(log(x/(1-x)))
}

ilogit <- function(x){
  return(exp(x)/(1+exp(x)))
}

blogit <- function(x, m, s){
  2/(1+exp(-(x-m)/s))-1
}

######################################

compile("./R/harps_and_hoods_population_model4.cpp",
        "-O1 -g", DLLFLAGS="")
dyn.load(dynlib("./R/harps_and_hoods_population_model4"))

###########
#Data part
###########
#Read in data
catch_data <- read.table("./Data/harpeast/catch_data.dat",header = FALSE)				#Catch data
pup_production <- read.table("./Data/harpeast/pup_production.dat",header = FALSE)			#Pup production estimates
fecundity <- read.table("./Data/harpeast/fecundity.dat",header = FALSE)				#Available fecundity data
##Pmatrix <- read.table("./Data/harpeast/wgharp.pma",header = FALSE)				#Birth ogives
## Or with collapsed maturity curve:
Pmatrix <- read.table("./Data/harpeast/collapsed.ogi",header = TRUE)				#Birth ogives
Pmatrix <- c(Pmatrix[,2], rep(1, 4))
priors <- read.table("./Data/harpeast/priors.dat",header = FALSE)					#Priors used

## Add extra priors:
priors <- rbind(priors, data.frame(V1=c(0.85, 0.5, -0.5, 0.1, 0.1), 
                                   V2=c(0.075, 0.25, 0.25, 0.01, 0.01)))
years.of.prediction = 0											#Number of years to run projections

## Read and prepare fish data
cap <- read.ICES()
## Read capelin index from cod stomach contents 
## prior to 1972 (From Marshall et al. 2000)
histCap <- read.csv('./Data/histCap.csv')
histCap$method <- rep('Rec', nrow(histCap))
cap$Data <- plyr::rbind.fill(histCap[-nrow(histCap),], cap$Data)
cap$Data$method[which(is.na(cap$Data$method))] <- 'Meas'
cap$Data$method <- as.factor(cap$Data$method)
cod <- read.ICES('./Data/_9841_9841.xml')

fish <- cod$Data[,c(1,6,9)]
names(fish)[c(2,3)] <- c('CodTB', 'CodSSB')
fish <- merge(fish, cap$Data[,c(1:4)], by='Year', all.x=T, all.y=F)
names(fish)[c(4,6)] <- c('CapTB', 'CapSSB')
fish <- fish[which(!is.na(fish$CodTB) & !is.na(fish$CapTB)),]
fish$suit <- (fish$CapTB/max(fish$CapTB, na.rm=T))-(fish$CodTB/max(fish$CodTB, na.rm=T))
fish$suit <- fish$suit-min(fish$suit)
fish$suit <- fish$suit/max(fish$suit)

#Prepare input to the model
data <- list()
data$Amax = 20													#Maximum age group
data$Cdata = as.matrix(catch_data)								#Add catch data
data$Nc = nrow(data$Cdata)									#Number of years with catch data
data$pup_production = as.matrix(pup_production)					#Pup production estimates
data$Np = nrow(data$pup_production)								#Number of years with pup production data
data$Fecundity = as.matrix(fecundity)							#Observed fecundity rates
data$Nf = nrow(fecundity)									#Number of years with fecundity data
##data$Pmat = as.matrix(Pmatrix)									#Maturity curves
## Or with collapsed maturity curve:
data$Pmat = Pmatrix									      #Maturity curve
data$Npred = years.of.prediction									#Number of years to run projections
data$priors = as.matrix(priors)									#Priors for estimated parameters
data$Npriors = nrow(priors)									#Number of priors
data$CQuota = c(0,0)											#Catch level in future projections
data$cap = fish$CapTB              #Capelin biomass
data$cod = fish$CodTB              #Cod biomass


#Initial values
Kinit = 2000000								#Initial population size
Minit = 0.09									#Natural adult mortality
M0init = 0.27								  #Natural pup mortality
finit = 0.85									#Mean fecundity rate (or pregnancy, more accurately)
##binit = 0.7									  #Intercept parameter
bCapinit = 0.5								#Slope parameter for capelin
bCodinit = -0.5								#Slope parameter for cod
Procinit = 0.1                #Process error
Obsinit = 0.1                 #Observation error
Prinit = 0.5						      #Mixture proportion parameter

#Transform some parameters to ensure that estimates are > 0
parameters <- list()
parameters$logK= log(Kinit)  			
parameters$Mtilde= logit(Minit)				
parameters$M0tilde= logit(M0init)				
parameters$ftilde = logit(finit)
##parameters$logBeta = log(binit)
parameters$logBetaCap = sign(bCapinit) * log(abs(bCapinit))
parameters$logBetaCod = sign(bCapinit) * log(abs(bCapinit))
parameters$logSdProc = log(Procinit)
parameters$logSdObs = log(Obsinit)
##parameters$logitp = log(Prinit)
parameters$F=rep(0,length.out = (data$Nc+data$Npred))


obj <- MakeADFun(data,parameters,DLL="harps_and_hoods_population_model4")

obj$fn()
obj$gr()

system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))

rep<-sdreport(obj, getJointPrecision=TRUE)
rep.matrix <- summary(rep)
rep.rnames = rownames(rep.matrix)
indN0 = which(rep.rnames=="N0");indN0 <- indN0[-1]
indN1 = which(rep.rnames=="N1");indN1 <- indN1[-1]
indFt = which(rep.rnames=="Ft");indFt <- indFt[-1]
yrs = 1946+c(0:(length(indFt)-1))

#Extract parameters
Kest = exp(opt$par[1])          
Mest = ilogit(opt$par[2])       
M0est = ilogit(opt$par[3])      
fest = ilogit(opt$par[4])       
bCapest = exp(opt$par[5])             
bCodest = exp(opt$par[6])    
Procest = exp(opt$par[7])
Obsest = exp(opt$par[8])

#allsd<-sqrt(diag(solve(rep$jointPrecision)))
#plsd <- obj$env$parList(par=allsd)
##par(mfrow=c(2,1))
plot(yrs,rep.matrix[indN1,1],type = "l",xlab = "Year",ylab = "Abundance",lwd = 2,col = 3)
lines(yrs,rep.matrix[indN0,1],col = "blue",lwd = 2)
points(data$pup_production[,1],data$pup_production[,2], pch=21, bg = 2, cex = 1.2)

plot(yrs,rep.matrix[indFt,1],type = "l",xlab = "Year",ylab = "Fecundity",lwd = 2,col = 4)
points(data$Fecundity[,1], data$Fecundity[,2], pch=21, bg=2)
