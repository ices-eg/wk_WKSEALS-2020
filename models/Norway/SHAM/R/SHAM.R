## Including capelin and cod biomass as covariates as covariates on fecundity 
## Martin Biuw and Tor Arne Øigård
## Modified from:
## Fitting state-space models to seal populations with scarce data
## Tor Arne Øigård and Hans J. Skaug
## ICES Journal of Marine Science, 2014
## Contact:  Martin Biuw (martin.biuw@hi.no) or Tor Arne Øigård (tor.arne.oigard@nr.no)
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
######################################

compile("R/SHAM.cpp")
dyn.load(dynlib("R/SHAM"))

###########
#Data part
###########
#Read in data
catch_data <- read.table("./Data/harpeast/catch_data.dat",header = FALSE)				#Catch data
pup_production <- read.table("./Data/harpeast/pup_production.dat",header = FALSE)			#Pup production estimates
fecundity <- read.table("./Data/harpeast/fecundity.dat",header = FALSE)				#Available fecundity data
#Pmatrix <- read.table("hessm.pma",header = FALSE)				#Birth ogives
#With collapsed maturity curves
Pmatrix <- read.table("./Data/harpeast/collapsed.ogi",header = TRUE)				#Birth ogives
Pmatrix <- c(Pmatrix[,2], rep(1, 4))
priors <- read.table("./Data/harpeast/priors.dat",header = FALSE)					#Priors used
years.of.prediction = 15											#Number of years to run projections

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
data$SSstart = 1985												#State-space process starts from this year
data$Cdata = as.matrix(catch_data)								#Add catch data
data$Nc = length(catch_data[,1])									#Number of years with catch data
data$pup_production = as.matrix(pup_production)					#Pup production estimates
data$Np = length(pup_production[,1])								#Number of years with pup production data
data$Fecundity = as.matrix(fecundity)							#Observed fecundity rates
data$Nf = length(fecundity$V1)									#Number of years with fecundity data
data$Pmat = as.matrix(Pmatrix)									#Birt ogives
data$Npred = years.of.prediction									#Number of years to run projections
data$priors = as.matrix(priors)									#Priors for estimated parameters
data$Npriors = length(priors$V1)									#Number of priors
data$xmax = data$Cdata[data$Nc,1]+data$Npred+1-data$SSstart+1	#Length of state-space process
data$CQuota = c(0,0)											#Catch level in future projections
data$cap = fish$CapTB/max(fish$CapTB)              #Capelin biomass
data$cod = fish$CodTB/max(fish$CodTB)              #Cod biomass

###############
#Parameter part
###############

#Initial values
Kinit = 2000000								#Initial population size
Minit = 0.09									#Natural adult mortality
M0init = 0.27								#Natural pup mortality
finit = 0.7									#Mean fecundity rate
ainit = 0.7									#AR(1) parameter
b1init = 0.4
b2init = 0.6
sigmainit = 0.85								#sigma

#Transform some parameters to ensure that estimates are > 0
parameters <- list()
parameters$logK= log(Kinit)  			
parameters$Mtilde= logit(Minit)				
parameters$M0tilde= logit(M0init)				
parameters$ftilde = logit(finit)
parameters$atilde = ainit
parameters$b1tilde = logit(b1init)
parameters$b2tilde = logit(b2init)
parameters$logSigma = log(sigmainit)
parameters$u=rep(0,length.out = (data$Nc+data$Npred-data$SSstart+data$Cdata[1,1]))


obj <- MakeADFun(data,parameters,random="u",DLL="SHAM")

obj$fn()
obj$gr()

system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))

rep<-sdreport(obj, getJointPrecision=TRUE)
rep.matrix <- summary(rep)
rep.rnames = rownames(rep.matrix)
indN0 = which(rep.rnames=="N0");indN0 <- indN0[-1]
indN1 = which(rep.rnames=="N1");indN1 <- indN1[-1]
indFt = which(rep.rnames=="Ft");indFt <- indFt[-1]
yrs = 1946:2035

N0 =rep.matrix[indN0,1]
N0sd = rep.matrix[indN0,2]
N1 = rep.matrix[indN1,1]
N1sd = rep.matrix[indN1,2]
Ft = rep.matrix[indFt,1]
Ftsd = rep.matrix[indFt,2]

#Extract parameters
Kest = exp(opt$par[1])          
Mest = ilogit(opt$par[2])       
M0est = ilogit(opt$par[3])      
fest = ilogit(opt$par[4])       
aest = opt$par[5]
b1est = ilogit(opt$par[6])
b2est = ilogit(opt$par[7])
sigmaest = exp(opt$par[8])    

save_results = FALSE
if(save_results){
#Save results
#Original model
#filename = "OriginalStateSpaceModel.Rdat"
#With prey data
  filename = "NewStateSpaceModel.Rdat"
  save(Kest = Kest,Mest = Mest,M0est=M0est,fest = fest,aest = aest,b1est = b1est,b2est=b2est,sigmaest = sigmaest,N1 = N1,N1sd = N1sd,N0 = N0,N0sd = N0sd,Ft = Ft,Ftsd = Ftsd,file = filename)
}
#allsd<-sqrt(diag(solve(rep$jointPrecision)))
#plsd <- obj$env$parList(par=allsd)
# X11("",9,7)
# plot(1946:2035,rep.matrix[indN1,1],type = "l",xlab = "Year",ylab = "Abundance",xlim = c(1946,2030),ylim = c(0,2000000),lwd = 3,col = "darkgreen")
# lines(1946:2035,rep.matrix[indN0,1],col = "blue",lwd = 3)
# lines(data$pup_production[,1],data$pup_production[,2],type = "p",col = "red",cex = 1.2,pch = 16)

#################
#Plot figures

#load("Data/Results/NewStateSpaceModel.Rdat")
dfnew=data.frame(Year=1946:2035,
              N1new=N1,
              N0new=N0,
              N1newl=N1-1.96*N1sd,
              N1newu=N1+1.96*N1sd,
              N0newl=N0-1.96*N0sd,
              N0newu=N0+1.96*N0sd,
              Fnew=Ft,
              Fnewl=Ft-1.96*Ftsd,
              Fnewu=F+1.96*Ftsd)

load("./Data/Results/OriginalStateSpaceModel.rdat")
dfold=data.frame(Year=1946:2035,
                 N1old=N1,
                 N0old=N0,
                 N1oldl=N1-1.96*N1sd,
                 N1oldu=N1+1.96*N1sd,
                 N0oldl=N0-1.96*N0sd,
                 N0oldu=N0+1.96*N0sd,
                 Fold=Ft,
                 Foldl=Ft-1.96*Ftsd,
                 Foldu=F+1.96*Ftsd)

df = data.frame(dfnew,dfold)

library(ggplot2)
theme_set(theme_minimal())

#Pup abundance
X11("",6,4)
ggplot(data=df,aes(x=Year))+
  geom_line(aes(y=N0new),colour="royalblue",size = 1.5)+
  geom_ribbon(aes(y = N0new,ymin = N0newl, ymax = N0newu),fill='lightblue',alpha=0.5)+
  #geom_line(aes(y=N0old),colour='red',size=1.5)+
  #geom_ribbon(aes(y = N0old,ymin = N0oldl, ymax = N0oldu),fill='lightpink',alpha=0.5)+
  labs(x='Year',y='Pup abundance',title='')

#Fecundity
X11("",9,7)
ggplot(data=df,aes(x=Year))+
  geom_line(aes(y=Fnew),colour="royalblue",size = 1.5)+
  #geom_ribbon(aes(y = Fnew,ymin = Fnewl, ymax = Fnewu),fill='lightblue',alpha=0.5)+
  #geom_line(aes(y=Fold),colour='red',size = 1.5)+
  #geom_ribbon(aes(y = Fold,ymin = Foldl, ymax = Foldu),fill='lightpink',alpha=0.5)+
  labs(x='Year',y='Modelled fecundity',title='')