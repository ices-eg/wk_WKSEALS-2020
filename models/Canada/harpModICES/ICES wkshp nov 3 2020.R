rm(list = ls())         #Delete all the current files in the workspace

##################################################################
#### Define the directory where the files are read and written ###
##################################################################

root <- "../wk_WKSEALS-2020/" # define the location of the folder with the wk-WKSEALS-2020 GitHub files

rep <- paste0(root, "\\models\\Canada\\harpModICES")
datDir <- paste0(root, "\\data\\Canada\\harpModICES\\")        
outdir<-paste(rep, "\\modelAMKcor_output\\", sep="") # set name of the output folder
if (!file.exists(substring(outdir,1, nchar(outdir)-1))) {dir.create(outdir)}  #create directory if it does not exist


##############################################################################################
###---------------------------------------------------------------------------------------####
###                      LOAD AND ORGANIZE THE MODEL PARAMETERS                           ####
###---------------------------------------------------------------------------------------####
##############################################################################################

start.time1 <- Sys.time()
start.time1

setwd(datDir)

### 1. SMOOTHED REPRODUCTIVE RATES
  # Load data
  data<-read.csv("temppr_2019.csv", sep=";", header=F)
  colnames(data) <- c("years", "age", "total", "npreg")
  data$total[data$total == 0] <- NA
  nyears<-length(unique(data$years))

  # create 2 empty matrices
  predictyears <- seq(1952, 2020, 1) # set years where smoothed data are needed
  span<-range(predictyears)[2]-range(predictyears)[1]+1
  preg.allyears.log <- matrix(NA, length(unique(data$age)), length(predictyears))
  colnames(preg.allyears.log) <- as.character(predictyears)
  se.allyears.log <- matrix(NA, length(unique(data$age)), length(predictyears))
  colnames(se.allyears.log) <- as.character(predictyears)

  # calculate smoothed data
  library(locfit)
  degoffree <- 2 # set degree of freedom
#  nearestn <- c(0.75,0.8,0.95,0.95,0.95)   # set nearest neighbor values chosen with LCV and AIC diagnostics measures (see Diagnostics part)
                                          # for age 4,5,6,7, and 8 respectively
  nearestn <- c(0.75,0.8,0.95,0.95,0.95)   # set nearest neighbor values chosen with LCV and AIC diagnostics measures (see Diagnostics part)
                                          # for age 4,5,6,7, and 8 respectively

  for (i in 1:length(unique(data$age))){
    subbyage <- subset (data, data$age == unique(data$age)[i])
    fit <- locfit(npreg~lp(years, deg = degoffree, nn = nearestn[i]), weights = total, data = subbyage, family = "binomial")
    preg.allyears.log[i,] <- preplot(fit, data.frame(years = predictyears))$fit
#    se.allyears.log[i,] <- preplot(fit, data.frame(years = predictyears), band="local")$se
    
    se.allyears.log[i,] <- preplot(fit, data.frame(years = predictyears), band="local")$se * 1.0
    
  }

  # for information give matrix of smoothed value from preg rates for each predicted year
  nb_ageplus <- 18  # set number of classes age+

  preg.young <- matrix(0,3,length(predictyears))
  preg.all <- expit(preg.allyears.log)
  preg.plus <- matrix(preg.all[length(unique(data$age)),],nb_ageplus,length(predictyears), byrow=T)

  preg.smooth <- rbind(preg.young, preg.all, preg.plus)

  preg.young.CIs <- matrix(0,3, length(predictyears))
  preg.all.CIupper <- expit(preg.allyears.log + se.allyears.log*1.96)
  preg.all.CIlower <- expit(preg.allyears.log - se.allyears.log*1.96)
  preg.plus.CIupper <- matrix(preg.all.CIupper[length(unique(data$age)),],nb_ageplus, length(predictyears), byrow=T)
  preg.plus.CIlower <- matrix(preg.all.CIlower[length(unique(data$age)),],nb_ageplus, length(predictyears), byrow=T)
  preg.smooth.CIupper <- rbind(preg.young.CIs, preg.all.CIupper, preg.plus.CIupper)
  preg.smooth.CIlower <- rbind(preg.young.CIs, preg.all.CIlower, preg.plus.CIlower)

  # Plots smoothed preg.rates and 95CIs for age classes 4 to 8+
  pdf(file=paste(outdir,"Smoothed preg rates and 95CIs.pdf",sep=""))
  par(mfrow=c(2,3))
  for (i in 4:8){
    plot(preg.smooth[i,]~colnames(preg.smooth), ylim = c(0,1), xlab="years", ylab="preg.rate", main=paste("age", i), type="l", col=1)
    lines(preg.smooth.CIupper[i,]~colnames(preg.smooth),col="blue", lty=2)
    lines(preg.smooth.CIlower[i,]~colnames(preg.smooth), col="blue", lty=2)
    points((data$npreg/data$total)[data$age == i]~data$years[data$age == i], pch=16, col="red", cex=0.7)
  }
  dev.off()


  # Create an output file including original data (n, n_preg) as well as mean preg rate and CIs evaluated by the smoother
  # for each predictyears and for each age class (note: values for "young" and "age plus" classes are not presented as they are equal
  # to 0 for "young" classes and to the values of the last age class with data (i.e. class "age 8") for the "age plus" classes.
  output <- data.frame()
  for (i in 1:nrow(preg.all)){
    output <- rbind (output, data.frame(Year = predictyears, Age = unique(data$age)[i], n = NA, n_preg = NA, Smooth_preg_rate = unname(preg.all[i,]), CI95lower = unname(preg.all.CIlower[i,]), CI95upper = unname(preg.all.CIupper[i,])))
    output$n[which(output$Year %in% data$years[data$age == unique(data$age)[i]])] <- data$total[data$age == unique(data$age)[i]]
    output$n_preg[which(output$Year %in% data$years[data$age == unique(data$age)[i]])] <- data$npreg[data$age == unique(data$age)[i]]
  }


### 2. ESTIMATE N FOR MISSING VALUES IN BINOMIAL DISTRIBUTION BASED ON SD FROM SMOOTHERS

    nresample <- 10000  # number of resampling (see below)

    # Progress bar
    library(tcltk)
    maxbar1 <- nresample
    pb1 <- tkProgressBar(title = "Resampling smoothed pregnancy rate distributions", min = 0, max = maxbar1, width = 400)


  rubix <- array(dim=c(26,span,nresample))        #  create an empty array to store the resampling results

  for(r in 1:nresample) {                       # this loops resamples 10,000 times within the smoothers' distribution to generate an SD

    setTkProgressBar(pb1, r, label=paste("Percent done",round(r/maxbar1*100), "%")) # update progress bar

    preg.allyears.loop <- matrix(NA, length(unique(data$age)), length(predictyears))
    colnames(preg.allyears.loop) <- as.character(predictyears)
    for (m in 1: length(unique(data$age))){
      preg.allyears.loop[m,] <- expit(rnorm(length(predictyears), mean = preg.allyears.log[m,], sd = se.allyears.log[m,]))
    }

    preg.plus.loop <- matrix(NA, nb_ageplus, length(predictyears))
    colnames(preg.plus.loop) <- as.character(predictyears)
    for (m in 1:18) {
      preg.plus.loop[m,] <- expit(rnorm (length(predictyears), mean = preg.allyears.log[length(unique(data$age)),], sd = se.allyears.log[length(unique(data$age)),]))
    }

     rubix[,,r] <- rbind(preg.young, preg.allyears.loop, preg.plus.loop)
  }
  close(pb1)

  # from this SD, we can estimate the size n of a binomial distribution binom(n,p) that would have the same SD
  sd.mat <- apply(rubix,c(1,2),sd)
  p.mat <- preg.smooth
  ncalc.mat <- round(p.mat*(1-p.mat)/(sd.mat^2))
  ncalc.mat[is.nan(ncalc.mat)] <- 1
  ncalc.mat[ncalc.mat<1] <- 1

  # Include the estimated size n in the output describing smoother predictions
  n_estim<-vector()
  for (i in (4:8)){
    n_estim <- c(n_estim, ncalc.mat[i,])
  }
  output <- data.frame(output, n_estim = as.numeric(t(n_estim)))

  write.csv(output, paste(outdir,"output_smoother_model.csv", sep=""), row.names = FALSE)

### 3.  SMOOTHED REPRODUCTIVE RATES (to use in the binomial distribution)
  npreg <- data[,4]
  npreg.mat <- matrix(npreg,5,nyears, byrow=T)
  colnames(npreg.mat) <- data[1:nyears,1]

  nsize <- data[,3]
  nsize.mat <- matrix(nsize,5,nyears, byrow=T)
  colnames(nsize.mat) <- data[1:nyears,1]

  preg.young <- matrix(0,3,nyears)
  preg.plus <- matrix(npreg.mat[5,],18,nyears, byrow=T)
  preg.all <- rbind(preg.young, npreg.mat, preg.plus)

  size.young <- matrix(0,3,nyears)
  size.plus <- matrix(nsize.mat[5,],18,nyears, byrow=T)
  size.all <- rbind(size.young, nsize.mat, size.plus)

  # create 2 empty matrices with all years and age classes
  preg.allyears <- matrix(NA, 26, span)
  colnames(preg.allyears) <- as.character(predictyears)
  size.allyears <- matrix(NA, 26, span)
  colnames(size.allyears) <- as.character(predictyears)

  # find which years have data
  match.years <- match(colnames(preg.all),colnames(preg.allyears))

  # replace values in allyears matrix by preg counts for years with data
  preg.allyears[,match.years] <- preg.all
  size.allyears[,match.years] <- size.all

  # transform pregnant counts into proportions
  preg.allyears <- preg.allyears/size.allyears

  # turn proportions of 100% into 95% (otherwise, there is no variability) # note that Binom(22,0.95): mean = 0.95, 2.5%-97.5% quantiles = 0.82-1.00

  preg.allyears[preg.allyears==1] <- 0.95

  # remove NaN's
  preg.allyears[is.nan(preg.allyears)]<-NA

  # find cells with NA or sample size < 40
  toreplace <- which(is.na(size.allyears) | size.allyears < 49)

  # find cells with NA or sample size < 40 for 8+ preg rate       
  smoothedvalues <- which(is.na(size.allyears[8,]) | size.allyears[8,] < 49)

  # replace these cells with smoothed preg rates and sample size = n as calculated above from SD of smoothed data
  preg.allyears[toreplace] <- preg.smooth[toreplace]
  size.allyears[toreplace] <- ncalc.mat[toreplace]

### Include all parameters in one object

  sealmod.params<-list()     #The parameters of the model are defined as a list (suite of objects)

  # Final age specific reproductive rates
  sealmod.params$preg.rate <- preg.allyears

  # Final age specific sample sizes for reproductive rates
  sealmod.params$preg.size <- (size.allyears )

  # Proportion of pups surviving an unusual mortality event (1952-2005)

  
  sealmod.params$prop.surv <- as.matrix(read.table("prop-surv-1952_updatesept2019f.csv", sep = ";", header = TRUE))  
  
   # food factor (1952-2005)
#  sealmod.params$FF <- as.matrix(read.table("FallCapelin.csv", sep = ";", header = TRUE))

  sealmod.params$FF <- as.matrix(read.table("CEIndex.csv", sep = ";", header = TRUE))
  
#  sealmod.params$FF <- as.matrix(read.table("CElindex.csv", sep = ";", header = TRUE))

  # Pup production mean estimates based on mark-recapture experiments and aerial surveys
  sealmod.params$pup.prod <- as.matrix(read.table("pup-prod-est-1952.csv", sep = ";", header = TRUE, na.strings = "NA"))
    #sealmod.params$pup.prod[57] <-1630000

  # Pup production SE estimates based on mark-recapture experiments and aerial surveys
  sealmod.params$pup.prod.se <- as.matrix(read.table("pup-prod-est-se-1952.csv", sep = ";", header = TRUE, na.strings = "NA"))
   # sealmod.params$pup.prod.se[57] <-110381

  # Initial vector of abundance (year 1952)
  sealmod.params$init.pop <- as.matrix(read.table("initial-pop-1952.csv", sep = ";", header = FALSE))

  # Observed reproductive rate for class 8+
  sealmod.params$obs.preg.8plus <- rep(NA, length(1952:2020))
  sealmod.params$obs.preg.8plus[match.years] <- npreg.mat[5,] / nsize.mat[5,]

  # Observed sample size for class 8+
  sealmod.params$obs.pregsize.8plus <- rep(NA, length(1952:2020))
  sealmod.params$obs.pregsize.8plus[match.years] <- nsize.mat[5,]

  # Natural mortality multiplier for first year seals
 # gamma.pup <- M
  M1<- 0.032
  raw.remov <- as.matrix(read.table("raw-removal-1952.csv",sep = ";", header = TRUE))
  prop.arctic.pup <- 0.034
  prop.green.pup <- scan("prop-green-pup-1952.csv",sep = ";",quiet=T)




##############################################################################################
###---------------------------------------------------------------------------------------####
###                           FIT TO PAST AND PRESENT DATA                                ####
###---------------------------------------------------------------------------------------####
##############################################################################################

start.time2 <- Sys.time()

sealmod.sim.f <- function(params, last.data.year, preg.rate, pup.prod){  # The result of the function is the pup production, the total population and the last year population vector

  years <- length(1952:last.data.year)

  ## Params(alpha, M, K)
  alpha <- params[1]
  M <- params[2]
  K <- params[3]


  ## Initial vector of abundance
  seal <- numeric(0)
  seal[2:26] <-  sealmod.params$init.pop[2:26] * alpha                 #seal will change at each loop
  seal[1] <- sum(preg.rate[2:26,1] * seal[2:26]/2)      #Pups in first year (1952)

  ## Setup the output simulation matrix
  seal.sim <- matrix(ncol= years+1, nrow= 26)
  seal.sim[,1] <- seal    #year 1952

  ## Initialize a vector for preg.rates estimated for classe 8+
  preg.8plus.sim <- rep(NA, years+1)

  ## Initialize a vector for pup.prod estimated for each age classe
  pup.prod.byage <- matrix(ncol= years+1, nrow= 6)

#  SandLadult<-c(2,2,2,1)

  for (year in 1:years) {

    SandLpup<-c(ifelse(year<=31,1/0.99,1/0.95),2,2,1)          #for 1960 model this value is set at 23 for 1983
    SandLadult<-c(ifelse(year<=31,1/0.6,2),2,2,1)
    prop.age.class<-prop.table(seal.sim[2:26,year])

    adult.kill <- sum(c(raw.remov[year,5],raw.remov[year,3]*(1-prop.green.pup[year]),raw.remov[year,2]*(1-prop.arctic.pup),raw.remov[year,6])*SandLadult)*prop.age.class
    pup.kill <- sum(c(raw.remov[year,4],raw.remov[year,3]*prop.green.pup[year],raw.remov[year,2]*prop.arctic.pup,raw.remov[year,7])*SandLpup)
    remov <- append(pup.kill,adult.kill) 

    ## Add density dependence effect on preg.rates of 8+
#    preg.8plus.sim[year+1] <- ifelse(0.88 * (1 - (sum(seal.sim[,year])/K^2.4))< 0, 0, (0.88 * (1 - (sum(seal.sim[,year])/K)^2.4))) # 0.88 is considered as max preg rate (can be defined as max preg rate when pop is low)
    
   preg.8plus.sim[year+1] <- ifelse(0.88 *(1 - (sum(seal.sim[,year])/(K*sealmod.params$FF[,year+1]))^2.4) < 0, 0,(0.88 * (1 - (sum(seal.sim[,year])/(K*sealmod.params$FF[,year+1]))^2.4)))      #  0.88 is considered as max preg rate (can be defined as max preg rate when pop is low)
      
 # preg.8plus.sim[year+1] <- ifelse(sealmod.params$FF[,year]*0.88 * (1 - (sum(seal.sim[,year])/K)^2.4) < 0,0,sealmod.params$FF[,year]*0.88 * (1 - (sum(seal.sim[,year])/K)^2.4))
 # 0.88 is considered as max preg rate (can be defined as max preg rate when pop is low)
    
#   preg.rate[8:26, year+1] <- rep (preg.8plus.sim[year+1], 19)   # Uncomment this line to allow the model to use preg.rate for 8+ modified by the density dependance formula in the pup prod fitting
#   preg.8plus.sim[year+1] <- preg.rate[8,year+1]  # When not using the two previous lines, uncomment this one


    # Seal numbers at age 1
    seal.age1 <- (seal[1] * (sealmod.params$prop.surv[,year]) - remov[1]) *( exp(- M) *(1-((sum(seal.sim[,year])/K)^2.4)) )     #seal[1,] = pups
#   seal.age1 <- (seal[1] * (sealmod.params$prop.surv[,year]^2) - remov[1]) *( exp(- M) *(1-((sum(seal.sim[,year])/(K*sealmod.params$FF[,year]))^2.4)) )
#   seal.age1 <- ((seal[1] * (sealmod.params$prop.surv[,year]^2) - remov[1]) *( exp(- M) *(1-((sum(seal.sim[,year])/(K*sealmod.params$FF[,year]))^2.4)) ) ) *  sample(rnorm(1,1,0.001),1,replace=TRUE)
    
#     seal.age1 <- (seal[1]  - remov[1]) *( exp(- M) *(1-((sum(seal.sim[,year])/K)^2.4)) ) 
     
#    seal.age1 <- (seal[1] * sealmod.params$prop.surv[year] - remov[1]) * exp(-M ) *(1-((sum(seal.sim[,(year)])/(K*sealmod.params$FF[,year]))^2.4))

    # Seal numbers for age > 1 and < 25
    seal.num <- (seal[2:24]*exp(-.03/2) - remov[2:24]) * exp(-.03/2)
    
#     seal.num <- (seal[2:24]*exp(-M/2) - remov[2:24]) * exp(-M/2)

    # Seal numbers for age 25+
    seal.old <- ((seal[26]+seal[25]) * exp(-.03/2) - remov[25] - remov[26]) * exp(-.03/2)


    # Pups produced
    pup.prod4 <-  preg.rate[4, year+1] * (c(seal.age1, seal.num, seal.old)/2)[3]   # pup production for 3 years old individuals
    pup.prod5 <-  preg.rate[5, year+1] * (c(seal.age1, seal.num, seal.old)/2)[4]
    pup.prod6 <-  preg.rate[6, year+1] * (c(seal.age1, seal.num, seal.old)/2)[5]
    pup.prod7 <-  preg.rate[7, year+1] * (c(seal.age1, seal.num, seal.old)/2)[6]
    pup.prod8plus <-  sum(preg.rate[8:26, year+1] * (c(seal.age1, seal.num, seal.old)/2)[7:25])
    pup.not8plus <- sum(pup.prod4, pup.prod5, pup.prod6, pup.prod7)
    pup.prod.byage[,(year+1)] <- c(pup.prod4, pup.prod5, pup.prod6, pup.prod7, pup.prod8plus, pup.not8plus)

#    seal.pups <- sum(preg.rate[2:26,year+1] * c(seal.age1, seal.num, seal.old)/2) *sqrt(sealmod.params$prop.surv[year+1])      # mortality at birth from poor ice
    
    seal.pups <- sum(preg.rate[2:26,year+1] * c(seal.age1, seal.num, seal.old)/2)
    
    seal.pups[seal.pups<0]=0.001

    # New vector of seals
    seal[1] <- seal.pups
    seal[2] <- seal.age1
    seal[3:25] <- seal.num
    seal[26] <- seal.old

    seal[!is.finite(seal)] <- 0.001
    seal[seal < 0.001]  <- 0.001

    seal.sim[,(year+1)] <- seal
  }
  last.pop.vec <- seal.sim[,years+1]
  sealmod.sim.f.res <- matrix(c(seal.sim[1,],apply(seal.sim,2,sum)),years+1,2)         #to get the total number of seals as well
  age.struct <- rbind(seal.sim[1:15,], apply(seal.sim[16:26,], 2, sum))
  list(sealmod.sim.f.res,last.pop.vec, preg.8plus.sim, pup.prod.byage, age.struct)
}


################################## OPTIMIZATION FUNCTION ##################################################################################

sealmod.fit.f <- function(params, last.data.year, preg.rate, pup.prod){           #You give initial parameters and the last data year and you get the value of the fonction which is the value of the SS for those 2 paramters

  if (all(params >= c(0.1,0.05,9e+6) & params <= c(0.4,0.7,1.4e+7))){  # check if parameters fall in a correct range

    sim.data <- sealmod.sim.f(params, last.data.year, preg.rate, pup.prod) ## call simulation part

    # Objective function 1  (fit pup population size)
    sim.pup.data <- sim.data [[1]][,1]   #Get the pup production from the sealmod.sim.f function
    obs.pup.data <- as.vector(pup.prod)[1:(length(1952:last.data.year)+1)]
    var.obs.pup.data <- as.vector((sealmod.params$pup.prod.se)^2)[1:(length(1952:last.data.year)+1)]

    objective1 <- sum((sim.pup.data - obs.pup.data)^2, na.rm=TRUE)

    # Objective function 2 (fit 8+ pregnancy rate)
    sim.preg.8plus <- sim.data [[3]]
    obs.preg.8plus <- preg.rate[8,]

    obs.preg.8plus[smoothedvalues] <- NA # remove smoothed values for the fitting part.

    obs.pregsize.8plus <- sealmod.params$obs.pregsize.8plus[1:(length(1952:last.data.year)+1)]
#   obs.preg.8plus[obs.preg.8plus == 1 & !is.na(obs.preg.8plus)] <- 0.999  #avoid problem on the logit scale (logit(1) = Inf)
#   sim.preg.8plus[sim.preg.8plus == 1 & !is.na(sim.preg.8plus)] <- 0.999
#   obs.preg.8plus[obs.preg.8plus == 0 & !is.na(obs.preg.8plus)] <- 0.001  #avoid problem on the logit scale (logit(0) = -Inf)
#   sim.preg.8plus[sim.preg.8plus == 0 & !is.na(sim.preg.8plus)] <- 0.001
#   library(boot)
#   objective2 <- sum((logit(sim.preg.8plus) - logit(obs.preg.8plus))^2, na.rm=TRUE)  # using logit because of the binomial characteristic of the pregnancy rate
    objective2 <- sum((sim.preg.8plus - obs.preg.8plus)^2, na.rm=TRUE)

    # Final objective value
#   objective3 <- objective1/var(obs.pup.data, na.rm = TRUE) + objective2/var(logit(obs.preg.8plus), na.rm = TRUE)
#     objective3 <- (objective1/var(obs.pup.data, na.rm = TRUE)) + objective2/var(obs.preg.8plus, na.rm = TRUE)
    objective3 <- (1*objective1/var(obs.pup.data, na.rm = TRUE)) + objective2/var(obs.preg.8plus, na.rm = TRUE)       # multiply by 2 on survey to weight it for heavy ice scenario 
#    objective3 <- objective1/var(obs.pup.data, na.rm = TRUE)
  }else{
    objective3 <- 1e+10
  }
  return(objective3)
}



#################### SAMPLING AND PARALLEL RUN ############################################################################

SealModel <- function(nbrep, params, last.data.year, preg.cor=0.90, NoConvLimit = 2, clust= min(c(nbrep, detectCores(logical=TRUE)-1)), ShowAgeStruct = c(TRUE, TRUE)){


  ## RESAMPLING AND OPTIMIZATION PROCESS ######################################################################

  optimization <- function(params, last.data.year, sealmod.params, preg.cor){

    ## Resampling pregnancy rates into binomials with defined correlation
    require(copula)
    norm.cop <- normalCopula(preg.cor, dim=8)
    rcop <- rCopula(1,norm.cop)
    preg.sample <- matrix(nrow=26,ncol=span)
    for (m in 1:span){
      preg.sample[1:8,m] <- t(sapply(seq(1,8), FUN = function(x){qbinom(rcop[,x],size=sealmod.params$preg.size[x,m], prob = sealmod.params$preg.rate[x,m])/sealmod.params$preg.size[x,m]}))
      preg.sample[9:26,m] <- preg.sample[rep(8,18),m]    }

    ## Resample the pup estimates
    pup.sample <- rnorm(ncol(sealmod.params$pup.prod),sealmod.params$pup.prod,sealmod.params$pup.prod.se)


    ## Optimization
    out<- optim(par=params, fn=function(x) sealmod.fit.f(x, last.data.year, preg.sample, pup.sample), method = "Nelder-Mead", hessian = F, , control=list(parscale=c(0.01,0.001,50000), maxit = 2000))   #  The optimisation, Nelder-Mead is the best and fastest optimisation method available in R but is a little different than Excel

    # Repeat optimisation if necessary to increase the fit (also avoid convergence on aberrant value)
    o<-1; go=-1
    while (go<0){
      prevout <- out$value
      out <- optim(par=out$par, fn=function(x) sealmod.fit.f(x, last.data.year, preg.sample, pup.sample), method = "Nelder-Mead", hessian = F, control=list(parscale=c(0.01,0.001,50000), maxit = 2000))
      difval <- prevout - out$value
      o <- o + 1
      go<-ifelse((out$value < 70 & difval < 0.01) | o > 50, 1,-1)
      if (out$value > 70) out$par <- c(rnorm(1,params[1], 0.01), rnorm(1, params[2], 0.01), rnorm(1, params[3], 1e+6))
    }

    if (out$convergence == 0){
      numbers<-sealmod.sim.f(out$par,last.data.year, preg.sample, pup.sample)   # create a population using the optimised parameters
      return(list(convergence = 0,
      pup = numbers[[1]][,1],       #  this line and below: save results
      tot = numbers[[1]][,2],       #
      last.pop.vec = numbers[[2]],	#
      pregsampl = preg.sample,
      preg8plus = numbers[[3]],
      pupsample = pup.sample,
      PupProdByAge = numbers[[4]],
      AgeStruct = numbers[[5]],
      AgeStructProp = t(apply(numbers[[5]], 1, function(x) x / apply(numbers[[5]],2,sum))),

      result.alpha = out$par[1],    #
      result.M = out$par[2],			  #
      result.K = out$par[3],			  #

      result.value = out$value))  			#  end saving results
    }else{
      return(list(preg.rate, pup.prod))
    }
  }


  ## PARALLEL RUN AND SAVE RESULTS ######################################################################

  ## Create object to save results
  result.alpha<-numeric(0)  #initialised the result vectors
  result.M<-numeric(0)      #idem
  result.K<-numeric(0)      #idem

  result.value<-numeric(0)  #idem
  
  pup<-matrix(NA,nbrep,length(1952:last.data.year)+1)  #idem
  tot<-matrix(NA,nbrep,length(1952:last.data.year)+1)  #idem
  last.pop.vec<-matrix(NA,26,nbrep)    #idem
  pregsampl <- list()
  preg8plus <- matrix(NA,nbrep,length(1952:last.data.year)+1)
  PupProdByAge <- array(NA, dim=c(6, length(1952:last.data.year)+1, nbrep))
  pupSampl<-matrix(NA,nbrep,length(1952:last.data.year)+1)  #  will be used to test the distribution of the residuals of the fitting
  preg8Sampl<-matrix(NA,nbrep,length(1952:last.data.year)+1)  #  will be used to test the distribution of the residuals of the fitting
  AgeStruct <- array(NA, dim=c(16, length(1952:last.data.year)+1, nbrep))
  AgeStructProp <- array(NA, dim=c(16, length(1952:last.data.year)+1, nbrep))

  ## Create a list that will contain parameters inducing non-convergence of the model
  sealmod.crash <<- list(preg.sample = list(), pup.sample = list(), crash.nb = 0, nbrep = nbrep, param = params, preg.cor = preg.cor, last.data.year = last.data.year)  # Note, this list is created in the base environment (i.e not in the function environment)
  
  
  ## Run parallel
  library(parallel)
  
  if (any(class(clust) == "numeric")){ t1 <- system.time(cl <- makeCluster(clust))
  }else if(
  any(class(clust) == "cluster")){ cl <- clust; t1 <- system.time(0)
  }else{
  stop("parameter clust should be numeric or an object of class cluster created by the function makeCluster (package parallel)")
  }

  trash <- clusterEvalQ(cl, rm(list=ls())) # remove all object in the cluster if any
  clusterExport(cl, c("params",            # send all the information (functions and parameters) to the cluster
                      "sealmod.fit.f",
                      "sealmod.sim.f",
                      "optimization",
                      "last.data.year",
                      "raw.remov",
                      "prop.green.pup",
                      "prop.arctic.pup",
                      "M1",
                      "smoothedvalues",
                      "sealmod.params",
                      "preg.cor",
                      "span"), envir=environment())


    run <- function(x){
      res<-clusterApplyLB(cl ,1:nbrep ,fun=function(i) optimization(params, last.data.year, sealmod.params, preg.cor))
      trial <- 0
      while (sum(unlist(lapply(res, function(x) x$convergence==0))) < nbrep & trial < NoConvLimit){  # repeat until the number of runs with convergence is obtained or until it attains the maximal number of trial allowed
        res <- c(res, clusterApplyLB(cl ,1:(nbrep-sum(unlist(lapply(res, function(x) x$convergence==0)))) ,fun=function(i) optimization(params, last.data.year, sealmod.params)))
        trial <- trial + 1
      }
      return(res)
    }
  
    t2<-system.time(tryCatch(resfinal <- run(), finally = if(any(class(clust) == "numeric")) stopCluster(cl)))  # All the parallel calculations are run here and results stored in the resfinal object
    
    
  ## Organize the results - distribute the values in the corresponding objects
     
    for (i in 1:length(resfinal)){
      out <- resfinal[[i]]
      if (out$convergence == 0){
        pup[i,]<-out$pup
        tot[i,]<-out$tot
        last.pop.vec[,i]<-out$last.pop.vec
        pregsampl[[i]] <- out$pregsampl
        preg8plus[i,] <- out$preg8plus
        PupProdByAge[,,i] <- out$PupProdByAge
        pupSampl[i,] <- c(out$pupsample, NA)
        preg8Sampl[i,] <- out$pregsampl[8,]
        preg8Sampl[i,smoothedvalues] <- NA    # save preg sample fitted (remove smoothed values) ... to be used for testing the distribution of residuals
        AgeStruct[,,i] <- out$AgeStruct
        AgeStructProp[,,i] <- out$AgeStructProp

        result.alpha[i]<-out$result.alpha
        result.M[i]<-out$result.M
        result.K[i]<-out$result.K

        result.value[i]<-out$result.value  			
        i<-i+1
      }else{
        sealmod.crash$crash.nb <<- sealmod.crash$crash.nb + 1    # save parameters that induced the non-convergence of the model
        sealmod.crash$pup.sample[[sealmod.crash$crash.nb]] <<- pup.sample
        sealmod.crash$preg.sample[[sealmod.crash$crash.nb]] <<- preg.sample
        if (sealmod.crash$crash.nb > NoConvLimit) { stop(paste("Optimization did not succeed to converge after", NoConvLimit, "trials", "\n\nNote: By default the number of trials is equal to half the number of model repetition", "\nIt can be set directly using the NoConvLimit argument of the function")) }
      }
    }
  result<-matrix(c(result.alpha,result.M,result.K, result.value),nbrep,4,dimnames=list(c(1:nbrep),c("Alpha","M","K","Value")))   # create the optimisation parameters result matrix
  
   
  ## Output parameters stats 
  
  ParamStats <- sapply(1:(ncol(result)-1), FUN = function(x){data.frame(Mean=mean(result[,x]), Median=median(result[,x]), Sd = sd(result[,x]), Q0.025 = quantile(result[,x], probs=0.025), Q0.975 = quantile(result[,x], probs=0.975))})
  colnames(ParamStats) <- c("Alpha", "M", "K")
  write.csv(ParamStats, paste(outdir,"ParamStats_for_", nbrep, "_iterations_InitParams_", paste(params, collapse="_"), "_PregCor_", preg.cor, ".csv", sep=""))


  ## Output Age structure

  if (ShowAgeStruct[1]) {
    # Age structure summary (seals numbers by age class)
    AgeStructSummary <- list(prob.025 = apply(AgeStruct, c(1,2), function(x) quantile(x, prob = 0.025)), prob.25 = apply(AgeStruct, c(1,2), function(x) quantile(x, prob = 0.25)), med = apply(AgeStruct, c(1,2), median),
                   prob.50 = apply(AgeStruct, c(1,2), function(x) quantile(x, prob = 0.50)), prob.75 = apply(AgeStruct, c(1,2), function(x) quantile(x, prob = 0.75)), prob.975 = apply(AgeStruct, c(1,2), function(x) quantile(x, prob = 0.975)))
    for(i in 1:length(AgeStructSummary)){ AgeStructSummary[[i]] <- data.frame(c(seq(0,14), "15+"), AgeStructSummary[[i]]) ; colnames(AgeStructSummary[[i]]) <- c("", seq(1952, last.data.year+1))}


    # Save raw numbers
  
    write.table(rbind("### Median"), paste(outdir,"AgeStruct_for_", nbrep, "_iterations.csv"), row.names=F, col.names=F)
    suppressWarnings(write.table(AgeStructSummary$med, paste(outdir,"AgeStruct_for_", nbrep, "_iterations.csv"), sep=",", append=TRUE, row.names=F))

    write.table(rbind("","### Prob .025"), paste(outdir,"AgeStruct_for_", nbrep, "_iterations.csv"), append=T, row.names=F, col.names=F)
    suppressWarnings(write.table(AgeStructSummary$prob.025, file=paste(outdir,"AgeStruct_for_", nbrep, "_iterations.csv"), sep=",", row.names=F, append = TRUE))

    write.table(rbind("", "### Prob .25"), paste(outdir,"AgeStruct_for_", nbrep, "_iterations.csv"), append=T, row.names=F, col.names=F)
    suppressWarnings(write.table(AgeStructSummary$prob.25, file=paste(outdir,"AgeStruct_for_", nbrep, "_iterations.csv"), sep=",", row.names=F, append = TRUE))

    write.table(rbind("", "### Prob .75"), paste(outdir,"AgeStruct_for_", nbrep, "_iterations.csv"), append=T, row.names=F, col.names=F)
    suppressWarnings(write.table(AgeStructSummary$prob.75, file=paste(outdir,"AgeStruct_for_", nbrep, "_iterations.csv"), sep=",", row.names=F, append = TRUE))

    write.table(rbind("", "### Prob .975"), paste(outdir,"AgeStruct_for_", nbrep, "_iterations.csv"), append=T, row.names=F, col.names=F)
    suppressWarnings(write.table(AgeStructSummary$prob.975, file=paste(outdir,"AgeStruct_for_", nbrep, "_iterations.csv"), sep=",", row.names=F, append = TRUE))
  }


  #########################################################################################################
  ################################## Plots ################################################################
  #########################################################################################################
  
  windows(18,12)
  par(mfrow = c(2,4), mar = c(4.5,5,3,1))
  paramName <- c("Alpha", "Adult mortality (M)", "Carrying capacity (K)")
  
  
  ## Plot the distribution of each parameter value along with basic statistics (mean, median, CI)
  for (i in 1:4){
    if (!is.na(result[1,i])& i<(ncol(result))){# & i<(ncol(result))){  # remove & i<(ncol(result)), in this line if you want to plot the graph for the objective function
      param.res = density(result[,i])
      hist(result[,i],  breaks = 20,
       main = ifelse(i==ncol(result), "Objective function", paste(paramName[i], "\nmean =", ifelse(mean(result[,i])<10, round(mean(result[,i]),4), round(mean(result[,i]),0)),
       "; median =", ifelse(median(result[,i])<10, round(median(result[,i]),4), round(median(result[,i]),0)),
       "\n95%CI", ifelse(quantile(result[,i], probs=0.025)<10, round(quantile(result[,i], probs=0.025),4), round(quantile(result[,i], probs=0.025),0)),
       "-", ifelse(quantile(result[,i], probs=0.975)<10, round(quantile(result[,i], probs=0.975),4), round(quantile(result[,i], probs=0.975),0)))),
       ylab = 'Density', cex.axis = 1.3, cex.lab = 1.3, probability = T)
      lines(param.res, col = "black", lwd = 3)
    }else{
      plot(1, type="n", xlab="", ylab="", bty="n", axes = FALSE)
    }
  }
  
  ## Plot the simulated 8+ preg rate vs Years
  plot(apply(preg8plus,2,median)~seq(1952,(last.data.year+1)), main="Preg 8+", xlab="years", ylab="rate", ylim=c(0,1), type="b")
  points(as.vector(sealmod.params$obs.preg.8plus)[1:(length(1952:last.data.year)+1)]~seq(1952,(last.data.year+1)), col="red", pch=16)
  points(apply(preg8plus,2,FUN=function(x){quantile(x, probs=0.025, na.rm=TRUE)})~seq(1952,(last.data.year+1)), type="l", col="blue")
  points(apply(preg8plus,2,FUN=function(x){quantile(x, probs=0.975, na.rm=TRUE)})~seq(1952,(last.data.year+1)), type="l", col="blue")
  legend("bottomleft", c("Simulated 8+ preg rate (median)", "0.025 and 0.975 quantiles", "Observations"), lty = c(1,1,-1), pch = c(1,-1,16), col=c("black","blue","red"))
  
  ## Plot the simulated pup production vs Years
  plot(apply(pup,2,median)~seq(1952,(last.data.year+1)), main = "Pup prod", xlab="years", ylab="n", ylim=c(0,2000000), type="b")
  points(as.vector(sealmod.params$pup.prod)[1:(length(1952:last.data.year)+1)]~seq(1952,(last.data.year+1)), col="red", pch=16)
  points(apply(pup,2,FUN=function(x){quantile(x, probs=0.025, na.rm=TRUE)})~seq(1952,(last.data.year+1)), type="l", col="blue")
  points(apply(pup,2,FUN=function(x){quantile(x, probs=0.975, na.rm=TRUE)})~seq(1952,(last.data.year+1)), type="l", col="blue")
  legend("topleft", c("Simulated pup prod total (median)",  "0.025 and 0.975 quantiles", "Observations"), lty = c(1,1,-1), pch = c(1,-1,16), col= c("black","blue","red"))
  
  ## Plot the simulated pup production for the 8+ age class and the other classes separatly
  plot(apply(PupProdByAge, c(2,1), median)[,5]~seq(1952,(last.data.year+1)), main="Pup prod by age", , xlab="years", ylab="n", ylim = c(min(PupProdByAge, na.rm=TRUE), max(PupProdByAge, na.rm=TRUE)), type="b", col=1)
  points(apply(PupProdByAge, c(2,1), FUN=function(x){quantile(x, probs=0.025, na.rm=TRUE)})[,5]~seq(1952,(last.data.year+1)), type="l", col="blue")
  points(apply(PupProdByAge, c(2,1), FUN=function(x){quantile(x, probs=0.975, na.rm=TRUE)})[,5]~seq(1952,(last.data.year+1)),type="l", col="blue")
  points(apply(PupProdByAge, c(2,1), median)[,6]~seq(1952,(last.data.year+1)), type="b", col=2)
  points(apply(PupProdByAge, c(2,1), FUN=function(x){quantile(x, probs=0.025, na.rm=TRUE)})[,6]~seq(1952,(last.data.year+1)), type="l", col="blue")
  points(apply(PupProdByAge, c(2,1), FUN=function(x){quantile(x, probs=0.975, na.rm=TRUE)})[,6]~seq(1952,(last.data.year+1)),type="l", col="blue")
  legend("topleft", c("Pup prod 8+ (median)", "Pup prod others (median)", "0.025 and 0.975 quantiles"), lty = c(1,1,1), pch = c(1,1,-1), col=c("black", "red", "blue"))
  
  ## Plot the simulated total population
  plot(apply(tot,2,median)~seq(1952,(last.data.year+1)), main = "Tot pop", xlab="years", ylab="n", ylim=c(1e+6, 1e+7), type="b")
  points(apply(tot,2,FUN=function(x){quantile(x, probs=0.025, na.rm=TRUE)})~seq(1952,(last.data.year+1)), type="l", col="blue")
  points(apply(tot,2,FUN=function(x){quantile(x, probs=0.975, na.rm=TRUE)})~seq(1952,(last.data.year+1)), type="l", col="blue")
  legend("topleft", c("Simulated population total (median)",  "0.025 and 0.975 quantiles"), lty = c(1,1), pch = c(1,-1), col= c("black","blue"))
  
  ## Save the plots as a pdf file
  dev.copy2pdf (file=paste(outdir,"Fitting_output_for_",nbrep,"_iterations_InitParams_", paste(params, collapse="_"), "_PregCor_", preg.cor, ".pdf", sep=""))
  #return(list(MakeCluster=t1,Calculation=t2))

  if (ShowAgeStruct[2]){

    # Age structure summary 2 (proportion of the total number in each age class)
    AgeStructPropSummary <- list(prob.025 = apply(AgeStructProp, c(1,2), function(x) quantile(x, prob = 0.025)), prob.25 = apply(AgeStructProp, c(1,2), function(x) quantile(x, prob = 0.25)), med = apply(AgeStructProp, c(1,2), median),
                   prob.50 = apply(AgeStructProp, c(1,2), function(x) quantile(x, prob = 0.50)), prob.75 = apply(AgeStructProp, c(1,2), function(x) quantile(x, prob = 0.75)), prob.975 = apply(AgeStructProp, c(1,2), function(x) quantile(x, prob = 0.975)))
    for(i in 1:length(AgeStructPropSummary)){ AgeStructPropSummary[[i]] <- data.frame(c(seq(0,14), "15+"), AgeStructPropSummary[[i]]) ; colnames(AgeStructPropSummary[[i]]) <- c("", seq(1952, last.data.year+1))}

    listPage <- vector()
    for ( i in 1:length(1952:(last.data.year + 1))){
      if(i == 1 | i%%20 == 0){
        windows(18,14)
        par(mfrow=c(4,5))
        listPage <- c(listPage, dev.cur())
      }
      barplot(AgeStructPropSummary$med[1:15,i+1], cex.names = 0.9, border="darkgrey", space=0, names.arg = AgeStructPropSummary$med[1:15,1], horiz = TRUE,las=1, main=colnames(AgeStructPropSummary$med)[i+1], xlim=c(0,0.25))
      mtext(paste("15+ =", round(AgeStructPropSummary$med[16,i+1],2)), side = 3, adj=-0.15, cex=0.7, col="darkgrey")
      mtext(paste("Age 1-14\n=", round(sum(AgeStructPropSummary$med[2:15,i+1]),2), "  "), side = 3, adj=0.8, line=-4, cex=0.7, col="darkgrey")
    }
    for (i in 1:length(listPage)){
      dev.set(which = listPage[i])
      dev.copy2pdf (file=paste(outdir,"AgeStructure_output_for_",nbrep,"_iterations_page", i,"-",length(listPage), ".pdf", sep=""))
    }
  }


  ## Final output for SealModel
  return(list("Params"=result,"Pupprod"=pup,"Totpop"=tot,"LastPopVect"=last.pop.vec, "Preg8plus"=preg8plus, "PregSample"=pregsampl, "PupSample"=pupSampl, "Preg8Sample" = preg8Sampl, "RunTime" = list("Make Cluster"=t1, "Parallel Calculations"=t2)))
  
}  
 




#############################################################################################
### Create cluster of computer processor threads for the calculation (Parallel computing) ###
#############################################################################################

createClust <- function(){
 library(parallel)
 
 if(exists("cltest")) {
    ans <- as.character(tkmessageBox(message = "A cluster already exists. Do you want to create a new cluster?\n(this will erase the information on the existing one)", icon = "question", type = "yesno", default = "no"))
    if (ans == "yes") stopCluster(cl = cltest)
 }else{
   ans <- "yes"
 }
 if (ans == "yes" | !exists("cltest")) {
      nbthreads <- ""
      cat("### Select the number of threads to be used in the cluster ###\n")
      while(nchar(nbthreads) == 0) nbthreads <- select.list(as.character(seq(1,detectCores(logical=TRUE)-1)), multiple=FALSE, title = "Select the number of threads")  # Choose the number of thread to be included in the cluster (keep one thread off to let windows run its jobs on this one)
      cltest <- makeCluster(as.numeric(nbthreads))
 }
 return(cltest)
}

## Run the function and print this run time

t3 <- system.time(cltest <- createClust())
cat("\n### External Cluster Creation Run Time ###\n")
print(t3)
flush.console()

###########################################################################################
###                             RUN THE SIMULATION                                      ###
###########################################################################################

## Run the function and print this run time
cat("\n### Simulation model Run Time ###\n")
print(system.time(resultSim <- SealModel(nbrep = 10000, params= c(.15,.4, 12000000), last.data.year=2019, preg.cor=0.85, NoConvLimit = 2, clust=cltest, ShowAgeStruct = c(TRUE, FALSE)))) # Use cluster created on the previous lines. Just rerun this line to make another run of the model

## Show internal run time
cat("\n### Internal SealModel Run Time ###\n")
print(resultSim$RunTime)

alarm()



#########################################################
## Details on the arguments of the SealModel function ###
#########################################################

## nbrep
# A numerical value corresponding to the number of model iterations (repetition) used to evaluate the parameters

## params
# A vector of three element corresponding respectively starting values for Alpha (the multiplying factor controlling the initial population size), the adult mortality (M), the carrying capacity (K)

## last.data.year
# A numerical value corresponding to the last year of data provided

## preg.cor
# A numerical value between 0 and 1 indicating the percentage of correlation expected between the pregnancy rates of seals from the different age classes (i.e. high value will imply that a good/bad year is good/bad for all age classes).

## NoConvLimit
# A numerical value corresponding to the number of model iterations allowed to let the model retry the fitting in case of non-convergence.

## clust ##
# The parameter clust define the number of threads used for the calculation
#   1- If it is not defined, the function will assume that it will use the number of threads available on that computer minus one (keep one thread off to let windows run its jobs on this one ... but maybe not needed)
#   2- If it is a number, the number will define the number of thread to be used (Note: if you have a dual core computer, windows create 4 threads)
#   3- You can also give to the function an object of class "cluster" created by the function makeCluster (package parallel) that define the threads already clustered and ready to run
#
#   Note: In the two first cases, the function SealModel will create the cluster before running the calculation. On windows, the time spend on this part increase exponentially with the number of thread included.
#         Moreover, the cluster created is stopped at the end of the function, so it is not usable for another calculation.
#
#         In the third case, you give a cluster to the function. So the function start parallel calculations immediatly (so it is the fastest way).
#         However, in order to obtain the cluster, you should use the function makeCluster before the first use of the SealModel function and thus it will spend some time creating the cluster (increasing exponentially with the number of threads included).
#         The advantage is this method is that this cluster will not be stopped at the end of the calculation, so you can re-use it for other runs.
#
#         Here is how to create that cluster:
#
#         library(parallel)
#         cltest <- makeCluster(3) # for example here, we create a cluster using 3 threads
#         SealModel(nbrep = 100, params= c(.2,.04, 9000000), last.data.year=2013, preg.cor=0.85, NoConvLimit = 2, clust = cltest)   # Here, we give the object created by the function makeCluster (i.e. cltest) to the SealModel function. 
#
#   Note2: The function will create outputs with a name specifying the parameters used, so on big computers you can run several R sessions at the same time running each a model with different parameters (you have to separate the number of threads available between the sessions (e.g. if you use the computer with 32 threads, you can run two models at the same time using 16 threads each)).
#                to find out how many threads available run the 2 lines below  
#               library(parallel)
#               detectCores(logical=TRUE)

## ShowAgeStruct ##
# The argument called ShowAgeStruct allows to extract and display the evolution of the population age structure
# (1) Number of individuals in each age classes for each year simulated is available in a text file created in the output directory
# (2) Plots proportion of the total number occurring in each age class for each year simulated. Plots are saved in the output directory.
# Default c(TRUE, TRUE) implies both (1) and (2). If you do not want data to be plotted just use c(TRUE, FALSE) for this argument.


# print(all<-system.time(resultSim <- SealModel(nbrep = 5000h, params= c(.2,.04, 9000000), last.data.year=2013, preg.cor=0.85, NoConvLimit = 2, clust=3)))


#######################################################
### Details on the output of the SealModel function ###
#######################################################

# The output of SealModel function can be used by the proj.sealmod.sim.f function (see below) to calculate projections of the population trajectory into the future

# Its a list of nine elements
# - Params: a data frame including the values obtained for Alpha, M, K, in each iteration of the model, along with the final value of the objective function issued from the fitting (fitting with the optim function)
# - Pupprod: a data frame including the pup production predicted for each modelled year in each iteration of the model
# - Totpop: a data frame including the total population size predicted for each modelled year in each iteration of the model
# - LastPopVect: a data frame including the number of individuals in each age class for the last year modeled in each iteration
# - Preg8plus: a data frame including the predicted pregnancy rate for the 8+ age class in each iteration
# - PregSample: a list providing for each iteration (each element of the list) a data frame including pregnancy rates resampled from observed values or smoothed ones for each age class in each modeled year
# - PupSample: a data frame including the pup production resampled from observed values for each modelled year in each iteration of the model
# - Preg8Sample: a data frame including the pregnancy rate for the 8+ age class resampled from observed values for each modeled year in each iteration
# - RunTime: a list of two elements including respectively the time spend in making the cluster (equal to 0 if a cluster already set is given to the function), and the time spend in parallel calculations (Note: this does not include the time spent plotting and saving the results)


################################################################################
#### Check Residual distribution around fitted values (Asked by a reviewer) ####
################################################################################
dev.new(); hist(as.vector(resultSim$Pupprod) - as.vector(resultSim$PupSample), 50, main="Residuals distribution around estimates of pup number") # Histogram of predicted vs observed values (note: observed values are resampled at each repetition)
PupSample <- resultSim$PupSample  
PupSample[,c(27:29,32)]<-NA
dev.new(); hist(as.vector(resultSim$Pupprod) - as.vector(PupSample), 50, main="Residuals distribution around estimates of pup number\nwithout pup estimates of 1978-1980 and 1983") # Histogram of predicted vs observed values (note: observed values are resampled at each repetition)
dev.new(); hist(as.vector(resultSim$Preg8plus) - as.vector(resultSim$Preg8Sample), 50, main="Residuals distribution around preg8+ estimates")




##############################################################################################
###---------------------------------------------------------------------------------------####
###                           PROJECTION IN THE FUTURE                                    ####
###---------------------------------------------------------------------------------------####
##############################################################################################


#   This function uses the result object created by the param.eval.rand function to project the population into the future
#   It uses the last year population vector as an initial population vector and its associated parameters to create the new population matrix from "last.year.data + 1" to "end simulation year
#   It projects by calculating mortality from quotas and others including variability in the mortality estimates
#   It then appened the pup production and total population vectors produced by the projection to the saved one (from 1952-2009)
#   It also creates graph to visualised the results and calculates basic statistics
#   The results are 2 vectors (final.sim.pup and final.sim.total) saved for each simulation
#   For more stability, always produce more or at least the same amount of parameters from the param.eval.rand function as the nb.proj requested in the projection function


#save.image("N:\\Arno\\ForMike\\Harp_seals\\Temp\\TestProj.RData")
#save.image("C:/Users/hammillm/Documents/Datas_seals/Harp/Populations/2015 Assessment/Simulations/PGFit1000030sept.RData")    
#rm(list = ls()); load("C:/Users/hammillm/Documents/Datas_seals/Harp/Populations/2015 Assessment/Simulations/pgtest.RData")
 
 


##################################################################
#### Define the directory where the files are read and written ###
##################################################################

  outdir<-paste(rep, "/modelAMKcor_output/", sep="") # set name of the output folder
  if (!file.exists(substring(outdir,1, nchar(outdir)-1))) {dir.create(outdir)}  #create directory if it does not exist




#########################################################################
#########################################################################
#########################################################################
#########################################################################
####                PROJECTION FUNCTION                        ##########
#########################################################################
#########################################################################
#########################################################################
#########################################################################

 
  proj.sealmod.sim.f <- function(result, last.data.year, end.year.sim, nb.proj, quota, preg.cor=0.85, fish=F){  # The result of the function is the pup production and total population, quota is a vector of quota which should be at least equal or longer than the number of projected years. You have to chose your correlation in pregnancy rates within years (the default=0.85).  Write fish=T if you want to save the projected population matrix
  
    ## Load the parameters
#    preg.rate.se <- sealmod.params$preg.rate.se
    preg.rate <- sealmod.params$preg.rate
    preg.size <- sealmod.params$preg.size
    pup.prod <- sealmod.params$pup.prod
    pup.prod.se <- sealmod.params$pup.prod.se
  
    ## Create the objects in which results will be stored
    final.sim.pup <- matrix(NA, nb.proj, end.year.sim-1952+1)
    final.sim.total <- matrix(NA, nb.proj, end.year.sim-1952+1)
    if (fish==T) {
      save.fish<-matrix(NA,26,nb.proj)                     #Initialize the matrix to save the population mean for the fish guy
      save.fish.sd<-matrix(NA,26,nb.proj)                  #Initialize the matrix to save the population standard deviation for the fish guy
    }                  
    library(MASS)
  
  
  
    ## Greenland catch model
    
     ## Greenland CASE 1 - PART A: You don't want to use a Greenland catch model. Uncomment the section for "Greenland CASE 1 - Part B" and Check that all the other CASES lines are commented
  
#    ## Greenland CASE 2 - PART A: Reevaluate the relation between Greenland catch and harp seal population size (based on new data for Greenland catch and/or new evaluation of harp seal population)
#    ##                  Need also to uncomment the section "Greenland CASE 2 - PART B" to be taken into account in the projection (check that all the other CASES lines are commented)
#      # merge info on greenland catch (provided in the raw-removal-1952.csv file) and seal pop size estimated from the simulation part (as mean size of the values estimated for each year in each iteration of the simulaiton model_
#      greencatch <- data.frame(catch = raw.remov[,3], pop=apply(result$Totpop, 2,mean)[1:length(1952:last.data.year)])  # if you want to run this outside of the function, replace "result" by the name of the object that contains the results of the simulations (e.g. resultSim). 
#
#      # Estimate breakpoint and first linear curve parameters
#      library(SiZer)
#      model <- piecewise.linear(greencatch$pop,greencatch$catch)
#      predFirstRange <- data.frame(pop=seq(min(greencatch$pop), as.numeric(model[1]), by=10000))
#      lmFirstRange <- lm(catch ~ pop,  data = greencatch[greencatch$pop <= model[1],])
#      pred <- predict(lmFirstRange, newdata = predFirstRange, interval="prediction")
#
#      # Estimate second linear curve parameters
#      #sdSecondRange <- sd(greencatch[greencatch$pop > model[1],]$catch)  # uncomment this line and the two following lines if you want to add those values to the graph below (these values are not used in the calculations)
#      #lwCI <- pred[nrow(pred),1] -  1.96 * sdSecondRange / sqrt(nrow(greencatch[greencatch$pop > model[1],]))
#      #upCI <- pred[nrow(pred),1] +  1.96 * sdSecondRange / sqrt(nrow(greencatch[greencatch$pop > model[1],]))
#      predMean <- as.numeric(predict(lmFirstRange, newdata = data.frame(pop= as.numeric(model[1]))))
#      predUp <- predMean + (max(greencatch[greencatch$pop > model[1],]$catch) - min(greencatch[greencatch$pop > model[1],]$catch)) / 2    # mean + (real range for Greenland catch when pop>breakpoint) / 2
#      predDown <- predMean - (max(greencatch[greencatch$pop > model[1],]$catch) - min(greencatch[greencatch$pop > model[1],]$catch)) / 2  # mean - (real range for Greenland catch when pop>breakpoint) / 2
#  
#      # Plot   # uncomment the following lines only if you want a graph of the relationship (piecewise regression) between greenland catches and the harp seal population size
#      plot(pred[,1] ~ predFirstRange$pop, type="l", xlim = range(greencatch$pop), ylim = range(greencatch$catch), xlab = "Pop", ylab = "Catch")
#      points(pred[,2] ~ predFirstRange$pop, type="l", col="red")
#      points(pred[,3] ~ predFirstRange$pop, type="l", col="red")
#      segments(x0=as.numeric(model[1]), y0=predMean, x1= max(greencatch$pop), y1 = predMean)
#      #segments(x0=as.numeric(model[1]), y0=lwCI, x1= max(greencatch$pop), y1 = lwCI, col="red")  #
#      #segments(x0=as.numeric(model[1]), y0=upCI, x1= max(greencatch$pop), y1 = upCI, col="red")
#      segments(x0=as.numeric(model[1]), y0=predUp, x1= max(greencatch$pop), y1 = predUp, col="red")
#      segments(x0=as.numeric(model[1]), y0=predDown, x1= max(greencatch$pop), y1 = predDown, col="red")
#      points(catch~pop, data = greencatch, col="blue")
  
  
#     ## Greenland CASE 3 - PART A: Create a function (GreencatchSample)that provide an estimate of the Greenland catch based on the harp seal population size.
#     ##                  Need also to uncomment the section "Greenland CASE 3 - PART B" to be taken into account in the projection (check that all the other CASES lines are commented)
#     ##                  Relationships used in the function are based on greenland catch data 1952-2012
#     ##                  and population size estimated by 1000 iterations of the simulation model using binomials for preg rates estimations and without effect of K on 8+ (only effect on age1; K set to 11024440)
#     ##                  The function give results based on a piecewise regression.
#     ##                  In the first linear part (i.e. before the break point), the function provide a random sample value extracted within the 95% prediction interval around the mean catch estimated by the piecewise regression
#     ##                  In the second linear part, the function provide a random sample value extracted within the interval defined as mean value estimated at the breakpoint + or - half the real range of Greenland catch values when the population > breakpoint
#
##      # Uncomment the following "Explanation" section only if you want to learn how the values used in the GreencatchSample function were calculated
##        ########### EXPLANATION ##############
##        # Data
##        catch <- c(16400, 16400, 19150, 15534, 10973, 12884, 16885, 8928, 16154, 11996, 8500, 10111, 9203, 9289, 7057, 4242, 7116, 6438, 6269, 5572, 5994, 9212, 7145, 6751.5, 11956, 12866, 16638, 17544.5,
##                  15255, 22973.5, 26926.5, 24784.5, 25828.5, 20785, 26098.5, 37859, 40414.75, 42970.5, 45526.25, 48082, 50637.75, 56319, 57373, 62749,73947, 68815.5, 81272, 93117, 98458.5, 85427.5, 66734.5, 66149,
##                  70585.5, 91695.5, 92210, 82778, 78897, 70680, 82843, 82843)
##        pop <- c(2307425, 2280462, 2296843, 2281399, 2202525, 2128308, 2161852, 2048831, 1987451, 1879289, 1960551, 1847559, 1746831, 1642978, 1644785, 1570090, 1484620, 1523909, 1374638, 1351427, 1376430,
##                 1498327, 1597074, 1671225, 1736373, 1825502, 1947993, 2113088, 2199440, 2332988, 2324212, 2484364, 2764362, 3036309, 3478664, 3751776, 4031376, 4345599, 4543594, 4838039, 5257876, 5314615,
##                 5896304, 6047371, 6346524, 6542766, 6594865, 6600070, 6961693, 7123091, 7244607, 7285746, 7036484, 7313658, 7219079, 7643993, 7789584, 7729284, 7409989, 7022857)
##        greencatch <- data.frame(catch = catch, pop=pop)
##  
##        # Piecewise Model (catch ~ -1.400414e+04 + 1.361350e-02 * population size; Breakpoint at pop = 7083436 and catch = 82426)
##        model <- piecewise.linear(greencatch$pop,greencatch$catch)
##        
##        # First linear part (pop <= 7083436):
##         # Calculate prediction interval (95% of predicted values are included in this interval)
##         ngreencatch <- nrow(greencatch[greencatch$pop <= model[1],]) # 51
##         MRSS <- 1/(ngreencatch - 2) * sum(residuals(lmFirstRange)^2) # 28447011 #mean residual sum of square
##         meanpop <- mean(greencatch$pop[greencatch$pop <= model[1]]) # 3161567
##         SSpop <- sum((greencatch$pop[greencatch$pop <= model[1]] - mean(greencatch$pop[greencatch$pop <= model[1]]))^2) # 1.762231e+14
##         dif <-  qt(0.025,(51-2))*sqrt(28447011*(1+1/51+((pop-3161567)^2/1.762231e+14)))
##         # example: at the break point
##         pop <- 7083436
##         dif <-  qt(0.025,(51-2))*sqrt(28447011*(1+1/51+((pop-3161567)^2/1.762231e+14)))
##         (-1.400414e+04 + 1.361350e-02 * pop) + abs(dif)  # upper prediction
##         (-1.400414e+04 + 1.361350e-02 * pop) - abs(dif)  # lower prediction
##        
##        # Second linear part (pop > 7083436): (mean Greenland catch = 82426; real range for Greenland catch when pop>breakpoint = 26061 
##         # example
##         pop <- 7100000
##         (-1.400414e+04 + 1.361350e-02 * pop) + diff(range(greencatch$catch[greencatch$pop > model[1]]))/2  # upper prediction
##         (-1.400414e+04 + 1.361350e-02 * pop) - diff(range(greencatch$catch[greencatch$pop > model[1]]))/2  # lower prediction
##        ########### END EXPLANATION ############## 
#
#  
#       # GreencatchSample function (do not forget to uncomment the section "Greenland CASE 3 - part B") 
#        GreencatchSample <- function(pop){     # This function allows to sample Greenland catch values into distributions detailed above
#          if (pop <= 7083436){
#            meancatch <- -1.400414e+04 + 1.361350e-02 * pop
#            dif <-  qt(0.025,(51-2))*sqrt(28447011*(1+1/51+((pop-3161567)^2/1.762231e+14))) # 95% prediction interval = meancatch + or - dif
#            samplecatch <- rnorm(1, mean = meancatch, sd = (abs(dif) / 1.96))
#          }else{
#            interval <- c(82426 - 26061 / 2, 82426 + 26061 / 2)
#            samplecatch <- runif(1, interval[1], interval[2])
#          }
#          samplecatch[samplecatch < 0] <- 0
#          return(round(samplecatch))
#        }



    ########### PROJECTION LOOP #################
  
    ## Progress bar
    maxbar1 <- nb.proj
    pb1 <- tkProgressBar(title = "Projections calculation progress", min = 0, max = maxbar1, width = 400)
  
    for (j in 1:nb.proj){
  
      setTkProgressBar(pb1, j, label=paste("Percent done",round(j/maxbar1*100), "%")) # update progress bar
  


      ## PARAMETERS
       ## There are 2 ways of selecting M and K for the projection in the future.  They are both based on the result object created by the function SealModel
#      ## 1. M and K are randomly sampled from normal distributions with a mean equal to the mean predicted values in the simulation the corresponding sd
#       M <- rnorm(1,mean(result[[1]][,2]),sd(result[[1]][,2]))
#       K <- rnorm(1,mean(result[[1]][,3]),sd(result[[1]][,3]))

       ## 2. M and K are selected as a pair from the result of the SealModel function.
        # To use this, you need to have at least the same amount of parameter pairs as the number of projection you are looking for  
        if (nrow(result[[1]]) < nb.proj){ stop("You need to have at least the same number of parameter combinations as the number of projection you are looking for")} # check if TRUE
        M <- result[[1]][j,2]
        K <- result[[1]][j,3]
        
      ## INITIAL POPULATION VECTOR
      init.sim <- proj.seal <- result[[4]][,j] # this is the population the year after the last data year. Note: the vector is matched with the parameters when the parameters (M, K) are seleted as pairs (second choice above)
      year.proj<-end.year.sim - last.data.year - 1
  
      ## SETUP THE OUTPUT MATRIX
      proj.seal.sim <- matrix(ncol= year.proj+1, nrow= 26)
      proj.seal.sim[,1] <- init.sim
  
      ## NEW REMOVAL VECTOR (includes all the catch from canada, greenland, arctic and bycatch)
      prop.pup.killed <- c(0.5,0.14,0.03,0.6)      #In order: canada, greenland, Arctic, bycatch,   we could easily add variance if needed, canada average of last 10 years
      SandLpup<-c(1/0.95,1/0.5,1/0.5,1)            #Vector of S&L for pup, values are for : (in order) canada, greenland, arctic, bycatch.  I could add variance on that as well
      SandLadult<-c(1/0.5,1/0.5,1/0.5,1)           #Vector of S&L for adult, values are for : (in order) canada, greenland, arctic, bycatch.  I could add variance on that as well
  
      
      ## ICE MORTALITY VS ICE COVERAGE  ### TO BE SIMPLIFIED

     #medianIce <- c(rep(8.5,5), 7, 4.5, 2.5, 1, 0.2, 0)  # approximate data from the September sea ice extent
     medianIce <- c(rep(3,4),3, 2.9,2.9,2.8,2.7,2.6,2.4, 2.2, 2.0, 1.9, 1.7, 1.5,1.4,1.3,1.1, 1.0, 0)  # approximate data from the paper winter ice extent  climate model 4.5
     
# medianIce <- c(rep(3,4),30, 28,26,24,22,20,18.5, 18, 16, 14, 12, 10,8,6,4, 2, 0)  # climate model 8.5
     
   #  medianIce <- c(rep(20.2,10), 20, 18.8, 19.1, 17.4, 10.2, 10,8.7,6,4, 3, 1)  # made up data
     
     set.seed = 205        # arnaud had a seed of 111, 165
     medianIce2 <- medianIce + rnorm(length(medianIce), 0, .1)  # Add noise
     set.seed = 155           #Arnaud had 223 140
#     medianIce3 <- medianIce + rnorm(length(medianIce), 0, 1.5)  # Add noise
     medianIce3 <- medianIce - rnorm(length(medianIce), 0, .1)  # Add noise
     dat <- data.frame(year = rep(seq(1900, 2100, by = 10), 3), ice=c(medianIce, medianIce2, medianIce3))       # arnaud did 100 y

     if(j == 1)  {
     
     modIce <- nls(ice ~ SSlogis(year, Asym, xmid, scal), data = dat)  # fit a logistic model on the data
     
      nls.control(maxiter = 100, tol = 1e-05, minFactor = 1/3000, printEval = FALSE, warnOnly = FALSE)         #1/1024
      }
#     plot(predict(modIce, newdata= data.frame(year = seq(1900, 2100))) ~ seq(1900, 2100) )

     # Extract parameters and add random normal noise with SD = to the SD estimated by the logistic model
     # This will allow to obtain a different logistic curve from each projection
     sel1 <- rnorm(1, coef(modIce)[1], summary(modIce)$parameters[, "Std. Error"][1])
     sel2 <- rnorm(1, coef(modIce)[2], summary(modIce)$parameters[, "Std. Error"][2])
     sel3 <- rnorm(1, coef(modIce)[3], summary(modIce)$parameters[, "Std. Error"][3])
     sel3 <- ifelse( sel3>0, -sel3, sel3)    # avoid positive scale that would reverse the logistic curve (i.e. increasing vs decreasing)

     # Add variability to the logistic curve predictions
     fun <- function(x, Asym, xmid, scal) { Asym / (1 + exp((xmid - x) / scal)) } # Logistic function
     pred <- fun(seq(last.data.year, end.year.sim) , sel1, sel2, sel3)
     set.seed = 200
     #set.seed = 333
     predNvar <- pred + rnorm(length(pred), 0, 0.3) # random normal noise SD estimated from the graph
     predNvar <- ifelse(predNvar < 0, 0, predNvar)
#      plot( fun(seq(1900, 2100) , sel1, sel2, sel3) + rnorm(length(seq(1900, 2100)), 0, 0.3) ~ seq(1900, 2100), type="l")
#     plot(predNvar ~ seq(last.data.year, end.year.sim))
     # iceProp <- predNvar / 8.5     # i replaced aranud's value of 12 with 8.5 for cosewic'
     iceProp <- predNvar / 3
  #   print(iceProp);flush.console()




      ## LOOP BY YEAR
      
      for (i in 1:year.proj) {           #Note that the loop goes 1 after the last year simulated
  
        ## GREENLAND REMOVALS
         # Greenland CASE 1 - Part B
          green.remov <-  runif(1,45000,55000)

#         # Greenland CASE 2 - Part B
#          if (sum(proj.seal) < model[1]){
#            p <- predict (lmFirstRange, newdata = data.frame(pop = sum(proj.seal)), interval="prediction")
#            green.remov <- rnorm(1, mean = p[1], sd = (p[1] - p[2]) / 1.96)
#          }else{
#            green.remov <- runif(1, predDown, predUp)
#          }
#
#         # Greenland CASE 3 - Part B
#          green.remov <- GreencatchSample(sum(proj.seal)) # use the function GreencatchSample defined above
  
        ## CATCHES FROM OTHER PLACES
         arctic.remov <- runif(1,999,1001) 
         bycatch.remov <- runif(1,1000,3000)
        # bycatch.remov <- runif(1,12289,12291)
  
        ## TOTAL REMOVALS
        #minister<-c(1.0,1.0,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.1)
        minister<-1
         proj.remov <- rbind(quota[i]*(sample(minister,1,replace=TRUE)) ,green.remov,arctic.remov,bycatch.remov)      #Minister factor:if a uniform distribution use(runif(1,1.0,1.1))instead of sample
         prop.age.class<-prop.table(proj.seal.sim[2:26,i])				#This creates a table of age class proportion to evaluate the kills
  
         # The next 3 lines create a vector of number of seal killed per age class which will be used as the remove vector in the projection model
         adult.kill <- sum((1-prop.pup.killed)*proj.remov*SandLadult)*prop.age.class
         pup.kill <- sum(proj.remov*prop.pup.killed*SandLpup)
         kill <- append(pup.kill,adult.kill)
  
        ## EVENTS INDUCING MORTALITY OR ACTING ON PREGNANCY RATES
         # Ice factor
#         icefactor <- 1 # uncomment this one and comment the following if you want to ignore the ice factor
#         icefactor <- sample(c(1,1,.85,.86,.88,0.714466707,0.35109587,0.76642337,0.772912202),1,replace=TRUE) # define the proportion of pup surviving due to ice condition   new quantitative measure
#         icefactor <- sample(c(.94,0.59,0.21,0.9,1),1) # define the proportion of pup surviving due to ice condition  -old measure
         iceThreshold = .7
         icefactor <- ifelse(iceProp[i] > iceThreshold, 1, (iceProp[i]))  # the proportion of pup surviving will be 1 if the iceProp is over the iceThreshold (to be defined) and will be equal to the proportion of ice remaining if below the threshold
         print(icefactor); flush.console()
         # Food factor
#         test <- c(1.76,1.80,0.2,0.45,1.92,0.19,0.01,0.27,0.56,1.47,0.24,0.93,1.32,0.65,0.76,1.14,1.19,0.67,1.32,
#         0.05,1.45,1.02,1.38,1.43,1.43,0.12,0.72,1.67,1.71) # create the vector
 
          
         foodfactor <- 1 # uncomment this one and comment the following if you want to ignore the ice factor
# #        foodfactor <- sample(test), 1) # define good or bad year for food.
                                                  # Act in reducing the effective pregnancy rate in a bad year. Here at mean 1 over 5 years is bad (10% of the normal pregnancy rates).
                                                  # The best year (1 over 5) induces a 10% increase over normal pregnancy rates (110% of the normal pregnancy rates)
#         foodfactor<- rnorm(100,1,0.2)
#           foodfactor <- sample(test)
        ## APPLY MORTALITY
         # Seal numbers at age 1
         proj.seal.age1 <- (proj.seal[1] * icefactor - kill[1]) * exp(- M)  *(1-(sum(proj.seal.sim[,(i)])/K)^2.4) # Include icefactor (see "Ice factor" section) and density dependence
         proj.seal.age1[is.nan(proj.seal.age1)]<-0.001
         proj.seal.age1[proj.seal.age1<0]<-0.001
  
         # Seal numbers for age > 1 and < 25
         proj.seal.num <-  (proj.seal[2:24]*exp(-M1/2) - kill[2:24]) * exp(-M1/2)
         proj.seal.num[proj.seal.num<0]<-0.001
  
         # Seal numbers for age 25+
         #proj.seal.old <- (proj.seal[25] * exp(-M/2) - kill[25]) * exp(-M/2)
         proj.seal.old <- ((proj.seal[26] + proj.seal[25]) * exp(-M1/2) - kill[25] - kill[26]) * exp(-M1/2)
         proj.seal.old[proj.seal.old<0]<-0.001
  
        ## PREGNANCY RATES OPTIONS
         # OPTION 1
 #        proj.preg.rate <- preg.rate[,length(preg.rate[1,])]  			# Is the last vector of pregnancy rates, therefore the pregnancy rate for future years is believed stable and similar to the one from the last data year.  this could be changed easily.
 #        proj.preg.size <- preg.size[,length(preg.size[1,])]
  
         # OPTION 2:  
#         proj.preg.rate <- c(0,0,0,  0.05,0.1,0.2,0.2,   rep(0.7,19))         #rep(0.7,19) means that 0.7 is repeated 19 times. It represents the pregnancy rates for ages 8 to 26 (i.e., 8+).
#         proj.preg.size <- c(1,1,1,  5,5,5,5, rep(40,19))                     # sample size 40 is repeated 19 times for the 8+ these values can be changed
  
         # OPTION 3: this might be the rpd up to 2013 
         #proj.preg.rate <- c(0,0,0,  0.05,0.1,0.2,0.3,   rep(runif(1,.2,.7),19))         #rep(0.7,19) means that 0.7 is repeated 19 times. It represents the pregnancy rates for ages 8 to 26 (i.e., 8+).
#          proj.preg.rate <- c(0,0,0,  # age 1 to 3
#                              sample(c(0.02,0.02,0.02,0.02,.02,.02,.02,.03,.03,.03),1,replace=TRUE),  #age 4
#                             sample(c(0.01,0.01,0.01,0.0,0.0,0.0,0.0,0.0,0.0,0.0),1,replace=TRUE),  #age 5
#                              sample(c(0.67,0.25,0.07,0.06,.05,.05,.04,.3,.02,.02),1,replace=TRUE),  # age 6  last 10 years sample(c(0.33,0.24,0.23,0.23,.23,.33,.23,.33,.24,.24),1),
#                              sample(c(.39,.37,.36,.35,.33,.32,.31,.30,.29,.27),1,replace=TRUE),  # age 7  last 10 years sample(c(0.33,0.32,.32,0.31,.31,.75,.3,.3,.50,.25),1),
#                              rep(sample(c(0.3,0.19,0.42,0.54,.86,.79,.74,.58,.74,.73),1,replace=TRUE),19))    # age 8+   last 10 years
  
#          proj.preg.size <- c(1,1,1,  3,4,4,4, rep(5,19))  # sample size 5 is repeated 19 times for the 8+ these values can be changed
#          rand.proj.preg.rate <- rbinom(26,size=proj.preg.size,prob=proj.preg.rate)/proj.preg.size

# OPTION 3:    this is rpd for last 10 years from 2019
#         proj.preg.rate <- c(0,0,0,  0.05,0.1,0.2,0.3,   rep(runif(1,.2,.7),19)) #rep(0.7,19) means that 0.7 is repeated 19 times. It represents the pregnancy rates for ages 8 to 26 (i.e., 8+).
          proj.preg.rate <- c(0,0,0,  # age 1 to 3
           sample(c(0.02,0.02,0.02,0.02,.02,.02,.02,.03,.03,.03),1,replace=TRUE),  #age 4
           sample(c(0.01,0.01,0.01,0.0,0.0,0.0,0.0,0.0,0.0,0.0),1,replace=TRUE),  #age 5
           sample(c(0.24,0.24,0.24,0.23,.23,.23,.23,.23,.3,.24),1,replace=TRUE),  # age 6  last 10 years sample(c(0.33,0.24,0.23,0.23,.23,.33,.23,.33,.24,.24),1),
           sample(c(.34,.33,.33,.32,.31,.31,.31,.30,.3,.3),1,replace=TRUE),  # age 7  last 10 years sample(c(0.33,0.32,.32,0.31,.31,.75,.3,.3,.50,.25),1),
          rep(sample(c(0.31,0.20,0.42,0.54,.86,.79,.74,.58,.74,.73),1,replace=TRUE),19))    # age 8+   last 10 years
  
          proj.preg.size <- c(1,1,1,  3,4,4,4, rep(5,19))  # sample size 5 is repeated 19 times for the 8+ these values can be changed
          rand.proj.preg.rate <- rbinom(26,size=proj.preg.size,prob=proj.preg.rate)/proj.preg.size

  
        # OPTION 4: RANDOM SAMPLE IN BINOMIALS based on last.data.year pregnancy rates (need the two lines of OPTION 1 activated)
#         rand.proj.preg.rate <- rbinom(26,size=proj.preg.size,prob=proj.preg.rate)/proj.preg.size
  
        # OPTION 5: CORRELATED RANDOM SAMPLE IN BINOMIALS based on last.data.year pregnancy rates (need the two lines of OPTION 1 activated) but modifed for K effect for age 8+
#         # K effect on preg rates of 8+ classes
#         proj.preg.rate[8:26] <- rep(ifelse(0.88 * (1 - (sum(proj.seal.sim[,i])/K)^2.4) < 0, 0, (0.88 * (1 - (sum(proj.seal.sim[,i])/K)^2.4))),19) # 0.88 is considered as max preg rate (can be defined as max preg rate when pop is low)
  
#         # Correlated random sample in binomials (note that proj.preg.rate was modified by previous lines
#         norm.cop <- normalCopula(preg.cor, dim=8)
#         rcop <- rcopula(norm.cop, 1)
#         rand.proj.preg.rate <- sapply(seq(1,8), FUN = function(x){qbinom(rcop[,x],size=proj.preg.size[x], prob = proj.preg.rate[x])/proj.preg.size[x]})
#         rand.proj.preg.rate <- c(rand.proj.preg.rate, rand.proj.preg.rate[rep(8, 18)])
  
        # OPTION 6: K effect on preg.rates for 8+ and preg.rates for classes 4 to 7 based on their relation with preg.rate of classe 8
##         rand.proj.preg.rate <- rep(0,26)
  
#         # K effect on preg rates of 8+ classes
##        rand.proj.preg.rate[8:26] <- rep(ifelse(0.88 * (1 - (sum(proj.seal.sim[,i])/K)^2.4) < 0, 0, (0.88 * (1 - (sum(proj.seal.sim[,i])/K)^2.4))),19) 
         # 0.88 is considered as max preg rate (can be defined as max preg rate when pop is low)
  
#          # This relation was considered as a logistic function fitted on raw pregrates for classes 4-8
#          pregraw <- read.csv("temppr_2014.csv", header=F); colnames(pregraw) <- c("years", "age", "total", "npreg")
#          pregraw$pregrate <- pregraw$npreg / pregraw$total
#          pregraw <- na.omit(pregraw)
#          logispreg <- nls(pregrate~SSlogis(age,a,m,s), weights = total, data=pregraw)   # weighted by sample size
#
#          # Sub-option 1 
#           rand.proj.preg.rate[4:7]<-SSlogis(4:7, rand.proj.preg.rate[8], coef(logispreg)[2], coef(logispreg)[3])  # the first parameter (the asymptotic value) is replaced by the preg.rate of classe 8 in order to estimate values for classes 4-7
#          # Sub-option 2
#           rand.proj.preg.rate[4:7]<-SSlogis(4:7, rand.proj.preg.rate[8], 5.368324, 0.5674322) # application of the previous line but avoiding recalculation of the logistic curve for each run (Warning, coefficient values could have changed with new data)
#          # Sub-option 3
#           rand.proj.preg.rate[4:7]<-SSlogis(4:7, rand.proj.preg.rate[8], rnorm(1,coef(logispreg)[2], summary(logispreg)$coefficient[2,2]), rnorm(1,coef(logispreg)[3], summary(logispreg)$coefficient[3,2]))  # same as previous line but allows variability for paramaters 2 and 3
#          # Sub-option 4
##           rand.proj.preg.rate[4:7]<-SSlogis(4:7, rand.proj.preg.rate[8], rnorm(1,mean = 5.368324, 0.13783), rnorm(1,mean = 0.5674322, 0.12531)) # application of the previous line but avoiding recalculation of the logistic curve for each run turn this on if running density dependence preg rates (Warning, coefficient values could have changed with new data)
  
        ## END PREGNANCY RATES OPTIONS
  
        
        ## APPLY NATALITY
         # Apply effect on natural event on pregnancy rates
         rand.proj.preg.rate <- ifelse(rand.proj.preg.rate > 0.88, 0.88,rand.proj.preg.rate) # no ice or food factor operating on rpd rate or pups born but lost immediately
  
         
#         rand.proj.preg.rate <- ifelse(icefactor *foodfactor* rand.proj.preg.rate > 0.88, 0.88, rand.proj.preg.rate*icefactor*foodfactor)
  
         # Calculate pups
         proj.seal.pups <- sum(rand.proj.preg.rate[2:26] * c(proj.seal.age1, proj.seal.num, proj.seal.old)/2)
         proj.seal.pups[proj.seal.pups < 0.001] <- 0.001
  
        ## NEW VECTOR OF SEALS
         proj.seal[1] <- proj.seal.pups
         proj.seal[2] <- proj.seal.age1
         proj.seal[3:25] <- proj.seal.num
         proj.seal[26] <- proj.seal.old
  
         proj.seal[!is.finite(proj.seal)] <- 0.001
         proj.seal[proj.seal < 0.001] <- 0.001
  
         proj.seal.sim[,(i+1)] <- proj.seal
         proj.seal.sim
      }
  
      if (fish==T) {
        data.fish<-matrix(NA,26,2,dimnames=list(c(1:26),c("mean","Stand Dev")))    ##########  Code to save mean and SD of seal population for the projected year for the fish guy
        save.fish[,j]<-apply(proj.seal.sim[,2:(year.proj+1)],1,mean)               ##  To activate, write "fish=T" while initializing the function
        save.fish.sd[,j]<-apply(proj.seal.sim[,2:(year.proj+1)],1,sd)              ##  The data are save in data.fish
        data.fish[,1]<-apply(save.fish,1,mean)                                     ##
        data.fish[,2]<-apply(save.fish.sd,1,mean)                                  ##
        data.fish<<-data.fish                                                      ##
      }                                                                            ##########  end of the fish code
  
      projection.result.pup <- append(result[[2]][j,],proj.seal.sim[1,2:(year.proj+1)])
      projection.result.total <- append(result[[3]][j,],apply(proj.seal.sim[,2:(year.proj+1)],2,sum))
      final.sim.pup[j,] <- projection.result.pup
      final.sim.total[j,] <- projection.result.total
    }
  
    close(pb1)
  
    ## Estimate the quantiles from all the trajectories
  
     # for the pups
     final.sim.pup.0.5 <- apply(final.sim.pup,2,quantile,0.5)
     final.sim.pup.0.025 <- apply(final.sim.pup,2,quantile,0.025)
     final.sim.pup.0.975 <- apply(final.sim.pup,2,quantile,0.975)
     final.sim.pup.sd <- apply(final.sim.pup,2,sd)
     final.sim.pup.mean <- apply(final.sim.pup,2,mean)
  
     # for the total
     final.sim.total.0.5 <- apply(final.sim.total,2,quantile,0.5)
     final.sim.total.0.025 <- apply(final.sim.total,2,quantile,0.025)
     final.sim.total.05 <- apply(final.sim.total,2,quantile,0.05)
     final.sim.total.10 <- apply(final.sim.total,2,quantile,0.1)
     final.sim.total.20 <- apply(final.sim.total,2,quantile,0.2)
     final.sim.total.25 <- apply(final.sim.total,2,quantile,0.25)
     final.sim.total.30 <- apply(final.sim.total,2,quantile,0.3)
     final.sim.total.40 <- apply(final.sim.total,2,quantile,0.4)
     final.sim.total.95 <- apply(final.sim.total,2,quantile,0.95)
     final.sim.total.0.975 <- apply(final.sim.total,2,quantile,0.975)
     final.sim.total.sd <- apply(final.sim.total,2,sd)
     final.sim.total.mean <- apply(final.sim.total,2,mean)
  
     test.table<-cbind(final.sim.pup.0.5,final.sim.pup.mean,final.sim.pup.sd, final.sim.pup.0.025,final.sim.pup.0.975,
     final.sim.total.0.5,final.sim.total.mean, final.sim.total.sd, final.sim.total.0.025, final.sim.total.05,final.sim.total.10, final.sim.total.20,
     final.sim.total.30, final.sim.total.40, final.sim.total.95, final.sim.total.0.975)
     write.csv(test.table, file= paste(outdir, "test_reference.csv", sep=""))
  
    ## Plot the pup distribution,
     dev.new()
     par(mfrow = c(2,1), mar = c(4.5,5,1,1))
  
     plot(1952:end.year.sim, pup.prod[1:(end.year.sim-1952+1)]/1e6, ylim=c(0,signif(max(final.sim.pup.0.975),2)/1e6),xlab="Years",ylab="Pup prod estimates (x 1e+06)",xaxp=c(1950,2100,30), cex.lab = 1.2, cex.axis = 1.2, pch = 19)
     polygon(x = c(1952:end.year.sim, rev(1952:end.year.sim)), y = c(final.sim.pup.0.025/1e6,rev(final.sim.pup.0.975/1e6)), col = "darkgrey", border = NA)
     abline(h = seq(0,signif(max(final.sim.total)/1e6,2),by = 0.1), col = "gray80") #plots the grey horizontal lines
     points(1952:end.year.sim, pup.prod[1:(end.year.sim-1952+1)]/1e6, pch = 19)
  
     for (n in 1952:end.year.sim) {
       segments(x0 = n, y0 = (pup.prod[n-1952+1]-pup.prod.se[n-1952+1])/1e6, x1 = n, y1 = (pup.prod[n-1952+1]+pup.prod.se[n-1952+1])/1e6)
       segments(x0 = n-0.3, y0 = (pup.prod[n-1952+1]+pup.prod.se[n-1952+1])/1e6, x1 = n+0.3, (pup.prod[n-1952+1]+pup.prod.se[n-1952+1])/1e6)
       segments(x0 = n-0.3, y0 = (pup.prod[n-1952+1]-pup.prod.se[n-1952+1])/1e6, x1 = n+0.3, (pup.prod[n-1952+1]-pup.prod.se[n-1952+1])/1e6)
     }
     lines(1952:end.year.sim,final.sim.pup.0.5/1e6)
  
  
    ## Plot the total number of seal
  
     plot(1952:end.year.sim, final.sim.total.0.5, ylim=c(0,signif(max(final.sim.total.0.975)/1e6,2)),xlab="Years",ylab="Seal numbers (x 1e+06)", xaxp=c(1950,2100,30), cex.lab = 1.2, cex.axis = 1.2, type = "l")
     polygon(x = c(1952:end.year.sim, rev(1952:end.year.sim)), y = c(final.sim.total.0.025/1e6,rev(final.sim.total.0.975)/1e6), col = "darkgrey", border = NA)
     abline(h = seq(0,signif(max(final.sim.total)/1e6,2),by = 1), col = "gray80") #plots the grey horizontal lines
     lines(1952:end.year.sim,final.sim.total.0.5/1e6)
  
     par(mfrow = c(1,1))
     dev.copy2pdf (file=paste(outdir,"Projection_output_for_",nb.proj,"_iterations.pdf", sep=""))
  
    ## This part of the code just gives out basic statistics about the distribution of the simulations
  
     pop.final <- final.sim.total[,(last.data.year-1952+2):(end.year.sim-1952+1)]
  
     Med <- ceiling(final.sim.total.0.5[(last.data.year-1952+2):(end.year.sim-1952+1)])
     N70 <- 0.70*max(final.sim.total.0.5[1:(last.data.year-1952+1)])
     N50 <- 0.50*max(final.sim.total.0.5[1:(last.data.year-1952+1)])
     N30 <- 0.30*max(final.sim.total.0.5[1:(last.data.year-1952+1)])
     respect.N70 <- apply(pop.final>N70,2,sum)/nb.proj
     respect.N50 <- apply(pop.final>N50,2,sum)/nb.proj
     respect.N30 <- apply(pop.final>N30,2,sum)/nb.proj
  
     respect.table <- matrix(c(Med,respect.N70,respect.N50,respect.N30),4,length((last.data.year+1):end.year.sim),byrow=TRUE,dimnames=list(c("Pop","N70","N50","N30"),c((last.data.year+1):end.year.sim)))
  
     ## Save the results of the simulation in a object in R
     final.sim.pup <<- final.sim.pup
     final.sim.total <<- final.sim.total
  
     return(respect.table)
  }


#################################################################################################
###                         RUN THE PROJECTION FUNCTION                                       ###
#################################################################################################

  ## Indicate scenario of removal -you can have the same harvest every year, or uncomment and have different harvest for each year of future simulation
    quota<-list()
#    quota[["60k per year"]] <- rep(60000, 21)
#    quota[["125k per year"]] <- rep(125000, 21)
#    quota[["150k per year"]] <- rep(150000, 21) # one for each projection year
    quota[["175k per year"]] <- rep(175000, 21)

#    quota[["200k per year"]] <- rep(200000, 21)
#    quota[["225k per year"]] <- rep(225000, 21)
#    quota[["250k per year"]] <- rep(250000, 21)
#   quota[["275k per year"]] <- rep(275000, 21)
#    quota[["300k per year"]] <- rep(300000, 21)
#    quota[["325k per year"]] <- rep(325000, 21)
#    quota[["350 per year"]] <- rep(350000, 21)
#   quota[["375k per year"]] <- rep(375000, 21)
#    quota[["300k per year"]] <- rep(300000, 21)
#   quota[["325k per year"]] <- rep(325000, 21)

#    quota[["400k per year"]] <- rep(400000, 21)
#    quota[["425k per year"]] <- rep(425000, 21)
#    quota[["450k per year"]] <- rep(450000, 21)
#    quota[["475k per year"]] <- rep(475000, 21)
#    quota[["500k per year"]] <- rep(500000, 21)
#    quota[["525k per year"]] <- rep(525000, 21)
#    quota[["550k per year"]] <- rep(550000, 21)
#    quota[["575k per year"]] <- rep(575000, 21)


#  quota[["Constant quota 400k"]] <- c(rep(500000,5),500000,165000,700000,165000,165000)
#quota[["Constant quota 400k"]] <- c(rep(200000,31)) 
#quota[["Constant quota 400k"]] <- c(rep(60000,81),60000,60000)
# quota[["Changing quotas 400k for 15 y then 100k for next 35y"]] <- c(rep(190000,15),rep(150000,15),rep(20000,40))  
 
  # define one element of the list giving it# a name that will be used in the output file.     this example runs one quota for 15 years, then different quoptas for each of## the remaining 5 years


  ## Run the projection
   resp.table <- lapply(quota, function(x) proj.sealmod.sim.f(result = resultSim, last.data.year = 2019, end.year.sim = 2040, nb.proj = 500, quota = x, fish=F))   # resultSim is the object created by the function SealModel (see Parallel_run script)

 

  ## Print the table indicating the respect of reference values
   print(resp.table)

  ## Save the table indicating the respect of reference values
   if (file.exists(paste(rep, "/modelAMKcor_output/respect_reference.csv", sep=""))){unlink(paste(rep, "/modelAMKcor_output/respect_reference.csv", sep=""))}  # remove the file if it already exists
   sapply(names(quota), function(x){
     write.table(paste("Scenario:", x), file= paste(outdir,"respect_reference.csv", sep=""), row.names = FALSE, col.names = FALSE, append=T)
     suppressWarnings(write.table(cbind(c("Pop", "N70", "N50", "N30"), resp.table[[x]]), file= paste(outdir,"respect_reference.csv", sep=""), append=T, sep=",", row.names=FALSE, col.names=T))
     }
   )

 start.time2 <- Sys.time() 
 
 endtime<-start.time2-start.time1
 
 endtime