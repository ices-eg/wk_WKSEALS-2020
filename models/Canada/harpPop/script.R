##########################################################
### This is a reworking of the Hammill Harp Seal model ###
### For the original see PG SIM MODEL Jan 4 2018       ###
##########################################################

# This script is made up of many different functions...
# The function have been seperated out into seperate R scripts that are then sourced when needed...
# The data prep code at the start of this script could also be pulled out...
# Then wrap the intro script into a markdown document...
# Run as a project so all data can be imported using dat/...
# Create folder to output data too and call it output/...
# There is a bug with the tkprogress button so replace with txtProgressBar

rm(list = ls())         #Delete all the current files in the workspace

##############################################################################################
###---------------------------------------------------------------------------------------####
###                      LOAD AND ORGANIZE THE MODEL PARAMETERS                           ####
###---------------------------------------------------------------------------------------####
##############################################################################################

start.time1 <- Sys.time()

### 1. SMOOTHED REPRODUCTIVE RATES

# Load data
data <- read.csv("dat/temppr_2017.csv", sep = ";", header = F)
colnames(data) <- c("years", "age", "total", "npreg")
data$total[data$total == 0] <- NA
nyears <- length(unique(data$years))

# create 2 empty matrices
predictyears <- seq(1952, 2018, 1) # set years where smoothed data are needed
span <- range(predictyears)[2] - range(predictyears)[1] + 1
preg.allyears.log <- matrix(NA, length(unique(data$age)), length(predictyears))
colnames(preg.allyears.log) <- as.character(predictyears)
se.allyears.log <- matrix(NA, length(unique(data$age)), length(predictyears))
colnames(se.allyears.log) <- as.character(predictyears)

# calculate smoothed data
#install.packages("locfit")
library(locfit)
degoffree <- 2 # set degree of freedom
nearestn <- c(0.75, 0.8, 0.95, 0.95, 0.95)   # set nearest neighbor values chosen with LCV and AIC diagnostics measures (see Diagnostics part)
# for age 4,5,6,7, and 8 respectively

for (i in 1:length(unique(data$age))){
  subbyage <- subset(data, data$age == unique(data$age)[i])
  fit <- locfit(npreg ~ lp(years, deg = degoffree, nn = nearestn[i]), weights = total, data = subbyage, family = "binomial")
  preg.allyears.log[i,] <- preplot(fit, data.frame(years = predictyears))$fit
  se.allyears.log[i,] <- preplot(fit, data.frame(years = predictyears), band = "local")$se
}

# for information give matrix of smoothed value from preg rates for each predicted year
nb_ageplus <- 18  # set number of classes age+

preg.young <- matrix(0, 3, length(predictyears))
preg.all <- expit(preg.allyears.log)
preg.plus <- matrix(preg.all[length(unique(data$age)),], nb_ageplus, length(predictyears), byrow = T)

preg.smooth <- rbind(preg.young, preg.all, preg.plus)

preg.young.CIs <- matrix(0, 3, length(predictyears))
preg.all.CIupper <- expit(preg.allyears.log + se.allyears.log * 1.96)
preg.all.CIlower <- expit(preg.allyears.log - se.allyears.log * 1.96)
preg.plus.CIupper <- matrix(preg.all.CIupper[length(unique(data$age)),],nb_ageplus, length(predictyears), byrow = T)
preg.plus.CIlower <- matrix(preg.all.CIlower[length(unique(data$age)),],nb_ageplus, length(predictyears), byrow = T)
preg.smooth.CIupper <- rbind(preg.young.CIs, preg.all.CIupper, preg.plus.CIupper)
preg.smooth.CIlower <- rbind(preg.young.CIs, preg.all.CIlower, preg.plus.CIlower)

# Plots smoothed preg.rates and 95CIs for age classes 4 to 8+
pdf(file = paste("output/", "Smoothed preg rates and 95CIs.pdf", sep = ""))
par(mfrow = c(2, 3))
for(i in 4:8){
  plot(preg.smooth[i,] ~ colnames(preg.smooth), ylim = c(0,1), xlab = "years", ylab = "preg.rate", main = paste("age", i), type = "l", col = 1)
  lines(preg.smooth.CIupper[i,] ~ colnames(preg.smooth), col = "blue", lty = 2)
  lines(preg.smooth.CIlower[i,] ~ colnames(preg.smooth), col = "blue", lty = 2)
  points((data$npreg/data$total)[data$age == i] ~ data$years[data$age == i], pch = 16, col = "red", cex = 0.7)
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

pb1 <- txtProgressBar(min = 0, max = nresample, style = 3) #Progress bar

rubix <- array(dim = c(26, span, nresample))        #  create an empty array to store the resampling results

for(r in 1:nresample) {                       # this loops resamples 10,000 times within the smoothers' distribution to generate an SD

  setTxtProgressBar(pb1, r) # update progress bar

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

# from this SD, we can estimate the size n of a binomial distribution binom(n,p) that would have the same SD
sd.mat <- apply(rubix, c(1, 2), sd)
p.mat <- preg.smooth
ncalc.mat <- round(p.mat * (1 - p.mat) / (sd.mat^2))
ncalc.mat[is.nan(ncalc.mat)] <- 1
ncalc.mat[ncalc.mat < 1] <- 1

# Include the estimated size n in the output describing smoother predictions
n_estim<-vector()

for (i in (4:8)){
  n_estim <- c(n_estim, ncalc.mat[i,])
}

output <- data.frame(output, n_estim = as.numeric(t(n_estim)))
write.csv(output, "output/output_smoother_model.csv", row.names = FALSE)


### 3.  SMOOTHED REPRODUCTIVE RATES (to use in the binomial distribution)

npreg <- data[,4]
npreg.mat <- matrix(npreg, 5, nyears, byrow = T)
colnames(npreg.mat) <- data[1:nyears, 1]

nsize <- data[,3]
nsize.mat <- matrix(nsize, 5, nyears, byrow = T)
colnames(nsize.mat) <- data[1:nyears, 1]

preg.young <- matrix(0, 3, nyears)
preg.plus <- matrix(npreg.mat[5,], 18, nyears, byrow = T)
preg.all <- rbind(preg.young, npreg.mat, preg.plus)

size.young <- matrix(0, 3, nyears)
size.plus <- matrix(nsize.mat[5,], 18, nyears, byrow = T)
size.all <- rbind(size.young, nsize.mat, size.plus)

# create 2 empty matrices with all years and age classes
preg.allyears <- matrix(NA, 26, span)
colnames(preg.allyears) <- as.character(predictyears)
size.allyears <- matrix(NA, 26, span)
colnames(size.allyears) <- as.character(predictyears)

# find which years have data
match.years <- match(colnames(preg.all), colnames(preg.allyears))

# replace values in allyears matrix by preg counts for years with data
preg.allyears[,match.years] <- preg.all
size.allyears[,match.years] <- size.all

# transform pregnant counts into proportions
preg.allyears <- preg.allyears/size.allyears

# turn proportions of 100% into 95% (otherwise, there is no variability) # note that Binom(22,0.95): mean = 0.95, 2.5%-97.5% quantiles = 0.82-1.00

preg.allyears[preg.allyears == 1] <- 0.95

# remove NaN's
preg.allyears[is.nan(preg.allyears)] <- NA

# find cells with NA or sample size < 40
toreplace <- which(is.na(size.allyears) | size.allyears < 40)

# find cells with NA or sample size < 40 for 8+ preg rate
smoothedvalues <- which(is.na(size.allyears[8,]) | size.allyears[8,] < 40)

# replace these cells with smoothed preg rates and sample size = n as calculated above from SD of smoothed data
preg.allyears[toreplace] <- preg.smooth[toreplace]
size.allyears[toreplace] <- ncalc.mat[toreplace]

### Include all parameters in one object
sealmod.params<-list()     #The parameters of the model are defined as a list (suite of objects)

# Final age specific reproductive rates
sealmod.params$preg.rate <- preg.allyears

# Final age specific sample sizes for reproductive rates
sealmod.params$preg.size <- size.allyears

# Proportion of pups surviving an unusual mortality event (1952-2005)
sealmod.params$prop.surv <- as.matrix(read.table("dat/prop-surv-1952_quantified.csv", sep = ";", header = TRUE))

# Pup production mean estimates based on mark-recapture experiments and aerial surveys
sealmod.params$pup.prod <- as.matrix(read.table("dat/pup-prod-est-1952.csv", sep = ";", header = TRUE, na.strings = "NA"))
#sealmod.params$pup.prod[57] <-1630000

# Pup production SE estimates based on mark-recapture experiments and aerial surveys
sealmod.params$pup.prod.se <- as.matrix(read.table("dat/pup-prod-est-se-1952.csv", sep = ";", header = TRUE, na.strings = "NA"))
# sealmod.params$pup.prod.se[57] <-110381

# Initial vector of abundance (year 1952)
sealmod.params$init.pop <- as.matrix(read.table("dat/initial-pop-1952.csv", sep = ";", header = FALSE))

# Observed reproductive rate for class 8+
sealmod.params$obs.preg.8plus <- rep(NA, length(1952:2014))
sealmod.params$obs.preg.8plus[match.years] <- npreg.mat[5,] / nsize.mat[5,]

# Observed sample size for class 8+
sealmod.params$obs.pregsize.8plus <- rep(NA, length(1952:2014))
sealmod.params$obs.pregsize.8plus[match.years] <- nsize.mat[5,]

# Natural mortality multiplier for first year seals
gamma.pup <- 3

raw.remov <- as.matrix(read.table("dat/raw-removal-1952.csv",sep = ";", header = TRUE))
prop.arctic.pup <- 0.034
prop.green.pup <- scan("dat/prop-green-pup-1952.csv",sep = ";",quiet=T)


##############################################################################################
###---------------------------------------------------------------------------------------####
###                           FIT TO PAST AND PRESENT DATA                                ####
###---------------------------------------------------------------------------------------####
##############################################################################################

start.time2 <- Sys.time()

#this is where the seal.mod.sim.f function was...
source("R/sealmodsim.r")

### OPTIMIZATION FUNCTION ###

# this is where the seal.mod.fit.f function was
source("R/sealmodfit.r")

#################### SAMPLING AND PARALLEL RUN ############################################################################

# this is where the SealModel function was
source("R/sealmodel.r")

#############################################################################################
### Create cluster of computer processor threads for the calculation (Parallel computing) ###
#############################################################################################

# Choose the number of thread to be included in the cluster (keep one thread off to let windows run its jobs on this one)
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
    while(nchar(nbthreads) == 0) nbthreads <- select.list(as.character(seq(1, detectCores(logical = TRUE)-1)), multiple = FALSE, title = "Select the number of threads")
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

print(
  system.time(
    resultSim <- SealModel(nbrep = 4000,
                           params= c(.1,.02, 9000000),
                           last.data.year = 2017,
                           preg.cor = 0.85,
                           NoConvLimit = 2,
                           clust = cltest,
                           ShowAgeStruct = c(TRUE, FALSE))
  )
) # Use cluster created on the previous lines. Just rerun this line to make another run of the model

#Takes about 30 minutes to run 4000 iterations

## Show internal run time
cat("\n### Internal SealModel Run Time ###\n")
print(resultSim$RunTime)
alarm()
dev.off()


################################################################################
#### Check Residual distribution around fitted values (Asked by a reviewer) ####
################################################################################
# Histogram of predicted vs observed values (note: observed values are resampled at each repetition)
dev.new(); hist(as.vector(resultSim$Pupprod) - as.vector(resultSim$PupSample), 50, main="Residuals distribution around estimates of pup number")
PupSample <- resultSim$PupSample
PupSample[,c(27:29,32)]<-NA
# Histogram of predicted vs observed values (note: observed values are resampled at each repetition)
dev.new(); hist(as.vector(resultSim$Pupprod) - as.vector(PupSample), 50, main="Residuals distribution around estimates of pup number\nwithout pup estimates of 1978-1980 and 1983")
dev.new(); hist(as.vector(resultSim$Preg8plus) - as.vector(resultSim$Preg8Sample), 50, main="Residuals distribution around preg8+ estimates")


#################################################################################################
###                         RUN THE PROJECTION FUNCTION                                       ###
#################################################################################################

source("R/projsealmodsim.r") # This is where the projection function was written

## Indicate scenario of removal -you can have the same harvest every year, or uncomment and have different harvest for each year of future simulation
quota<-list()
#   quota[["500k per year"]] <- rep(500000, 50) # one for each projection year
#   quota[["400k per year"]] <- rep(400000, 5)
#   quota[["500k per year"]] <- rep(500000, 5)
#  quota[["Constant quota 400k"]] <- c(rep(500000,5),500000,165000,700000,165000,165000)

quota[["Constant quota 400k"]] <- c(rep(400000,4),400000,400000)

# quota[["Changing quotas 400k for 15 y then 100k for next 35y"]] <- c(rep(190000,15),rep(150000,15),rep(20000,40))

# define one element of the list giving it a name that will be used in the output file
# this example runs one quota for 15 years, then different quoptas for each of## the remaining 5 years

## Run the projection
resp.table <- lapply(quota, function(x)
  proj.sealmod.sim.f(result = resultSim,
                     last.data.year = 2017,
                     end.year.sim = 2023,
                     nb.proj = 4000,
                     quota = x,
                     fish = F)
)   # resultSim is the object created by the function SealModel (see Parallel_run script)

## Print the table indicating the respect of reference values
print(resp.table)

## Save the table indicating the respect of reference values
if (file.exists("output/respect_reference.csv")){unlink("output/respect_reference.csv")}  # remove the file if it already exists
sapply(names(quota), function(x){
  write.table(paste("Scenario:", x), file= "output/respect_reference.csv", row.names = FALSE, col.names = FALSE, append=T)
  suppressWarnings(write.table(cbind(c("Pop", "N70", "N50", "N30"), resp.table[[x]]), file = "output/respect_reference.csv", append = T, sep = ",", row.names = FALSE, col.names = T))
})

start.time2 <- Sys.time()
endtime <- start.time2 - start.time1
endtime

#ends
