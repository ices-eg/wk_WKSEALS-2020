#' Load and prepare data
#'
#' Load and prepare data to run the assessment model using TMB requirements.
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param Amax Maximum age group. Default is 20 years.
#' @param years_of_prediction Number of years in the future to project the model. Default is 15 years.
#' @param Fpast Which fecundity rate to use in past estimations. Fproj = "mean" uses mean value of observed fecundity rates. Otherwise available estimates (linearly interpolated) are used.
#' @param Fproj Which fecundity rate to use in future projections. Fproj = "mean" uses mean value of observed fecundity rates. Otherwise a fixed Fproj can be set.
#' @param catch_quota Catch quota for 0 and 1+animals to be used in future projections. catch_quota = "mean" uses mean value of observed catch for the last 5 years.
#' @param return_fec Return entire fecundity table, including information on years of sampling. Useful for plotting etc.
#' @return data List of loaded data ready for TMB.
#' @keywords input, data
#' @export
#' @examples
#' load.data(population = "harpeast")

load.data <- function(population = "harpeast",Amax = 20,years_of_prediction = 15,
                      Fpast='est', Fproj = 'mean', catch_quota=c(0,0), return_fec=F)
{
  # Read in data ---------------

  # Catch data
  catch_data <- read.table(paste(population,"/catch_data.dat",sep = ""),header = FALSE)
  # Pup production estimates
  pup_production <- read.table(paste(population,"/pup_production.dat",sep = ""),header = FALSE)
  # Available fecundity data
  fecundity <- read.table(paste(population,"/fecundity.dat",sep = ""),header = FALSE)
  # Birth ogives for different periods
  Pdat <- read.table(paste(population,"/wgharp.ogi",sep = ""),sep = "",header = TRUE)
  # Which periods the various birth ogives applies to
  Pper <- read.table(paste(population,"/wgharp.ogp",sep = ""),header = TRUE)
  #Priors used
  priors <- read.table(paste(population,"/priors.dat",sep = ""),header = FALSE)					#Priors used

  years <- c(catch_data[1,1],catch_data[dim(catch_data)[1],1])
  #MAYBE ADD STEPWISE CHANGES IN FECUNDITY AND BIRTH OGIVES INSTEAD OF LINEAR TRANSITION
  if(Fpast=='mean') {
    fecundity[,2] <- rep(mean(fecundity[,2]), nrow(fecundity))
    fecundity[,3] <- rep(mean(fecundity[,3]), nrow(fecundity))
  }
  FecAndP = build.PandF(Fdat = fecundity,Fproj = Fproj,Pdat = Pdat,Pper = Pper,years = years)
  Fdt = FecAndP$Fdt       #SJEKK VERDIEN PÅ DEN SISTE OG SAMMENLIKN MED WGHARP
  Pmat = FecAndP$Pmatrix
  
  if(any(catch_quota=='mean')) {
    catch_quota <- round(apply(tail(catch_data, 5)[,-1], 2, mean))
  }
  
  # Prepare input to the model ------------
  data <- list()
  data$Amax = Amax													#Maximum age group
  data$Cdata = as.matrix(catch_data)								#Add catch data
  data$Nc = length(catch_data[,1])									#Number of years with catch data
  data$pupProductionData = as.matrix(pup_production)					#Pup production estimates
  data$Np = length(pup_production[,1])								#Number of years with pup production data
  data$Ftmp = as.vector(Fdt)                      # Fecundity rates
  data$Pmat = as.matrix(Pmat)									#Birt ogives
  data$Npred = years_of_prediction									#Number of years to run projections
  data$priors = as.matrix(priors)									#Priors for estimated parameters
  data$Npriors = length(priors$V1)									#Number of priors
  data$CQuota = catch_quota											#Catch level in future projections
  data$fecundity = fecundity
  
  if('Pper' %in% names(FecAndP)) {
    data$Pper = FecAndP$Pper
  }
  if(return_fec) {
    data$fecundity <- fecundity
  }
  
  return(data)
}

#' Download Repository Data
#'
#' Downloading the most updated data on harp seals in the West Ice (along the coast of Greenland) and in the East Ice (the White Sea) and the Hooded seal population in the West Ice from the Seal Data repository.
#' @param url The URL of the data repository.
#' @param chooseFolder Logical parameter if you want to specify to which folder to download the data.
#' @return A folder containing the Data for the various populations and a Scripts folder for an example of complete analysis of the data.
#' @keywords input, data
#' @export
#' @examples
#' downloadData(url = "")
downloadData <- function(url = "https://github.com/NorskRegnesentral/HarpAndHoodSealData/archive/master.zip",
                         chooseFolder = TRUE)
{
  if(chooseFolder){
    wdFolder = choose.dir()
    if(!is.na(wdFolder)) setwd(wdFolder)
  }
  
  cat("\n Downloading data from repository...")
  download.file(url = url,
                destfile = "main.zip")
  unzip("main.zip")
  
  #Cleaning up
  cat("\n Cleaning up files....\n")
  #currentFolder = paste(getwd(),"/HarpAndHoodedSealData-master/",sep = "")
  #newFolder = getwd()
  
  #newFolder = ""
  setwd(paste0(getwd(),"/HarpAndHoodSealData-master/"))
  currentFolder = getwd()
  setwd('..')
  newFolder = getwd()
  #newFolder = "../"
  
  file.copy(from = file.path(currentFolder,
                             list.files(currentFolder)),
            to = newFolder,
            overwrite = TRUE,
            recursive = TRUE,
            copy.mode = TRUE)
  
  if (file.exists("main.zip")) 
    #Delete file if it exists
    file.remove("main.zip")
  
  if (file.exists("HarpAndHoodSealData.Rproj")) 
    #Delete file if it exists
    file.remove("HarpAndHoodSealData.Rproj")

    if (file.exists("README.md")) 
    #Delete file if it exists
    file.remove("README.md")

  if (file.exists("HarpAndHoodSealData-master")) 
  #Delete file if it exists
  unlink("HarpAndHoodSealData-master",recursive = TRUE)
  
  
}




#' Construct time varying fecundity rates and birth ogives
#'
#' Prepare the fecundity vector and the birth ogive matrix used by the model. A linear
#' transition between years with missing data is used.
#' @param Fdat Observed fecundity rates.
#' @param Fproj Which fecundity rate used in projections. Default is mean value over all observed fecundity rates ("mean"). If not a specified fecundity rate can be used (just insert number).
#' @param Pdat Observed birth ogives.
#' @param Pper For which time periods the birth ogives are valid.
#' @param years Start and stop year for the model (without projections).
#' @param return.periods Whether to return data frame containing information 
#' about time periods of P and F sampling. NOTE! This is to allow plotting fecundity curves,
#' distinguishing between real sampled values and interpolated between sampling periods. 
#' The returned object, Pper, must be omitted from data that are input into model.
#' 
#' @return Pmatrix Time varying birth ogives.
#' @keywords Fecundity, birth ogive
#' @export
#' @examples
#' FvecandPmat <- build.PandF(fecundity,catch_data,Fproj,Pdat,Pper)

#MAYBE ADD STEPWISE CHANGES IN FECUNDITY AND BIRTH OGIVES INSTEAD OF LINEAR TRANSITION
build.PandF <- function(Fdat = fecundity,Fproj = Fproj,Pdat = Pdat,Pper = Pper,years = years,return.periods=T)
{

  yr1 = years[1]
  yr2 = years[2]
  nyr = yr2 - yr1 + 1

  Fvec = array(0,nyr)
  Fvec[1:(Fdat$V1[1] - yr1 + 1)] = Fdat$V2[1]

  ## M. Biuw 2019/06/11: Added 'if' below statement to avoid error
  ## when using only one fecundity estimate (as for hoods)
  
  if(length(Fdat$V1)>1) {
    for (i in 2:length(Fdat$V1)){
      if(Fdat$V1[i]-Fdat$V1[i-1]==1){
        Fvec[Fdat$V1[i]-yr1+1] = Fdat$V2[i]
      } else {
        i1 = Fdat$V1[i-1] - yr1 + 1
        i2 = Fdat$V1[i] - yr1 + 1
        Fvec[i1:i2] = seq(Fdat$V2[i-1],Fdat$V2[i],length.out = (i2-i1+1))
      }
    }
    Fvec[i2:nyr] = Fdat$V2[i]
  } else {  
    i2 <- Fdat$V1 - yr1 + 1
    Fvec[i2:nyr] = Fdat$V2
  }	  
  
  if(class(Fproj) == "character") Fvec[(yr2-yr1+1):nyr] = mean(Fdat$V2)

  if(class(Fproj)=="numeric") Fvec[(yr2-yr1+1):nyr] = Fproj


  P = matrix(0,nyr,max(Pdat$Age))

  for(i in 1:(Pper$Pstop[1]-yr1+1)){
    P[i,] = Pdat[,2]
  }

  for(i in 2:length(Pper$Pstart)){
    i1 = Pper$Pstop[i-1]-yr1+1
    i2 = Pper$Pstart[i]-yr1+1
    i3 = Pper$Pstop[i]-yr1+1

    for(j in i2:i3){
      P[j,] = Pdat[,i+1]
    }

    for(age in 1:max(Pdat$Age)){
      P[i1:i2,age] = seq(P[i1,age],P[i2,age],length.out = (i2-i1+1))
    }

  }

  for(i in i3:nyr){
    P[i,] = Pdat[,length(Pper$Pstart)+1]
  }

  ## M. Biuw 2019/08/13: Added option for returning data frame with info on 
  ## time periods for P and F sampling.
  
  if(return.periods) {
    return(list(Fdt = Fvec,Pmatrix = P, Pper = Pper))
  } else {
    return(list(Fdt = Fvec,Pmatrix = P))
  }
}

#' Load initial values
#'
#' Load initial values for which to start the optimization from. These are the values to be
#' estimated from the model.
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param fromFile Load initial values from file for specified population or set initial values manually.
#' @param Kinit Initial value for population size.
#' @param Minit Initial value for mortality of 1+ population, i.e., seals of age 1 or more.
#' @param M0init Initial value for pup mortality.
#' @return parameters List of parameters with initial values to be estimated fromt the model.
#' @keywords input, parameters, initial values
#' @export
#' @examples
#' load.initial.values(population = "harpeast")

load.initial.values <- function(population = "harpeast",fromFile = TRUE,Kinit = 2000000,Minit=0.09,M0init=0.27)
{

  if(fromFile == TRUE){
    initial_values <- read.table(paste(population,"/initial_values.dat",sep = ""),header = FALSE)

    #Initial values
    Kinit = initial_values[1,]								#Initial population size
    Minit = initial_values[2,]									#Natural adult mortality
    M0init = initial_values[3,]								#Natural pup mortality
  }


  #Transform some parameters to ensure that estimates are > 0
  parameters <- list()
  parameters$logK= log(Kinit) #THESE SHOULD BE BOUNDED PARAMETERS, FIX THAT LATER
  parameters$Mtilde= logit(Minit)
  parameters$M0tilde= logit(M0init)

  return(parameters)
  }


plot.catch <- function(pop='harpwest', unit=10000) {
  dat <- load.data(pop)$Cdata
  require(RColorBrewer)
  theCols <- brewer.pal(ncol(dat), 'Dark2')
  
  add.alpha <- function(col, alpha=1){
    if(missing(col))
      stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
            rgb(x[1], x[2], x[3], alpha=alpha))  
  }
  
  unit.form <- format(unit, big.mark=' ')
  if(unit!=1) {
    matplot(dat[,1], dat[,-1]/unit, type='n', 
            xlab='', ylab=paste0('Catch level (in ', unit.form, ')'))
  } else {
    matplot(dat[,1], dat[,-1]/unit, type='n', 
            xlab='', ylab='Catch level')
  }
  
  polygon(c(dat[,1], rev(dat[,1])), c(dat[,2]/unit, rep(0, nrow(dat))),
          border=theCols[2], col=add.alpha(theCols[2], 0.7)) 
  
  polygon(c(dat[,1], rev(dat[,1])), c(dat[,3]/unit, rep(0, nrow(dat))),
          border=theCols[1], col=add.alpha(theCols[1], 0.7)) 
  
  legend('topright', col=theCols[c(1:2)], lwd=2, c('1+ catches', 'Pup catches'), bty='n')
}

plot.fecundity <- function(pop='harpwest') {
  dat <- load.data(pop, return_fec=T)
  theXlim=range(dat$Cdata[,1])
  fec <- dat$fecundity
  theYlim=c(0.1*floor(min(fec[,2])*10),
            0.1*ceiling(max(fec[,2])*10))
  require(RColorBrewer)
  fApp <- approx(fec[,1], fec[,2], dat$Cdata[,1], rule=2)
  plot(fec[,1], fec[,2], xlim=theXlim, ylim=theYlim,
       type='n', xlab='', ylab='F')
  lines(fApp$x, fApp$y, lty=2, col=4)
  segments(fec[,1], fec[,2]-fec[,3],
           fec[,1], fec[,2]+fec[,3])
  segments(fec[,1]-0.5, fec[,2]-fec[,3],
           fec[,1]+0.5, fec[,2]-fec[,3])
  segments(fec[,1]-0.5, fec[,2]+fec[,3],
           fec[,1]+0.5, fec[,2]+fec[,3])
  points(fec[,1], fec[,2], pch=21, bg=2)
  legend('bottomleft', pch=c(21, NA), lty=c(NA, 2), pt.bg=c(2, NA), c('Historical fecundities', 'Linear interpolation between periods'), bty='n', col=c(1, 4))
}

plot.ogives <- function(pop='harpwest') {
  dat <- load.data(pop)
  ogi <- dat$Pmat[match(dat$Pper$Pstart, dat$Cdata[,1]),]
  matplot(c(1:dim(dat$Pmat)[2]), t(dat$Pmat), 
          type='l', lty=1, col='lightgrey',
          xlab='Age (years)', ylab='Proportion mature females')
  matlines(c(1:dim(ogi)[2]), t(ogi), lty=1, col=c(2:(dim(ogi)[1]+1)))
  
  leg.txt <- dat$Pper$Pstart
  leg.txt[which(dat$Pper$Pstart!=dat$Pper$Pstop)] <- 
    paste(leg.txt[which(dat$Pper$Pstart!=dat$Pper$Pstop)], 
          dat$Pper$Pstop[which(dat$Pper$Pstart!=dat$Pper$Pstop)], sep='-')
  legend('bottomright', col=c(2:(dim(ogi)[1]+1)), leg.txt, bty='n', lwd=1)
}
