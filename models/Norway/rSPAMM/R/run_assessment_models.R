#' Load the population model for harp seals and hooded seals
#'
#' Load the model object to be optimized.
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param Amax Maximum age group. Default is 20 years.
#' @param years_of_prediction Number of years in the future to project the model. Default is 15 years.
#' @param Fproj Which fecundity rate to use in future projections. Fproj = "mean" uses mean value of observed fecundity rates. Otherwise a fixed Fproj can be set.
#' @return obj Model object to be optimized.
#' @keywords population model
#' @export
#' @examples
#' load.model.object(population = "harpeast")

load.model.object <- function(dat = data,par = parameters,template='harps_and_hoods_population_model2')
{
  require(TMB)
  cat('Compiling model dll, be patient....\n\n')
  flush.console()
  
  compile(paste0("R/", template, ".cpp"),"-O1 -g",DLLFLAGS="")
 
  cat('\n\nDone!\n')
  flush.console()
  
  cat('\nLoading model dll....\n')
  flush.console()


  dyn.load(dynlib(paste0("R/", template)))
  #dyn.unload(dynlib(paste0("R/", template)))
  #dyn.unload(dynlib("R/harps_and_hoods_population_model2"))
  
  cat('\n\nDone!\n')
  flush.console()
  
  cat('\nBuilding model object....\n')
  flush.console()
  
  ## M. Biuw 2019/08/13: Added safety catch to omit data frame with 
  ## info on time periods for P and F sampling, from data used to build model object.
  
  if('Pper' %in% names(dat)) {
    dat <- dat[-match('Pper', names(dat))]
  }
  
  obj <- MakeADFun(data=dat,parameters=par,DLL=template,silent = TRUE)
  
  cat('Done!\n')
  flush.console()
  
  return(obj)
}

#' Run the loaded population model for harp seals and hooded seals
#'
#' Optimize the model object.
#' @param dat The data used for optimizing the model.
#' @param par The initial values of the parameters to be estimated.
#' @param print2screen Logical parameter to decide wether to output results to screen.
#' @return Output from optimized model.
#' @keywords population model
#' @export
#' @examples
#' run.model()
run.model <- function(data = data,par = parameters,print2screen = TRUE)
{
  
  #Load C part---------------------
  tmbDir <- system.file("libs", package = "rSPAMM")
  if(Sys.info()["sysname"] =="Windows")dyn.load(paste(tmbDir,"/x64/rSPAMM",sep = ""))
  if(Sys.info()["sysname"] =="Linux")dyn.load(paste(tmbDir,"/rSPAMM.so",sep = ""))
  #-----------------------------------
  
  if('Pper' %in% names(data)) {
    data <- data[-match('Pper', names(data))]
  }

    
  obj <- TMB::MakeADFun(data=data,parameters=par,DLL="rSPAMM",silent = TRUE)
  
  opt = nlminb(obj$par,obj$fn,obj$gr)
  
  #Print relevant output to screen
  if(print2screen){
    cat('\n--------------------------------------------------\n')
    if(opt$convergence== 0){
    cat(paste('\n Optimization converged: ',opt$message,'\n'))
    cat('\n\n Parameter estimates\n')
    cat(' -------------------\n')
    cat(paste(' Initial population size: K = ',round(exp(opt$par[1])),'\n'))
    cat(paste(' Pup mortality:          M0 = ',round(ilogit(opt$par[3]),2),'\n'))
    cat(paste(' 1+ mortality:            M = ',round(ilogit(opt$par[2]),2),'\n'))
    } else cat('\nOptimization did not converge \n')
    cat('\n--------------------------------------------------\n')
  }
  retlist = list()
  retlist$obj = obj
  retlist$opt = opt
  return(retlist)
  
}


#' Get results from optimized population model for harp seals and hooded seals
#'
#' Get model results.
#' @param data Data object that was used to fit the model
#' @param optobject The TMB model object and the output from the optimization (a list)
#' @return results Results returned to the workspace or saved to file.
#' @keywords population model
#' @export
#' @examples
#' load.model.object(population = "harpeast")
model.results <- function(data, optobject) 
{
  object = optobject$obj
  opt = optobject$opt
  
  rep=TMB::sdreport(object, getJointPrecision=TRUE)
  rep.matrix = summary(rep)
  rep.rnames = rownames(rep.matrix)
  indN0 = which(rep.rnames=="N0");indN0 <- indN0[-1]
  indN1 = which(rep.rnames=="N1");indN1 <- indN1[-1]
  indNTot = which(rep.rnames=="NTot");indNTot <- indNTot[-1]
  indD1 = which(rep.rnames=="D1");
  indDNmax = which(rep.rnames=="DNmax");   
  indNTotmax = which(rep.rnames == "NTotmax")
##  indN0Current = which(rep.rnames=="N0CurrentYear");
##  indN1Current = which(rep.rnames=="N1CurrentYear");
##  indNTotCurrent = which(rep.rnames=="NTotCurrentYear");
## Above lines replaced, seems TMB has wrong indexes for 
## N0Current and NTotCurrent
  
  yrs = c(min(data$Cdata[,1]):(max(data$Cdata[,1])+data$Npred+1))
  
  curYr <- match(max(data$Cdata[,1]), yrs)
  indN0Current <- indN0[curYr]
  indN1Current <- indN1[curYr]
  indNTotCurrent <- indNTot[curYr]
  
  #Extract parameters
  Kest = exp((rep.matrix[1,1]))
  Kll = exp(rep.matrix[1,1]-1.96*rep.matrix[1,2])
  Kest.sd = (Kest-Kll)/1.96
  
  Mest = ilogit((rep.matrix[2,1]))
  Mll = ilogit(rep.matrix[2,1]-1.96*rep.matrix[2,2])
  Mest.sd = (Mest-Mll)/1.96
  
  
  M0est = ilogit((rep.matrix[3,1]))
  M0ll = ilogit(rep.matrix[3,1]-1.96*rep.matrix[3,2])
  M0est.sd = (M0est-Mll)/1.96
  
  
  D1 = rep.matrix[indD1,1]
  DNmax = rep.matrix[indDNmax,1]
  N0Current = rep.matrix[indN0Current,1]
  N1Current = rep.matrix[indN1Current,1]
  NTotCurrent = rep.matrix[indNTotCurrent,1]
  NTotmax = rep.matrix[indNTotmax,1]
  D1.sd = rep.matrix[indD1,2]
  DNmax.sd = rep.matrix[indDNmax,2]
  N0Current.sd = rep.matrix[indN0Current,2]
  N1Current.sd = rep.matrix[indN1Current,2]
  NTotCurrent.sd = rep.matrix[indNTotCurrent,2]
  
  ## Some suggestions:
  cur.yr <- dim(data$Cdata)[1]
  N1Current <- rep.matrix[indN1[cur.yr],1]
  N1Current.sd <- rep.matrix[indN1[cur.yr],2]
    
  res = list(rep=rep, rep.matrix=rep.matrix, rep.rnames=rep.rnames, indN0=indN0,
       indN1=indN1, indNTot=indNTot, indD1=indD1, indNmax=indDNmax, indN0Current=indN0Current,
       indN1Current=indN1Current, indNTotCurrent=indNTotCurrent, 
       years=yrs, Kest=Kest,Kest.sd = Kest.sd, Mest=Mest,Mest.sd=Mest.sd, M0est=M0est, M0est.sd = M0est.sd, 
       D1=D1, DNmax=DNmax, N0Current=N0Current, N1Current=N1Current, 
       NTotCurrent=NTotCurrent, NTotmax = NTotmax, D1.sd=D1.sd, DNmax.sd=DNmax.sd,
       N0Current.sd=N0Current.sd, N1Current.sd=N1Current.sd, NTotCurrent.sd=NTotCurrent.sd)
  return(res)
  }


