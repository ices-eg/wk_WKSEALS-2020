#' Internal helper function for finding equilibrium quota
#'
#' Finding D-based equilibrium quota
#' @param Tot Total annual catch for population projections 
#' @param dataD The data to be analyzed
#' @param parametersD Parameters used for the model
#' @param quota Proportional catch of 0 and 1+ animals
#' @return minD Minimum found for difference between the depletion coefficient D and 1
#' @keywords population model
#' @export

eq.quota.helper.D <- function(Tot,dataD,parametersD,quota = c(0,1))
{
  # Function called from find.eq.quota, that should be used!
  #dataD <- load.data(population=population, catch_quota=Tot*quota)
  #parametersD <- load.initial.values(population)
  dataD$CQuota = Tot*quota
  #objD <- load.model.object(dataD, parametersD)
  objD = TMB::MakeADFun(data=dataD,parameters=parametersD,DLL="rSPAMM",silent = TRUE)
  
  optD = nlminb(objD$par,objD$fn,objD$gr)
  repD = TMB::sdreport(objD, getJointPrecision=TRUE)
  abs(1-repD$value[match('D1', names(repD$value))])
}


#' Internal helper function for finding equilibrium quota
#'
#' Finding N-based equilibrium quota
#' @param Tot Total annual catch for population projections 
#' @param dataD The data to be analyzed
#' @param parametersD Parameters used for the model
#' @param quota Proportional catch of 0 and 1+ animals
#' @return minN Minimum found for difference between current and projected 1+ population size
#' @keywords population model
#' @export

eq.quota.helper.N1 <- function(Tot,dataN1,parametersN1,quota)
{
  # Function called from find.eq.quota, that should be used!
  #dataN1 <- load.data(population=population, catch_quota=Tot*quota)
  #parametersN1 <- load.initial.values(population)
  dataN1$CQuota = Tot*quota
  
  #objN1 <- load.model.object(dataN1, parametersN1)  
  objN1 = TMB::MakeADFun(data=dataN1,parameters=parametersN1,DLL="rSPAMM",silent = TRUE)

  
  optN1 = nlminb(objN1$par,objN1$fn,objN1$gr)
  repN1 = TMB::sdreport(objN1, getJointPrecision=TRUE)
  N1Cur = repN1$value[match("N1CurrentYear", names(repN1$value))]
  N1Proj = repN1$value[max(which(names(repN1$value) == "N1"))]
  #N1Proj = repN1$value[match("DNmax", names(repN1$value))]
  abs(N1Cur-N1Proj)
}


#' Main function for finding equilibrium quota
#'
#' Finding D-based or N-based equilibrium quota
#' @param MIN Minimum extent of search window for total catch 
#' @param MAX Maximum extent of search window for total catch 
#' @param data The data to be analyzed
#' @param parameters Parameters used for the model
#' @param quota Proportional catch of 0 and 1+ animals
#' @param method Set whether D-based (depletion coefficient D) or N-based (1+ population size) criterion should be used for optimisation  (Dbased/Nbased)
#' @return qEq Optimum quota for achieving equilibrium projected population size
#' @keywords population model
#' @export

find.eq.quota <- function(MIN=1000,
                          MAX=40000,
                          quota=c(0,1),
                          data = data,
                          parameters = parameters,
                          method = "Dbased")
{
  # Function to find equilibrium quota
  quota = quota/sum(quota)
  if (method == "Dbased"){
    tmp = optimize(eq.quota.helper.D,lower=MIN,upper=MAX,dataD = data, parametersD = parameters,quota=quota,tol=50)
  } else
  {
    tmp = optimize(eq.quota.helper.N1,lower=MIN,upper=MAX,dataN1 = data, parametersN1 = parameters,quota=quota,tol=50)
  }
  cat("-----------------------------------------------------\n\n")
  cat("Estimated equilibrium:\n")
  cat("Pups:  ",round(tmp$minimum*quota)[1],"\n")
  cat("Adults:",round(tmp$minimum*quota)[2],"\n")
  cat("Total :",sum(round(tmp$minimum*quota)),"\n\n")
  cat("-----------------------------------------------------")
  invisible(tmp$minimum*quota)
}


#' Internal helper function for finding N70 quota
#'
#' Finding N-based N70 quota
#' @param Tot Total annual catch for population projections 
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param quota Proportional catch of 0 and 1+ animals
#' @param predYears Define how many years the reduction will be carried out on. Default (10)
#' @return minN70 Minimum found for difference between N70 and the projected total population size
#' @keywords population model
#' @export

N70.helper.Nmax <- function(Tot,dataNmax,parametersNmax,quota,predYears)
{

  dataNmax$CQuota = Tot*quota
  objNmax = TMB::MakeADFun(data=dataNmax,parameters=parametersNmax,DLL="rSPAMM",silent = TRUE)
  optNmax = nlminb(objNmax$par,objNmax$fn,objNmax$gr)
  repNmax = TMB::sdreport(objNmax, getJointPrecision=TRUE)

  indNTot = which(names(repNmax$value)=="NTot")
  indNTot = indNTot[-1]
  indCur = diff(range(dataNmax$Cdata[,1]))+1
  NTot = repNmax$value[indNTot]
  NTotSD = repNmax$sd[indNTot]
  N70 = 0.7*max(NTot[c(1:indCur)])
  
  #Lower limit of 80 percent confidence interval
  Npred = NTot[indCur+predYears]-qnorm(1-0.1)*NTotSD[indCur+predYears]
  
  #if(Npred>0) {
  return(abs(N70-Npred))
  #} else {
  #  99999
  #}  
}


#' Internal helper function for finding N70 quota
#'
#' Finding D-based N70 quota
#' @param Tot Total annual catch for population projections 
#' @param dataD The data to be analyzed
#' @param parametersD Parameters used for the model
#' @param quota Proportional catch of 0 and 1+ animals
#' @param predYears Define how many years the reduction will be carried out on. Default (10)
#' @return minD07 Minimum found for difference between the depletion coefficient D and 0.7
#' @keywords population model
#' @export

N70.helper.D <- function(Tot,dataD,parametersD,quota,predYears)
{
  # Function called from find.N70 that should be used!
  #dataD <- load.data(population=population, catch_quota=Tot*quota)
  #parametersD <- load.initial.values(population)
  dataD$Npred = predYears
  dataD$CQuota = Tot*quota
  
  #objD <- load.model.object(dataD, parametersD)
  objD <- TMB::MakeADFun(data=dataD,parameters=parametersD,DLL="rSPAMM",silent = TRUE)
  
  optD <- nlminb(objD$par,objD$fn,objD$gr)
  repD <- TMB::sdreport(objD, getJointPrecision=TRUE)

  DNmax = repD$value[match('DNmax', names(repD$value))]
  DNmaxSD = repD$sd[match('DNmax', names(repD$value))]

  Dpred = DNmax-qnorm(1-0.1)*DNmaxSD
  
  #if(Dpred>0) {
  return(abs(0.7-Dpred))
  #} else {
  #  99999
  #}  
}


#' Main function for finding N70 quota
#'
#' Finding the catch level that brings the population down to N70 with probability 0.8 in a 10 year period
#' @param MIN Minimum extent of search window for total catch 
#' @param MAX Maximum extent of search window for total catch 
#' @param quota Proportional catch of 0 and 1+ animals
#' @param data The data to be analyzed
#' @param parameters Parameters used for the model
#' @param predYears Define how many years the reduction will be carried out on. Default (15)
#' @param method Set whether D-based (depletion coefficient D) or N-based (total population size) criterion should be used for optimisation (Dbased,Nbased)
#' @return q70 Optimum quota for achieving projected size of 70% of maximum population size
#' @keywords population model
#' @export

find.N70.quota <- function(MIN=100,
                           MAX=50000,
                           quota=c(0,1),
                           data = data,
                           parameters = parameters,
                           predYears = 15,
                           method = "Nbased")
{
  #Check if current population is below N70
  #If current population is below N70 return
  objTest = TMB::MakeADFun(data=data,parameters=parameters,DLL="rSPAMM",silent = TRUE)
  optTest = nlminb(objTest$par,objTest$fn,objTest$gr)
  repTest = TMB::sdreport(objTest, getJointPrecision=TRUE)
  
  indNTot = which(names(repTest$value)=="NTot")
  indNTot = indNTot[-1]
  indCur = diff(range(data$Cdata[,1]))+1
  NTot = repTest$value[indNTot]
  NTotSD = repTest$sd[indNTot]
  NTotCur = NTot[indCur]
  N70 = 0.7*max(NTot[c(1:indCur)])
  Npred = NTot[indCur+predYears]-qnorm(1-0.1)*NTotSD[indCur+predYears]
  
  if(Npred>N70) isAbove = TRUE else isAbove = FALSE
  
  
  # Function to find 70% quota
  if(isAbove){
    quota = quota/sum(quota)
    if (method == "Dbased"){
      tmp = optimize(N70.helper.D,lower=MIN,upper=MAX,dataD = data,parametersD = parameters,quota=quota,predYears = predYears,tol=5)
      }
    if (method == "Nbased") {
      tmp = optimize(N70.helper.Nmax,lower=MIN,upper=MAX,dataNmax=data,parametersNmax=parameters,quota=quota,predYears = predYears,tol=5)
    }
    #cat("N70 quota for",population,"(pups,adults,total):",round(tmp$minimum*quota),sum(round(tmp$minimum*quota)),"\n")
    cat("-----------------------------------------------------\n\n")
    cat("Estimated N70 quota:\n")
    cat("Pups:  ",round(tmp$minimum*quota)[1],"\n")
    cat("Adults:",round(tmp$minimum*quota)[2],"\n")
    cat("Total :",sum(round(tmp$minimum*quota)),"\n\n")
    cat("-----------------------------------------------------")
    tmp$minimum*quota
  } else{
    cat("\n ---------------------------------------\n\n")
    cat(" Current population size is already within\n")
    cat(" the 80 percent confidence interval of\n")
    cat(" the ",predYears," year prediction.\n")
    cat(" Hence, no catch level will be estimated.\n")
    cat("\n ---------------------------------------\n\n")
    return(NA)
  }
}



#' Function for finding PBR quota
#'
#' Finding Potential Biological Removal (PBR) quotas
#' @param n0 Estimated current population size of 0 animals
#' @param n1 Estimated current population size of 1+ animals
#' @param se0 Standard error of estimate for n0
#' @param se1 Standard error of estimate for n1
#' @param rMax Maximum rate of population increase (by default 0.12, commonly used for pinnipeds)
#' @param Fr Assumed recovery factor
#' @param quota Proportional catch of 0 and 1+ animals
#' @return Returns Minimum projected population size (Nmin) and Potential Biological Removal (PBR),
#' along with the PBR divided by 0 and 1+ animals given specified quota 
#' @keywords population model
#' @export

PBR <- function(n0=n0, 
                n1=n1, 
                se0=se0, 
                se1=se1,
                rMax=0.12, 
                Fr=0.5, 
                quota=c(0.14,1-0.14), 
                cv=NA) {
  
  if(is.na(cv)) cv <- sqrt((se0^2) + (se1^2)+(2*se0*se1))/(n0+n1)
  
  Nmin <- round((n0+n1) / exp(0.842*sqrt(log(1+cv^2))))
  pbr <- round(0.5 * rMax * Fr * Nmin)
  
  quota <- as.vector(quota)
  
  list(Nmin=Nmin, CV=cv, PBR=pbr, n0catch=round(pbr*quota[1]), n1catch=round(pbr*quota[2]))
  
}
