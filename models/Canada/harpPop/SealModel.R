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


# print(all<-system.time(resultSim <- SealModel(nbrep = 20, params= c(.2,.04, 9000000), last.data.year=2013, preg.cor=0.85, NoConvLimit = 2, clust=3)))


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

SealModel <- function(nbrep,
                      params,
                      last.data.year,
                      preg.cor = 0.90,
                      NoConvLimit = 2,
                      clust= min(c(nbrep, detectCores(logical=TRUE)-1)),
                      ShowAgeStruct = c(TRUE, TRUE)){


  ## RESAMPLING AND OPTIMIZATION PROCESS ######################################################################

  optimization <- function(params, last.data.year, sealmod.params, preg.cor){

    ## Resampling pregnancy rates into binomials with defined correlation
    #  install.packages("copula")
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
                      "gamma.pup",
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
  write.csv(ParamStats, paste("output/ParamStats_for_", nbrep, "_iterations_InitParams_", paste(params, collapse="_"), "_PregCor_", preg.cor, ".csv", sep=""))


  ## Output Age structure

  if (ShowAgeStruct[1]) {
    # Age structure summary (seals numbers by age class)
    AgeStructSummary <- list(prob.025 = apply(AgeStruct, c(1,2), function(x) quantile(x, prob = 0.025)), prob.25 = apply(AgeStruct, c(1,2), function(x) quantile(x, prob = 0.25)), med = apply(AgeStruct, c(1,2), median),
                             prob.50 = apply(AgeStruct, c(1,2), function(x) quantile(x, prob = 0.50)), prob.75 = apply(AgeStruct, c(1,2), function(x) quantile(x, prob = 0.75)), prob.975 = apply(AgeStruct, c(1,2), function(x) quantile(x, prob = 0.975)))
    for(i in 1:length(AgeStructSummary)){ AgeStructSummary[[i]] <- data.frame(c(seq(0,14), "15+"), AgeStructSummary[[i]]) ; colnames(AgeStructSummary[[i]]) <- c("", seq(1952, last.data.year+1))}


    # Save raw numbers

    write.table(rbind("### Median"), paste("output/AgeStruct_for_", nbrep, "_iterations.csv"), row.names=F, col.names=F)
    suppressWarnings(write.table(AgeStructSummary$med, paste("output/AgeStruct_for_", nbrep, "_iterations.csv"), sep=",", append=TRUE, row.names=F))

    write.table(rbind("","### Prob .025"), paste("output/AgeStruct_for_", nbrep, "_iterations.csv"), append=T, row.names=F, col.names=F)
    suppressWarnings(write.table(AgeStructSummary$prob.025, file=paste("output/AgeStruct_for_", nbrep, "_iterations.csv"), sep=",", row.names=F, append = TRUE))

    write.table(rbind("", "### Prob .25"), paste("output/AgeStruct_for_", nbrep, "_iterations.csv"), append=T, row.names=F, col.names=F)
    suppressWarnings(write.table(AgeStructSummary$prob.25, file=paste("output/AgeStruct_for_", nbrep, "_iterations.csv"), sep=",", row.names=F, append = TRUE))

    write.table(rbind("", "### Prob .75"), paste("output/AgeStruct_for_", nbrep, "_iterations.csv"), append=T, row.names=F, col.names=F)
    suppressWarnings(write.table(AgeStructSummary$prob.75, file=paste("output/AgeStruct_for_", nbrep, "_iterations.csv"), sep=",", row.names=F, append = TRUE))

    write.table(rbind("", "### Prob .975"), paste("output/AgeStruct_for_", nbrep, "_iterations.csv"), append=T, row.names=F, col.names=F)
    suppressWarnings(write.table(AgeStructSummary$prob.975, file=paste("output/AgeStruct_for_", nbrep, "_iterations.csv"), sep=",", row.names=F, append = TRUE))
  }


  #########################################################################################################
  ################################## Plots ################################################################
  #########################################################################################################

  quartz(18,12)
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
  dev.copy2pdf (file=paste("output/Fitting_output_for_",nbrep,"_iterations_InitParams_", paste(params, collapse="_"), "_PregCor_", preg.cor, ".pdf", sep=""))
  #return(list(MakeCluster=t1,Calculation=t2))

  if (ShowAgeStruct[2]){

    # Age structure summary 2 (proportion of the total number in each age class)
    AgeStructPropSummary <- list(prob.025 = apply(AgeStructProp, c(1,2), function(x) quantile(x, prob = 0.025)), prob.25 = apply(AgeStructProp, c(1,2), function(x) quantile(x, prob = 0.25)), med = apply(AgeStructProp, c(1,2), median),
                                 prob.50 = apply(AgeStructProp, c(1,2), function(x) quantile(x, prob = 0.50)), prob.75 = apply(AgeStructProp, c(1,2), function(x) quantile(x, prob = 0.75)), prob.975 = apply(AgeStructProp, c(1,2), function(x) quantile(x, prob = 0.975)))
    for(i in 1:length(AgeStructPropSummary)){ AgeStructPropSummary[[i]] <- data.frame(c(seq(0,14), "15+"), AgeStructPropSummary[[i]]) ; colnames(AgeStructPropSummary[[i]]) <- c("", seq(1952, last.data.year+1))}

    listPage <- vector()
    for ( i in 1:length(1952:(last.data.year + 1))){
      if(i == 1 | i%%20 == 0){
        quartz(18,14)
        par(mfrow=c(4,5))
        listPage <- c(listPage, dev.cur())
      }
      barplot(AgeStructPropSummary$med[1:15,i+1], cex.names = 0.9, border="darkgrey", space=0, names.arg = AgeStructPropSummary$med[1:15,1], horiz = TRUE,las=1, main=colnames(AgeStructPropSummary$med)[i+1], xlim=c(0,0.25))
      mtext(paste("15+ =", round(AgeStructPropSummary$med[16,i+1],2)), side = 3, adj=-0.15, cex=0.7, col="darkgrey")
      mtext(paste("Age 1-14\n=", round(sum(AgeStructPropSummary$med[2:15,i+1]),2), "  "), side = 3, adj=0.8, line=-4, cex=0.7, col="darkgrey")
    }
    for (i in 1:length(listPage)){
      dev.set(which = listPage[i])
      dev.copy2pdf (file=paste("output/AgeStructure_output_for_",nbrep,"_iterations_page", i,"-",length(listPage), ".pdf", sep=""))
    }
  }


  ## Final output for SealModel
  return(list("Params"=result,"Pupprod"=pup,"Totpop"=tot,"LastPopVect"=last.pop.vec, "Preg8plus"=preg8plus, "PregSample"=pregsampl, "PupSample"=pupSampl, "Preg8Sample" = preg8Sampl, "RunTime" = list("Make Cluster"=t1, "Parallel Calculations"=t2)))

}
