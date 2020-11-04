# OPTIMIZATION FUNCTION #

sealmod.fit.f <- function(params, last.data.year, preg.rate, pup.prod){           #You give initial parameters and the last data year and you get the value of the fonction which is the value of the SS for those 2 paramters

  if (all(params >= c(0.09,0.01,4e+6) & params <= c(1,.06,1.2e+7))){  # check if parameters fall in a correct range

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
    objective3 <- objective1/var(obs.pup.data, na.rm = TRUE) + objective2/var(obs.preg.8plus, na.rm = TRUE)
    #    objective3 <- objective1/var(obs.pup.data, na.rm = TRUE)
  }else{
    objective3 <- 1e+10
  }
  return(objective3)
}

