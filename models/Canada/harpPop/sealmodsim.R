# The result of the function is the pup production, the total population and the last year population vector

sealmod.sim.f <- function(params,
                          last.data.year,
                          preg.rate,
                          pup.prod){

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
    preg.8plus.sim[year+1] <- ifelse(0.88 * (1 - (sum(seal.sim[,year])/K)^2.4) < 0, 0, (0.88 * (1 - (sum(seal.sim[,year])/K)^2.4))) # 0.88 is considered as max preg rate (can be defined as max preg rate when pop is low)
    #   preg.rate[8:26, year+1] <- rep (preg.8plus.sim[year+1], 19)   # Uncomment this line to allow the model to use preg.rate for 8+ modified by the density dependance formula in the pup prod fitting
    #   preg.8plus.sim[year+1] <- preg.rate[8,year+1]  # When not using the two previous lines, uncomment this one


    # Seal numbers at age 1
    seal.age1 <- (seal[1] * sealmod.params$prop.surv[year] - remov[1]) * exp(-gamma.pup * M) *(1-((sum(seal.sim[,(year)])/K)^2.4))      #seal[1,] = pups

    # Seal numbers for age > 1 and < 25
    seal.num <- (seal[2:24]*exp(-M/2) - remov[2:24]) * exp(-M/2)

    # Seal numbers for age 25+
    seal.old <- ((seal[26]+seal[25]) * exp(-M/2) - remov[25] - remov[26]) * exp(-M/2)


    # Pups produced
    pup.prod4 <-  preg.rate[4, year+1] * (c(seal.age1, seal.num, seal.old)/2)[3]   # pup production for 3 years old individuals
    pup.prod5 <-  preg.rate[5, year+1] * (c(seal.age1, seal.num, seal.old)/2)[4]
    pup.prod6 <-  preg.rate[6, year+1] * (c(seal.age1, seal.num, seal.old)/2)[5]
    pup.prod7 <-  preg.rate[7, year+1] * (c(seal.age1, seal.num, seal.old)/2)[6]
    pup.prod8plus <-  sum(preg.rate[8:26, year+1] * (c(seal.age1, seal.num, seal.old)/2)[7:25])
    pup.not8plus <- sum(pup.prod4, pup.prod5, pup.prod6, pup.prod7)
    pup.prod.byage[,(year+1)] <- c(pup.prod4, pup.prod5, pup.prod6, pup.prod7, pup.prod8plus, pup.not8plus)

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

