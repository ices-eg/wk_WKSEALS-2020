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

#outdir <- paste(rep, "/modelAMKcor_output/", sep="") # set name of the output folder
#  if (!file.exists(substring(outdir,1, nchar(outdir)-1))) {dir.create(outdir)}  #create directory if it does not exist




#########################################################################
#########################################################################
#########################################################################
#########################################################################
####                PROJECTION FUNCTION                        ##########
#########################################################################
#########################################################################
#########################################################################
#########################################################################


proj.sealmod.sim.f <- function(result,
                               last.data.year,
                               end.year.sim,
                               nb.proj,
                               quota,
                               preg.cor = 0.85,
                               fish = F){
  # The result of the function is the pup production and total population
  # quota is a vector of quota which should be at least equal or longer than the number of projected years.
  # You have to chose your correlation in pregnancy rates within years (the default = 0.85).
  # Write fish = T if you want to save the projected population matrix



  ## Load the parameters
  # preg.rate.se <- sealmod.params$preg.rate.se
  preg.rate <- sealmod.params$preg.rate
  preg.size <- sealmod.params$preg.size
  pup.prod <- sealmod.params$pup.prod
  pup.prod.se <- sealmod.params$pup.prod.se

  ## Create the objects in which results will be stored
  final.sim.pup <- matrix(NA, nb.proj, end.year.sim - 1952 + 1)
  final.sim.total <- matrix(NA, nb.proj, end.year.sim - 1952 + 1)

  if(fish == T) {
    save.fish <- matrix(NA, 26, nb.proj)     # Initialize the matrix to save the population mean for the fish guy
    save.fish.sd <- matrix(NA, 26, nb.proj)  # Initialize the matrix to save the population standard deviation for the fish guy
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

  pb1 <- txtProgressBar(min = 0, max = nb.proj, style = 3) #Progress bar

  for (j in 1:nb.proj){

    setTxtProgressBar(pb1, j) # update progress bar

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
    proj.seal.sim <- matrix(ncol = year.proj + 1, nrow = 26)
    proj.seal.sim[,1] <- init.sim

    ## NEW REMOVAL VECTOR (includes all the catch from canada, greenland, arctic and bycatch)
    prop.pup.killed <- c(0.95, 0.14, 0.03, 0.6)     #In order: canada, greenland, Arctic, bycatch,   we could easily add variance if needed, canada average of last 10 years
    SandLpup <- c(1/0.95, 1/0.5, 1/0.5, 1)          #Vector of S&L for pup, values are for : (in order) canada, greenland, arctic, bycatch.  I could add variance on that as well
    SandLadult <- c(1/0.5, 1/0.5, 1/0.5, 1)         #Vector of S&L for adult, values are for : (in order) canada, greenland, arctic, bycatch.  I could add variance on that as well


    ## LOOP BY YEAR

    for (i in 1:year.proj) {           #Note that the loop goes 1 after the last year simulated

      ## GREENLAND REMOVALS
      # Greenland CASE 1 - Part B
      green.remov <-  runif(1, 66000, 92000)

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
      arctic.remov <- runif(1, 999, 1001)
      bycatch.remov <- runif(1, 12289, 12291)

      ## TOTAL REMOVALS
      #minister<-c(1.0,1.0,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.1)
      minister <- 1
      proj.remov <- rbind(quota[i] * (sample(minister, 1, replace = TRUE)), green.remov,arctic.remov, bycatch.remov)      #Minister factor:if a uniform distribution use(runif(1,1.0,1.1))instead of sample
      prop.age.class <- prop.table(proj.seal.sim[2:26, i])				#This creates a table of age class proportion to evaluate the kills

      # The next 3 lines create a vector of number of seal killed per age class which will be used as the remove vector in the projection model
      adult.kill <- sum((1 - prop.pup.killed) * proj.remov * SandLadult) * prop.age.class
      pup.kill <- sum(proj.remov * prop.pup.killed * SandLpup)
      kill <- append(pup.kill, adult.kill)

      ## EVENTS INDUCING MORTALITY OR ACTING ON PREGNANCY RATES
      # Ice factor
      #         icefactor <- 1 # uncomment this one and comment the following if you want to ignore the ice factor
      icefactor <- sample(c(1, 1, .85, .86, .88, 0.714466707, 0.35109587, 0.76642337, 0.772912202), 1, replace = TRUE) # define the proportion of pup surviving due to ice condition   new quantitative measure
      #         icefactor <- sample(c(.94,0.59,0.21,0.9,1),1) # define the proportion of pup surviving due to ice condition  -old measure

      # Food factor
      #         test <- c(1.76,1.80,0.2,0.45,1.92,0.19,0.01,0.27,0.56,1.47,0.24,0.93,1.32,0.65,0.76,1.14,1.19,0.67,1.32,
      #         0.05,1.45,1.02,1.38,1.43,1.43,0.12,0.72,1.67,1.71) # create the vector


      foodfactor <- 1 # uncomment this one and comment the following if you want to ignore the ice factor
      #        foodfactor <- sample(test), 1) # define good or bad year for food.
      # Act in reducing the effective pregnancy rate in a bad year. Here at mean 1 over 5 years is bad (10% of the normal pregnancy rates).
      # The best year (1 over 5) induces a 10% increase over normal pregnancy rates (110% of the normal pregnancy rates)

      #           foodfactor <- sample(test)
      ## APPLY MORTALITY
      # Seal numbers at age 1
      proj.seal.age1 <- (proj.seal[1] * icefactor - kill[1]) * exp(-gamma.pup * M)  *(1-(sum(proj.seal.sim[,(i)])/K)^2.4) # Include icefactor (see "Ice factor" section) and density dependence
      proj.seal.age1[is.nan(proj.seal.age1)] <- 0.001
      proj.seal.age1[proj.seal.age1 < 0] <- 0.001

      # Seal numbers for age > 1 and < 25
      proj.seal.num <-  (proj.seal[2:24]*exp(-M/2) - kill[2:24]) * exp(-M/2)
      proj.seal.num[proj.seal.num < 0] <- 0.001

      # Seal numbers for age 25+
      #proj.seal.old <- (proj.seal[25] * exp(-M/2) - kill[25]) * exp(-M/2)
      proj.seal.old <- ((proj.seal[26] + proj.seal[25]) * exp(-M/2) - kill[25] - kill[26]) * exp(-M/2)
      proj.seal.old[proj.seal.old < 0] <- 0.001

      ## PREGNANCY RATES OPTIONS
      # OPTION 1
      #         proj.preg.rate <- preg.rate[,length(preg.rate[1,])]  			# Is the last vector of pregnancy rates, therefore the pregnancy rate for future years is believed stable and similar to the one from the last data year.  this could be changed easily.
      #         proj.preg.size <- preg.size[,length(preg.size[1,])]

      # OPTION 2:
      #         proj.preg.rate <- c(0,0,0,  0.05,0.1,0.2,0.2,   rep(0.7,19))         #rep(0.7,19) means that 0.7 is repeated 19 times. It represents the pregnancy rates for ages 8 to 26 (i.e., 8+).
      #         proj.preg.size <- c(1,1,1,  5,5,5,5, rep(40,19))                     # sample size 40 is repeated 19 times for the 8+ these values can be changed

      # OPTION 3:
      #         proj.preg.rate <- c(0,0,0,  0.05,0.1,0.2,0.3,   rep(runif(1,.2,.7),19))         #rep(0.7,19) means that 0.7 is repeated 19 times. It represents the pregnancy rates for ages 8 to 26 (i.e., 8+).
      proj.preg.rate <- c(0, 0, 0,  # age 1 to 3
                          sample(c(.02, .02, .02, .02, .02, .02, .02, .03, .03, .03), 1, replace = TRUE),  #age 4
                          sample(c(.14, .13, .12, .11, .11, .10, .09, .08, .08, .07), 1, replace = TRUE),  #age 5
                          sample(c(.24, .23, .22, .22, .21, .20, .20, .19, .18, .17), 1, replace = TRUE),  # age 6  last 10 years sample(c(0.12,0.13,0.14,0.16,.18,.19,.21,.22,.24,.26),1),
                          sample(c(.39, .37, .36, .35, .33, .32, .31, .30, .29, .27), 1, replace = TRUE),  # age 7  last 10 years sample(c(0.33,0.34,.35,0.35,.36,.37,.38,.39,.64,.41),1),
                          rep(sample(c(.55, .65, .41, .64, .56, .76, .74, .55, .29, .20), 1, replace = TRUE), 19))    #   age 8+   last 10 years

      proj.preg.size <- c(1, 1, 1, 3, 4, 4, 4, rep(5, 19))  # sample size 5 is repeated 19 times for the 8+ these values can be changed
      rand.proj.preg.rate <- rbinom(26, size = proj.preg.size, prob = proj.preg.rate)/proj.preg.size

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
      #         rand.proj.preg.rate <- rep(0,26)

      #         # K effect on preg rates of 8+ classes
      #         rand.proj.preg.rate[8:26] <- rep(ifelse(0.88 * (1 - (sum(proj.seal.sim[,i])/K)^2.4) < 0, 0, (0.88 * (1 - (sum(proj.seal.sim[,i])/K)^2.4))),19)
      # 0.88 is considered as max preg rate (can be defined as max preg rate when pop is low)

      #          # This relation was considered as a logistic function fitted on raw pregrates for classes 4-8
      #          pregraw <- read.csv("temppr_2014.csv", header=F); colnames(pregraw) <- c("years", "age", "total", "npreg")
      #          pregraw$pregrate <- pregraw$npreg / pregraw$total
      #          pregraw <- na.omit(pregraw)
      #          logispreg <- nls(pregrate~SSlogis(age,a,m,s), weights = total, data=pregraw)   # weighted by sample size
      #
      #          # Sub-option 1
      #          rand.proj.preg.rate[4:7]<-SSlogis(4:7, rand.proj.preg.rate[8], coef(logispreg)[2], coef(logispreg)[3])  # the first parameter (the asymptotic value) is replaced by the preg.rate of classe 8 in order to estimate values for classes 4-7
      #          # Sub-option 2
      #           rand.proj.preg.rate[4:7]<-SSlogis(4:7, rand.proj.preg.rate[8], 5.368324, 0.5674322) # application of the previous line but avoiding recalculation of the logistic curve for each run (Warning, coefficient values could have changed with new data)
      #          # Sub-option 3
      #           rand.proj.preg.rate[4:7]<-SSlogis(4:7, rand.proj.preg.rate[8], rnorm(1,coef(logispreg)[2], summary(logispreg)$coefficient[2,2]), rnorm(1,coef(logispreg)[3], summary(logispreg)$coefficient[3,2]))  # same as previous line but allows variability for paramaters 2 and 3
      #          # Sub-option 4
      #           rand.proj.preg.rate[4:7]<-SSlogis(4:7, rand.proj.preg.rate[8], rnorm(1,mean = 5.368324, 0.13783), rnorm(1,mean = 0.5674322, 0.12531)) # application of the previous line but avoiding recalculation of the logistic curve for each run turn this on if running density dependence preg rates (Warning, coefficient values could have changed with new data)

      ## END PREGNANCY RATES OPTIONS


      ## APPLY NATALITY
      # Apply effect on natural event on pregnancy rates
      rand.proj.preg.rate <- ifelse(foodfactor * rand.proj.preg.rate > 0.88, 0.88, foodfactor * rand.proj.preg.rate) #(see section "Food factor" to check which kind of food factor is used)

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

    if (fish == T) {
      data.fish <- matrix(NA, 26, 2, dimnames = list(c(1:26), c("mean", "Stand Dev")))  ##########  Code to save mean and SD of seal population for the projected year for the fish guy
      save.fish[,j] <- apply(proj.seal.sim[,2:(year.proj + 1)], 1, mean)                ##  To activate, write "fish=T" while initializing the function
      save.fish.sd[,j] <- apply(proj.seal.sim[,2:(year.proj + 1)], 1, sd)               ##  The data are save in data.fish
      data.fish[,1] <- apply(save.fish, 1, mean)                                        ##
      data.fish[,2] <- apply(save.fish.sd, 1, mean)                                     ##
      data.fish <<- data.fish                                                           ##
    }                                                                                   ##########  end of the fish code

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

  test.table <- cbind(final.sim.pup.0.5, final.sim.pup.mean,final.sim.pup.sd, final.sim.pup.0.025,final.sim.pup.0.975,
                      final.sim.total.0.5,final.sim.total.mean, final.sim.total.sd, final.sim.total.0.025, final.sim.total.05,final.sim.total.10, final.sim.total.20,
                      final.sim.total.30, final.sim.total.40, final.sim.total.95, final.sim.total.0.975)
  write.csv(test.table, "output/test_reference.csv")

    ## Plot the pup distribution,
  dev.new()
  par(mfrow = c(2,1), mar = c(4.5,5,1,1))

  plot(1952:end.year.sim, pup.prod[1:(end.year.sim-1952+1)]/1e6, ylim=c(0,signif(max(final.sim.pup.0.975),2)/1e6),xlab="Years",ylab="Pup prod estimates (x 1e+06)", cex.lab = 1.2, cex.axis = 1.2, pch = 19)
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

  plot(1952:end.year.sim, final.sim.total.0.5, ylim=c(0,signif(max(final.sim.total.0.975)/1e6,2)),xlab="Years",ylab="Seal numbers (x 1e+06)", cex.lab = 1.2, cex.axis = 1.2, type = "l")
  polygon(x = c(1952:end.year.sim, rev(1952:end.year.sim)), y = c(final.sim.total.0.025/1e6,rev(final.sim.total.0.975)/1e6), col = "darkgrey", border = NA)
  abline(h = seq(0,signif(max(final.sim.total)/1e6,2),by = 1), col = "gray80") #plots the grey horizontal lines
  lines(1952:end.year.sim,final.sim.total.0.5/1e6)

  par(mfrow = c(1,1))
  dev.copy2pdf(file = paste("output/Projection_output_for_", nb.proj, "_iterations.pdf", sep = ""))

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

